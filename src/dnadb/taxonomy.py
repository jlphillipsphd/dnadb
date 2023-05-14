from dataclasses import dataclass, field, replace
from functools import cached_property
from itertools import chain
import io
import json
from lmdbm import Lmdb
import numpy as np
import numpy.typing as npt
from pathlib import Path
import re
from sortedcontainers import SortedDict
from typing import Dict, Generator, Iterable, List, Optional, overload, Tuple, Union

from .db import DbFactory, DbWrapper
from .utils import open_file

RANKS = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
RANK_PREFIXES = ''.join(rank[0] for rank in RANKS).lower()

# Utility Functions --------------------------------------------------------------------------------

def split_taxonomy(taxonomy: str) -> Tuple[str, ...]:
    """
    Split taxonomy label into a tuple
    """
    # print(f'Split: {tuple(re.findall(r"\w__([^;]+)", taxonomy))}')
    return tuple(re.findall(r"\w__([^;]+)", taxonomy))


def join_taxonomy(taxonomy: Union[Tuple[str, ...], List[str]]) -> str:
    """
    Merge a taxonomy tuple into a string format
    """
    depth = len(taxonomy)
    assert depth >= 1 and depth <= len(RANKS), "Invalid taxonomy"
    taxonomy = tuple(taxonomy) + ("",) * (len(RANKS) - depth)
    return "; ".join([f"{RANK_PREFIXES[i]}__{taxon}" for i, taxon in enumerate(taxonomy)])

# Taxonomy TSV Utilities ---------------------------------------------------------------------------

@dataclass(frozen=True, order=True)
class TaxonomyEntry:
    __slots__ = ("identifier", "label")
    identifier: str
    label: str

    @classmethod
    def deserialize(cls, entry: bytes) -> "TaxonomyEntry":
        return cls.from_str(entry.decode())

    @classmethod
    def from_str(cls, entry: str) -> "TaxonomyEntry":
        """
        Create a taxonomy entry from a string
        """
        identifier, taxonomy = entry.rstrip().split('\t')
        return cls(identifier, taxonomy)

    def taxons(self) -> Tuple[str, ...]:
        return split_taxonomy(self.label)

    def serialize(self) -> bytes:
        return str(self).encode()

    def __str__(self):
        return f"{self.identifier}\t{self.label}"


class TaxonomyDbFactory(DbFactory):
    """
    A factory for creating LMDB-backed databases of taxonomy entries.

    [index to label]
    0 -> k__bacteria;...
    1 -> ...
    ...

    [label to index]
    k__bacteria;... -> 0
    ... -> 1
    ...

    [label counts]
    0_count -> 2
    1 -> 1
    ...

    [label index to fasta id]
    0_0 -> abc
    0_1 -> def
    1_0 -> efg
    ...

    [fasta_id to label index]
    abc -> 0
    def -> 0
    efg -> 1
    ...
    """
    __slots__ = ("num_entries",)

    def __init__(self, path: Union[str, Path], chunk_size: int = 10000):
        super().__init__(path, chunk_size)
        self.num_entries = np.int32(0)

    def write_entry(self, entry: TaxonomyEntry):
        """
        Create a new taxonomy LMDB database from taxonomy entries.
        """
        if not self.contains(entry.label):
            # index -> label, label -> index
            self.write(str(self.num_entries), entry.label.encode())
            self.write(entry.label, self.num_entries.tobytes())
            self.write(f"count_{self.num_entries}", np.int32(0).tobytes())
            self.num_entries += 1
        index: np.int32 = np.frombuffer(self.read(entry.label), dtype=np.int32, count=1)[0]
        count: np.int32 = np.frombuffer(self.read(f"count_{index}"), dtype=np.int32, count=1)[0]
        self.write(f"{index}_{count}", entry.identifier.encode())
        self.write(f">{entry.identifier}", index.tobytes())
        self.write(f"count_{index}", (count + 1).tobytes())

    def write_entries(self, entries: Iterable[TaxonomyEntry]):
        for entry in entries:
            self.write_entry(entry)

    def before_close(self):
        self.write("length", self.num_entries.tobytes())
        super().before_close()


class TaxonomyDb(DbWrapper):
    __slots__ = ("length",)

    def __init__(self, taxonomy_db_path: Union[str, Path]):
        super().__init__(taxonomy_db_path)
        self.length = np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]

    def contains_id(self, fasta_identifier: str) -> bool:
        """
        Check if a FASTA identifier exists in the database.
        """
        return f">{fasta_identifier}" in self.db

    def contains_label(self, label: str) -> bool:
        """
        Check if a taxonomy label exists in the database.
        """
        return label in self.db

    def count(self, label_index: int) -> int:
        """
        Get the number of sequences with a given label index.
        """
        return int(np.frombuffer(self.db[f"count_{label_index}"], dtype=np.int32, count=1)[0])

    def counts(self) -> Generator[int, None, None]:
        """
        Get the number of sequences for each label index.
        """
        for i in range(self.length):
            yield self.count(i)

    def id_to_index(self, fasta_identifier: str) -> int:
        """
        Get the taxonomy index for a given FASTA identifier.
        """
        return int(np.frombuffer(self.db[f">{fasta_identifier}"], dtype=np.int32, count=1)[0])

    def id_to_label(self, fasta_identifier: str) -> str:
        """
        Get the taxonomy label for a given FASTA identifier.
        """
        return self.label(self.id_to_index(fasta_identifier))

    def label(self, index: int) -> str:
        """
        Get the taxonomy label for a given index.
        """
        return self.db[str(index)].decode()

    def labels(self) -> Generator[str, None, None]:
        """
        Get the taxonomy labels.
        """
        for i in range(self.length):
            yield self.label(i)

    def label_to_index(self, label: str) -> np.int32:
        """
        Get the taxonomy index for a given label.
        """
        return np.frombuffer(self.db[label], dtype=np.int32, count=1)[0]

    # @cached_property
    # def hierarchy(self):
    #     return TaxonomyHierarchy.deserialize(self.db["hierarchy"])

    def __len__(self):
        return self.length

    def __iter__(self):
        for i in range(len(self)):
            yield self.db[str(i)].decode()

# Taxonomy Hierarchy -------------------------------------------------------------------------------

@dataclass(frozen=True, order=True)
class Taxon:
    # __slots__ = ("name", "rank", "parent", "children")
    name: str
    rank: int
    parent: Optional["Taxon"] = field(default=None, hash=False)
    children: Dict[str, "Taxon"] = field(default_factory=dict, init=False, hash=False, repr=False)

    def __post_init__(self):
        if self.parent is not None:
            self.parent.children[self.name.casefold()] = self

    def __contains__(self, other: Union[str, "Taxon"]):
        if isinstance(other, Taxon):
            other = other.name
        return other.casefold() in self.children

    def __iter__(self):
        yield from self.children.values()

    def __eq__(self, other: Union[str, "Taxon"]):
        if isinstance(other, str):
            return self.name.casefold() == other.casefold()
        return self.rank == other.rank and self.name.casefold() == other.name.casefold()

    def __getitem__(self, key: Union[str, "Taxon"]):
        if isinstance(key, Taxon):
            key = key.name
        return self.children[key.casefold()]

class TaxonomyHierarchy:

    class TaxonDict(SortedDict):
        """
        A dictionary of Taxon instances with alphabetically-sorted, case-insensitive keys.
        """
        def insert(self, taxon: Taxon):
            """
            Insert a Taxon instance into the dictionary using the Taxon's name as the key.
            """
            self[taxon.name] = taxon

        def __contains__(self, key: Union[str, Taxon]) -> bool:
            if isinstance(key, Taxon):
                key = key.name
            return super().__contains__(key.casefold())

        def __getitem__(self, key: Union[str, Taxon]) -> Taxon:
            if isinstance(key, Taxon):
                key = key.name
            return super().__getitem__(key.casefold())

        def __setitem__(self, key: Union[str, Taxon], value: Taxon):
            if isinstance(key, Taxon):
                key = key.name
            super().__setitem__(key.casefold(), value)

        def __delitem__(self, key: Union[str, Taxon]):
            if isinstance(key, Taxon):
                key = key.name
            super().__delitem__(key.casefold())

        def keys(self) -> chain[str]:
            return super().keys()

        def values(self) -> chain[Taxon]:
            return super().values()

        def __iter__(self) -> chain[str]:
            return super().__iter__()

    @classmethod
    def merged(
        cls,
        hierarchies: Iterable["TaxonomyHierarchy"],
        depth: Optional[int] = None
    ) -> "TaxonomyHierarchy":
        """
        Merge multiple taxonomy hierarchies into one.
        """
        depth = depth if depth is not None else max(hierarchy.depth for hierarchy in hierarchies)
        merged = cls(depth)
        for hierarchy in hierarchies:
            for taxon in hierarchy:
                if taxon.rank >= depth:
                    break
                merged.taxons[taxon.rank][taxon.name] = Taxon(
                    taxon.name,
                    taxon.rank,
                    parent=(merged.taxons[taxon.rank - 1][taxon.parent.name]
                            if taxon.parent else None))
        return merged

    def __init__(self, depth: int):
        self.depth = depth
        self.taxons = tuple(TaxonomyHierarchy.TaxonDict() for _ in range(depth))
        self.__taxon_to_id_map: Optional[Tuple[Dict[Taxon, int], ...]] = None
        self.__id_to_taxon_map: Optional[Tuple[Dict[int, Taxon], ...]] = None

    def add_entries(self, entries: Iterable[TaxonomyEntry]):
        """
        Add taxonomy entries to the hierarchy.

        Args:
            entries (Iterable[TaxonomyEntry]): The taxonomy entries to add.
        """
        for entry in entries:
            self.add_entry(entry)

    def add_entry(self, entry: TaxonomyEntry):
        """
        Add a taxonomy entry to the hierarchy.

        Args:
            entry (TaxonomyEntry): The taxonomy entry to add.
        """
        self.add_taxonomy(entry.label)

    def add_taxonomies(self, taxonomies: Iterable[str]):
        """
        Add taxonomy labels to the hierarchy.

        Args:
            taxonomies (Iterable[str]): The taxonomy labels to add (e.g. "k__Bacteria; ...").
        """
        for taxonomy in taxonomies:
            self.add_taxonomy(taxonomy)

    def add_taxonomy(self, taxonomy: str):
        """
        Add a taxonomy to the hierarchy.

        Args:
            taxonomy (str): The taxonomy label to add (e.g. "k__Bacteria; ...").
        """
        return self.add_taxons(split_taxonomy(taxonomy))

    def add_taxons(self, taxonomy: Tuple[str, ...]):
        """
        Add a taxonomy in the form of a taxon tuple to the hierarchy.

        Args:
            taxonomy (Tuple[str, ...]): The taxon tuple to add.
        """
        if isinstance(taxonomy, str):
            taxonomy = split_taxonomy(taxonomy)
        parent: Optional[Taxon] = None
        for rank, name in enumerate(taxonomy):
            if not self.has_taxon(name, rank):
                self.__taxon_to_id_map = None
                self.__id_to_taxon_map = None
                self.taxons[rank].insert(Taxon(name, rank, parent))
            parent = self.taxons[rank][name]

    def has_entry(self, taxonomy: TaxonomyEntry, strict: bool = False) -> bool:
        """
        Check if the hierarchy has the given taxonomy.

        Args:
            taxonomy (str): The taxonomy label to check (e.g. "k__Bacteria; ...").
            strict (bool, optional): Ensure every taxon exists. Defaults to False.

        Returns:
            bool: The hierarchy contains the taxonomy.
        """
        return self.has_taxonomy(taxonomy.label, strict)

    def has_taxonomy(self, taxonomy: str, strict: bool = False) -> bool:
        """
        Check if the hierarchy has the given taxonomy.

        Args:
            taxonomy (TaxonomyEntry): The taxonomy to check.
            strict (bool, optional): Ensure every taxon exists. Defaults to False.

        Returns:
            bool: The hierarchy contains the taxonomy.
        """
        if isinstance(taxonomy, TaxonomyEntry):
            taxonomy = taxonomy.label
        taxons = split_taxonomy(taxonomy)
        return self.has_taxons(taxons, strict)

    def has_taxon(self, taxon: str, rank: int) -> bool:
        """
        Check if the given taxon is present within the hierarchy.

        Args:
            taxon (str): The taxon name to check.
            rank (int): The rank of the taxon.

        Returns:
            bool: The taxon is present within the hierarchy.
        """
        return (rank < self.depth) and taxon.casefold() in self.taxons[rank]

    def has_taxons(self, taxons: Tuple[str, ...], strict: bool = False) -> bool:
        """
        Check if the hierarchy has a list of taxons.

        Note: When strict == False, a valid hierarchy is assumed.
        """
        if not strict:
            return self.has_taxon(taxons[-1], len(taxons) - 1)
        return all(taxon.casefold() in self.taxons[rank] for rank, taxon in enumerate(taxons))

    def reduce_entry(self, taxonomy: TaxonomyEntry, strict: bool = False) -> TaxonomyEntry:
        """
        Reduce the taxonomy to a valid known taxonomy label in the hierarchy.

        Args:
            taxonomy (TaxonomyEntry): The taxonomy entry to reduce.
            strict (bool, optional): Ensure every taxon exists. Defaults to False.

        Returns:
            TaxonomyEntry: A new TaxonomyEntry instance containing the reduced taxonomy label.
        """
        return replace(taxonomy, label=self.reduce_taxonomy(taxonomy.label, strict))

    def reduce_taxonomy(self, taxonomy: str, strict: bool = False) -> str:
        """
        Reduce the taxonomy to a valid known taxonomy label in the hierarchy.

        Args:
            taxonomy (str): The taxonomy label to reduce (e.g. "k__Bacteria; ...").
            strict (bool, optional): Ensure every taxon exists. Defaults to False.

        Returns:
            str: The reduced taxonomy label.
        """
        return join_taxonomy(self.reduce_taxons(split_taxonomy(taxonomy), strict))

    def reduce_taxons(self, taxons: Tuple[str, ...], strict: bool = False) -> Tuple[str, ...]:
        """
        Reduce the taxonomy tuple to a valid known taxonomy in the hierarchy.

        Args:
            taxons (Tuple[str, ...]): The taxonomy tuple to reduce.
            strict (bool, optional): Ensure every taxon exists. Defaults to False.

        Returns:
            Tuple[str, ...]: The reduced taxonomy tuple.
        """
        if strict:
            for rank, taxon in enumerate(taxons):
                if taxon not in self.taxons[rank]:
                    return taxons[:rank]
            return taxons
        for i in range(len(taxons) - 1, -1, -1):
            if taxons[i] in self.taxons[i]:
                return taxons[:i + 1]
        return tuple()

    def tokenize(self, taxonomy: str, pad: bool = False) -> npt.NDArray[np.int64]:
        """
        Tokenize the taxonomy label into a a tuple of taxon integer IDs

        Args:
            taxonomy (str): The taxonomy label to tokenize (e.g. "k__Bacteria; ...").
            pad (bool): Pad the taxonomy with -1s to the depth of the hierarchy. Defaults to False.

        Returns:
            Tuple[str, ...]: The tokenized taxonomy.
        """
        return self.tokenize_taxons(split_taxonomy(taxonomy), pad)

    def tokenize_taxons(self, taxons: Tuple[str, ...], pad: bool = False) -> npt.NDArray[np.int64]:
        """
        Tokenize the taxonomy tuple into a a tuple of taxon integer IDs

        Args:
            taxons (Tuple[str, ...]): The taxonomy tuple to tokenize.
            pad (bool): Pad the taxonomy with -1s to the depth of the hierarchy. Defaults to False.

        Returns:
            Tuple[str, ...]: The tokenized taxonomy.
        """
        result = np.empty(len(taxons), np.int64) if not pad else np.full(self.depth, -1, np.int64)
        for rank, taxon in enumerate(taxons):
            result[rank] = self.taxon_to_id_map[rank][self.taxons[rank][taxon]]
        return result

    def detokenize(self, taxon_tokens: npt.NDArray[np.int64]) -> str:
        """
        Detokenize the taxonomy tokens into a taxonomy label.

        Args:
            taxon_tokens (npt.NDArray[np.int64]): The taxonomy tokens.

        Returns:
            str: The detokenized taxonomy label.
        """
        return join_taxonomy(self.detokenize_taxons(taxon_tokens))

    def detokenize_taxons(self, taxon_tokens: npt.NDArray[np.int64]) -> Tuple[str, ...]:
        """
        Detokenize the taxonomy tokens into a taxonomy tuple.

        Args:
            taxon_tokens (npt.NDArray[np.int64]): The taxonomy tokens.

        Returns:
            Tuple[str, ...]: The detokenized taxonomy tuple.
        """
        result = tuple()
        for rank, token in enumerate(taxon_tokens):
            if token < 0:
                break
            result += (self.id_to_taxon_map[rank][token].name,)
        return result

    @property
    def taxon_to_id_map(self) -> Tuple[Dict[Taxon, int], ...]:
        """
        A mapping of taxon instances to taxon IDs.
        """
        if self.__taxon_to_id_map is None:
            self.__taxon_to_id_map = tuple(
                {taxon: i for i, taxon in enumerate(taxons.values())}
                for taxons in self.taxons)
        return self.__taxon_to_id_map

    @property
    def id_to_taxon_map(self) -> Tuple[Dict[int, Taxon], ...]:
        """
        A mapping of taxon IDs to Taxon instances.
        """
        if self.__id_to_taxon_map is None:
            self.__id_to_taxon_map = tuple(
                {i: taxon for i, taxon in enumerate(taxons.values())}
                for taxons in self.taxons)
        return self.__id_to_taxon_map

    def __iter__(self) -> Generator[Taxon, None, None]:
        """
        Breadth-first iteration over the taxonomy hierarchy.
        """
        for rank in range(self.depth):
            yield from self.taxons[rank].values()


def entries(
    taxonomy: Union[io.TextIOBase, Iterable[TaxonomyEntry], str, Path]
) -> Iterable[TaxonomyEntry]:
    """
    Create an iterator over a taxonomy file or iterable of taxonomy entries.
    """
    if isinstance(taxonomy, (str, Path)):
        with open_file(taxonomy, 'r') as buffer:
            yield from read(buffer)
    elif isinstance(taxonomy, io.TextIOBase):
        yield from read(taxonomy)
    else:
        yield from taxonomy


def read(buffer: io.TextIOBase) -> Generator[TaxonomyEntry, None, None]:
    """
    Read taxonomies from a tab-separated file (TSV)
    """
    for line in buffer:
        identifier, taxonomy = line.rstrip().split('\t')
        yield TaxonomyEntry(identifier, taxonomy)


def write(buffer: io.TextIOBase, entries: Iterable[TaxonomyEntry]):
    """
    Write taxonomy entries to a tab-separate file (TSV)
    """
    for entry in entries:
        buffer.write(f"{entry.identifier}\t{entry.label}\n")
