from dataclasses import dataclass, field, replace
import heapq
import io
import json
import numpy as np
import numpy.typing as npt
from pathlib import Path
import re
from typing import Dict, Generator, Iterable, List, Optional, Tuple, Union

from .db import DbFactory, DbWrapper
from .types import int_t
from .utils import open_file

RANKS = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
RANK_PREFIXES = ''.join(rank[0] for rank in RANKS).lower()

# Utility Functions --------------------------------------------------------------------------------

def split_taxonomy(taxonomy: str, keep_empty: bool = False) -> Tuple[str, ...]:
    """
    Split taxonomy label into a tuple
    """
    return tuple(re.findall(r"\w__([^;]*)" if keep_empty else r"\w__([^;]+)", taxonomy))


def join_taxonomy(taxonomy: Union[Tuple[str, ...], List[str]], depth: Optional[int] = None) -> str:
    """
    Merge a taxonomy tuple into a string format
    """
    if depth is None:
        depth = len(taxonomy)
    assert depth >= 1 and depth <= len(RANKS), "Invalid taxonomy"
    taxonomy = tuple(taxonomy) + ("",)*(depth - len(taxonomy))
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

    def contains_fasta_id(self, fasta_identifier: str) -> bool:
        """
        Check if a FASTA identifier exists in the database.
        """
        return f">{fasta_identifier}" in self.db

    def contains_label(self, label: str) -> bool:
        """
        Check if a taxonomy label exists in the database.
        """
        return label in self.db

    def count(self, label_index: int_t) -> int:
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

    def fasta_id_with_label(self, label_index: int_t, fasta_index: int_t) -> str:
        """
        Get the FASTA identifier for a given label and index.
        """
        return self.db[f"{label_index}_{fasta_index}"].decode()

    def fasta_ids_with_label(self, label_index: int_t) -> Generator[str, None, None]:
        """
        Get the FASTA identifiers for a given label.
        """
        for i in range(self.count(label_index)):
            yield self.fasta_id_with_label(label_index, i)

    def fasta_id_to_index(self, fasta_identifier: str) -> int:
        """
        Get the taxonomy index for a given FASTA identifier.
        """
        return int(np.frombuffer(self.db[f">{fasta_identifier}"], dtype=np.int32, count=1)[0])

    def fasta_id_to_label(self, fasta_identifier: str) -> str:
        """
        Get the taxonomy label for a given FASTA identifier.
        """
        return self.label(self.fasta_id_to_index(fasta_identifier))

    def label(self, label_index: int_t) -> str:
        """
        Get the taxonomy label for a given index.
        """
        return self.db[str(label_index)].decode()

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

    def __len__(self):
        return self.length

    def __iter__(self):
        for i in range(len(self)):
            yield self.db[str(i)].decode()

# Taxonomy Hierarchy -------------------------------------------------------------------------------

@dataclass(frozen=True, order=True)
class Taxon:
    # __slots__ = ("name", "rank", "parent", "children")
    rank: int
    name: str
    parent: Optional["Taxon"] = field(default=None, repr=False)
    children: Dict[str, "Taxon"] = field(default_factory=dict, init=False, hash=False, repr=False)

    def add_child(self, name: str):
        assert name not in self, f"Attempted to add duplicate child taxon: {name}"
        return Taxon(self.rank + 1, name, self)

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
        return self.rank == other.rank \
            and self.name.casefold() == other.name.casefold() \
            and self.parent == other.parent

    def __getitem__(self, key: Union[str, "Taxon"]):
        if isinstance(key, Taxon):
            key = key.name
        return self.children[key.casefold()]

class TaxonomyHierarchy:

    @classmethod
    def deserialize(cls, hierarchy_json_bytes: Union[str, bytes, bytearray]) -> "TaxonomyHierarchy":
        """
        Deserialize a taxonomy hierarchy from a bytes object.
        """
        def recursive_add(taxon_json: Dict, parent: Taxon):
            taxon = parent.add_child(taxon_json["name"])
            for child_json in taxon_json["children"]:
                recursive_add(child_json, taxon)
        hierarchy_json = json.loads(hierarchy_json_bytes)
        hierarchy = cls(hierarchy_json["depth"])
        for taxon_json in hierarchy_json["taxons"]:
            recursive_add(taxon_json, hierarchy.taxon_tree_head)
        return hierarchy

    @classmethod
    def from_dbs(cls, dbs: Iterable[TaxonomyDb], depth: int = 6) -> "TaxonomyHierarchy":
        """
        Create a taxonomy hierarchy from multiple taxonomy databases.
        """
        hierarchy = cls(depth)
        for db in dbs:
            hierarchy.add_taxonomies(db)
        return hierarchy

    @classmethod
    def merged(
        cls,
        hierarchies: Iterable["TaxonomyHierarchy"],
        depth: Optional[int] = None
    ) -> "TaxonomyHierarchy":
        """
        Merge multiple taxonomy hierarchies into one.
        """
        depth = depth if depth is not None else min(hierarchy.depth for hierarchy in hierarchies)
        def recursive_insert(head: Taxon, taxon: Taxon):
            if taxon.rank >= depth:
                return
            if taxon not in head:
                head.add_child(taxon.name)
            for child in taxon:
                recursive_insert(head[taxon.name], child)
        merged = cls(depth)
        for hierarchy in hierarchies:
            for child in hierarchy.taxon_tree_head:
                recursive_insert(merged.taxon_tree_head, child)
        return merged

    def __init__(self, depth: int = 6):
        self.depth = depth
        self.taxon_tree_head = Taxon(-1, "_root_")
        self._taxon_to_id_map: Optional[Tuple[Dict[Taxon, int], ...]] = None
        self._id_to_taxon_map: Optional[Tuple[List[Taxon], ...]] = None

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

    def add_taxons(self, taxons: Tuple[str, ...]):
        """
        Add a taxonomy in the form of a taxon tuple to the hierarchy.

        Args:
            taxonomy (Tuple[str, ...]): The taxon tuple to add.
        """
        head = self.taxon_tree_head
        for taxon in taxons[:self.depth]:
            if taxon not in head:
                self._id_to_taxon_map = None
                self._taxon_to_id_map = None
                head = head.add_child(taxon)
            else:
                head = head[taxon]

    def has_entry(self, entry: TaxonomyEntry) -> bool:
        """
        Check if the hierarchy has the given taxonomy.

        Args:
            entry (str): The taxonomy entry to check

        Returns:
            bool: The hierarchy contains the taxonomy.
        """
        return self.has_taxonomy(entry.label)

    def has_taxonomy(self, taxonomy: str) -> bool:
        """
        Check if the hierarchy has the given taxonomy.

        Args:
            taxonomy (TaxonomyEntry): The taxonomy to check.

        Returns:
            bool: The hierarchy contains the taxonomy.
        """
        return self.has_taxons(split_taxonomy(taxonomy))

    def has_taxons(self, taxons: Tuple[str, ...]) -> bool:
        """
        Determine if the given taxons are present in the hierarchy.

        Args:
            taxons (Tuple[str, ...]): The taxon hierarchy to check.

        Returns:
            bool: The presence of the taxons in the hierarchy.
        """
        if len(taxons) > self.depth and taxons[self.depth] != "":
            return False
        head = self.taxon_tree_head
        for taxon in taxons[:self.depth]:
            if taxon == "":
                return True
            if taxon not in head:
                return False
            head = head[taxon]
        return True

    def reduce_entry(self, entry: TaxonomyEntry) -> TaxonomyEntry:
        """
        Reduce the taxonomy to a valid known taxonomy label in the hierarchy.

        Args:
            entry (TaxonomyEntry): The taxonomy entry to reduce.

        Returns:
            TaxonomyEntry: A new TaxonomyEntry instance containing the reduced taxonomy label.
        """
        return replace(entry, label=self.reduce_taxonomy(entry.label))

    def reduce_taxonomy(self, taxonomy: str) -> str:
        """
        Reduce the taxonomy to a valid known taxonomy label in the hierarchy.

        Args:
            taxonomy (str): The taxonomy label to reduce (e.g. "k__Bacteria; ...").

        Returns:
            str: The reduced taxonomy label.
        """
        return join_taxonomy(self.reduce_taxons(split_taxonomy(taxonomy)), depth=self.depth)

    def reduce_taxons(self, taxons: Tuple[str, ...]) -> Tuple[str, ...]:
        """
        Reduce the taxonomy tuple to a valid known taxonomy in the hierarchy.

        Args:
            taxons (Tuple[str, ...]): The taxonomy tuple to reduce.

        Returns:
            Tuple[str, ...]: The reduced taxonomy tuple.
        """
        result: Tuple[str, ...] = tuple()
        head = self.taxon_tree_head
        for taxon in taxons[:self.depth]:
            if taxon not in head:
                break
            result += (taxon,)
            head = head[taxon]
        return result

    def tokenize(
        self,
        taxonomy: str,
        pad: bool = False,
        include_missing: bool = False
    ) -> npt.NDArray[np.int32]:
        """
        Tokenize the taxonomy label into a a tuple of taxon integer IDs

        Args:
            taxonomy (str): The taxonomy label to tokenize (e.g. "k__Bacteria; ...").
            pad (bool): Pad the taxonomy with -1s to the depth of the hierarchy. Defaults to False.
            include_missing (bool): Assign missing taxons in the tokenized taxonomy to 0. Defaults to False.

        Returns:
            np.ndarray[np.int32]: The tokenized taxonomy.
        """
        return self.tokenize_taxons(split_taxonomy(taxonomy), pad, include_missing)

    def tokenize_taxons(
        self,
        taxons: Tuple[str, ...],
        pad: bool = False,
        include_missing: bool = False
    ) -> npt.NDArray[np.int32]:
        """
        Tokenize the taxonomy tuple into a a tuple of taxon integer IDs

        Args:
            taxons (Tuple[str, ...]): The taxonomy tuple to tokenize.
            pad (bool): Pad the taxonomy with -1s to the depth of the hierarchy. Defaults to False.
            include_missing (bool): Assign missing taxons in the tokenized taxonomy to 0. Defaults to False.

        Returns:
            np.ndarray[np.int32]: The tokenized taxonomy.
        """
        taxons = taxons[:self.depth] # todo should we trim silently or throw error?
        result = np.empty(len(taxons), np.int32) if not pad else np.full(self.depth, -1, np.int32)
        head = self.taxon_tree_head
        for taxon in taxons:
            head = head[taxon]
            result[head.rank] = self.taxon_to_id_map[head.rank][head]
        if include_missing:
            result += 1
        return result

    def detokenize(self, taxon_tokens: npt.NDArray[np.int32], include_missing: bool = False) -> str:
        """
        Detokenize the taxonomy tokens into a taxonomy label.

        Args:
            taxon_tokens (npt.NDArray[np.int64]): The taxonomy tokens.
            include_missing (bool): Assign missing taxons in the tokenized taxonomy to 0. Defaults to False.

        Returns:
            str: The detokenized taxonomy label.
        """
        return join_taxonomy(self.detokenize_taxons(taxon_tokens, include_missing), depth=self.depth)

    def detokenize_taxons(
        self,
        taxon_tokens: npt.NDArray[np.int32],
        include_missing: bool = False
    ) -> Tuple[str, ...]:
        """
        Detokenize the taxonomy tokens into a taxonomy tuple.

        Args:
            taxon_tokens (npt.NDArray[np.int64]): The taxonomy tokens.
            include_missing (bool): Assign missing taxons in the tokenized taxonomy to 0. Defaults to False.

        Returns:
            Tuple[str, ...]: The detokenized taxonomy tuple.
        """
        if include_missing:
            taxon_tokens -= 1
        result: Tuple[str, ...] = tuple()
        token: np.int32
        for rank, token in enumerate(taxon_tokens):
            if token < 0:
                break
            result += (self.id_to_taxon_map[rank][token].name,)
        return result

    def serialize(self) -> bytes:
        """
        Serialize the taxonomy hierarchy as a JSON string.
        """
        def dfs_serialize(head: Taxon):
            return {
                "name": head.name,
                "children": [dfs_serialize(child) for child in head]
            }
        taxons = []
        for head in self.taxon_tree_head:
            taxons.append(dfs_serialize(head))
        return json.dumps({
            "depth": self.depth,
            "taxons": taxons
        }).encode()

    @property
    def taxons(self) -> Tuple[Tuple[Taxon, ...], ...]:
        """
        A tuple of tuples of taxons at each rank in the hierarchy.
        """
        return tuple(tuple(taxons.keys()) for taxons in self.taxon_to_id_map)

    @property
    def taxon_counts(self) -> Tuple[int, ...]:
        """
        The number of taxons at each rank in the hierarchy.
        """
        return tuple(len(taxons) for taxons in self.taxon_to_id_map)

    @property
    def taxon_to_id_map(self) -> Tuple[Dict[Taxon, int], ...]:
        """
        A mapping of taxon instances to taxon IDs.
        """
        if self._taxon_to_id_map is None:
            self._taxon_to_id_map = tuple({} for _ in range(self.depth))
            for taxon in self:
                self._taxon_to_id_map[taxon.rank][taxon] = len(self._taxon_to_id_map[taxon.rank])
        return self._taxon_to_id_map

    @property
    def id_to_taxon_map(self) -> Tuple[List[Taxon], ...]:
        """
        A mapping of taxon IDs to Taxon instances.
        """
        if self._id_to_taxon_map is None:
            self._id_to_taxon_map = tuple([] for _ in range(self.depth))
            for taxon in self:
                # taxon_id = self.taxon_to_id_map[taxon.rank][taxon]
                # assert len(self._id_to_taxon_map[taxon.rank]) == taxon_id
                self._id_to_taxon_map[taxon.rank].append(taxon)
        return self._id_to_taxon_map

    def __eq__(self, other: "TaxonomyHierarchy"):
        """
        Check if two taxonomy hierarchies are equal.
        """
        for taxon, other_taxon in zip(self, other):
            if taxon != other_taxon:
                return False
        return True

    def __iter__(self) -> Generator[Taxon, None, None]:
        """
        Breadth-first sorted iteration over the taxonomy hierarchy.
        """
        q: list[Taxon] = sorted(self.taxon_tree_head.children.values())
        while len(q) > 0:
            taxon = q.pop(0)
            yield taxon
            q += sorted(taxon.children.values())

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
