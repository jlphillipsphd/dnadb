from dataclasses import dataclass, replace
import io
from itertools import chain
import json
import numpy as np
import numpy.typing as npt
from pathlib import Path
import re
from typing import Dict, Generator, Iterable, List, Optional, overload, Tuple, Union

from .db import DbFactory, DbWrapper
from .types import int_t
from .utils import open_file, sort_dict

RANKS = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
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

# Taxonomy ID Map ----------------------------------------------------------------------------------

class TaxonomyIdMap:
    """
    A bidirectional map between taxonomy labels and integer IDs.
    """
    @classmethod
    def deserialize(cls, id_map_json_bytes: Union[str, bytes, bytearray]) -> "TaxonomyIdMap":
        """
        Deserialize a taxonomy ID map from a bytes object.
        """
        id_map = json.loads(id_map_json_bytes)
        return cls(id_map)

    @classmethod
    def from_db(cls, db: Union[TaxonomyDb, Iterable[TaxonomyDb]]) -> "TaxonomyIdMap":
        """
        Create a taxonomy ID map from the given taxonomy database(s).
        """
        if isinstance(db, TaxonomyDb):
            db = [db]
        return cls(chain(*(d.labels() for d in db)))

    def __init__(self, taxonomy_labels: Optional[Iterable[str]] = None):
        self.id_to_label_map: List[str] = []
        self.label_to_id_map: Dict[str, int] = {}
        if taxonomy_labels:
            return self.add_taxonomies(taxonomy_labels)

    def add_taxonomy(self, taxonomy_label: str):
        """
        Add a taxonomy label to the ID map.
        """
        if taxonomy_label in self.label_to_id_map:
            return
        self.label_to_id_map[taxonomy_label] = len(self.label_to_id_map)
        self.id_to_label_map.append(taxonomy_label)

    def add_taxonomies(self, taxonomy_labels: Iterable[str]):
        """
        Add a set of taxonomy labels
        """
        for label in taxonomy_labels:
            self.add_taxonomy(label)

    def labels(self) -> Generator[str, None, None]:
        """
        Get the taxonomy labels.
        """
        yield from self.label_to_id_map

    def has_label(self, label: str) -> bool:
        """
        Check if the ID map has a given taxonomy label.
        """
        return label in self.label_to_id_map

    def id_to_label(self, label_id: int) -> str:
        """
        Get the taxonomy label for a given label ID.
        """
        return self.id_to_label_map[label_id]

    def label_to_id(self, label: str) -> int:
        """
        Get the taxonomy label ID for a given label.
        """
        return self.label_to_id_map[label]

    def __eq__(self, other: "TaxonomyIdMap"):
        """
        Check if two taxonomy ID maps are equal.
        """
        return self.id_to_label_map == other.id_to_label_map

    @overload
    def __getitem__(self, key: str) -> int:
        ...

    @overload
    def __getitem__(self, key: int) -> str:
        ...

    def __getitem__(self, key: Union[str, int]) -> Union[str, int]:
        if isinstance(key, str):
            return self.label_to_id(key)
        return self.id_to_label(key)

    def __iter__(self) -> Iterable[str]:
        return iter(self.label_to_id_map)

    def __len__(self) -> int:
        return len(self.id_to_label_map)

    def serialize(self) -> bytes:
        return json.dumps(self.id_to_label_map).encode()

    def display(self, max: Optional[int] = None):
        print(str(self))
        if len(self) == 0:
            print("  Empty")
            return
        n = len(self) if max is None else min(len(self), max)
        spacing = int(np.log10(n)) + 1
        for i, label in enumerate(self.id_to_label_map[:n]):
            print(f"  {i:>{spacing}}: {label}")

    def __str__(self) -> str:
        return f"TaxonomyIdMap({len(self)})"

    def __repr__(self) -> str:
        return str(self)

# Taxonomy Hierarchy -------------------------------------------------------------------------------

class TaxonomyIdTree:
    @classmethod
    def deserialize(cls, taxonomy_id_tree_bytes: Union[str, bytes]):
        deserialized = json.loads(taxonomy_id_tree_bytes)
        tree = cls(deserialized["depth"], deserialized["pad"])
        tree.tree = deserialized["tree"]
        return tree

    @classmethod
    def from_db(cls, db: TaxonomyDb, depth: int = 7) -> "TaxonomyIdTree":
        """
        Create a taxonomy hierarchy from a taxonomy database.
        """
        tree = cls(depth)
        tree.add_labels(db)
        return tree

    def __init__(self, depth: int = 7, pad: bool = True):
        self.depth = depth
        self.pad = pad
        self.tree = {}
        self._id_to_taxons_map: Union[Tuple[List[Tuple[str, ...]], ...], None]
        self._taxons_to_id_map: Union[Dict[Tuple[str, ...], int], None]
        self._mark_dirty()

    def add_taxons(self, taxons: Tuple[str, ...]):
        """
        Add a taxon hierarchy to the tree.
        """
        self._mark_dirty()
        head = self.tree
        if len(taxons) < self.depth and self.pad:
            taxons += ('',)*(self.depth - len(taxons))
        for taxon in taxons[:self.depth]:
            if not self.pad and taxon == "":
                break
            if taxon not in head:
                head[taxon] = {}
            head = head[taxon]

    def add_label(self, label: str):
        """
        Add a taxonomy label to the tree.
        """
        self.add_taxons(split_taxonomy(label, keep_empty=self.pad))

    def add_labels(self, labels: Iterable[str]):
        """
        Add a list of taxonomy labels to the tree.
        """
        for label in labels:
            self.add_label(label)

    def add_entry(self, entry: TaxonomyEntry):
        """
        Add a taxonomy entry to the tree.
        """
        self.add_label(entry.label)

    def add_entries(self, entries: Iterable[TaxonomyEntry]):
        """
        Add taxonomy entries to the tree.
        """
        for entry in entries:
            self.add_entry(entry)

    def has_entry(self, entry: TaxonomyEntry) -> bool:
        """
        Check if the tree contains the given taxonomy entry.
        """
        return self.has_label(entry.label)

    def has_label(self, label: str) -> bool:
        """
        Check if the tree contains the given taxonomy label.
        """
        return self.has_taxons(split_taxonomy(label, keep_empty=self.pad))

    def has_taxons(self, taxons: Tuple[str, ...]) -> bool:
        """
        Check if the tree contains the given hierarchy of taxons.
        """
        return taxons in self.taxons_to_id_map

    def id_to_label(self, depth: int, identifier: int) -> str:
        """
        Map an identifier to a taxonomy label.
        """
        return join_taxonomy(self.id_to_taxons(depth, identifier))

    def id_to_taxons(self, depth: int, identifier: int) -> Tuple[str, ...]:
        """
        Map a taxon identifier to its taxon hierarchy tuple.
        """
        return self.id_to_taxons_map[depth][identifier]

    def label_to_id(self, label: str) -> int:
        """
        Map the given taxonomy label to a hierarchy integer identifier.
        """
        return self.taxons_to_id(split_taxonomy(label, keep_empty=False))

    def taxons_to_id(self, taxons: Tuple[str, ...]) -> int:
        """
        Map the given taxonomy hierarchy tuple to an integer identifier for the hierarchy.
        """
        return self.taxons_to_id_map[taxons]

    def reduce_entry(self, entry: TaxonomyEntry) -> TaxonomyEntry:
        """
        Reduce the given taxonomy entry to a valid label in this tree.
        """
        return replace(entry, label=self.reduce_taxonomy(entry.label))

    def reduce_label(self, label: str) -> str:
        """
        Reduce the given taxonomy label to a valid label in this tree.
        """
        taxons = split_taxonomy(label, keep_empty=self.pad)
        return join_taxonomy(self.reduce_taxons(taxons), depth=len(taxons))

    def reduce_taxons(self, taxons: Tuple[str, ...]) -> Tuple[str, ...]:
        """
        Reduce the given taxon hierarchy to a valid taxon hierarchy in this tree.
        """
        head = self.tree
        result: Tuple[str, ...] = ()
        for taxon in taxons:
            if taxon not in head:
                break
            result += (taxon,)
            head = head[taxon]
        return result

    def tokenize_label(self, taxonomy: str) -> npt.NDArray[np.int32]:
        """
        Tokenize the taxonomy label into a a tuple of taxon integer IDs

        Args:
            taxonomy (str): The taxonomy label to tokenize (e.g. "k__Bacteria; ...").

        Returns:
            np.ndarray[np.int32]: The tokenized taxonomy.
        """
        return self.tokenize_taxons(split_taxonomy(taxonomy, keep_empty=True)[:self.depth])

    def tokenize_taxons(self, taxons: Tuple[str, ...]) -> npt.NDArray[np.int32]:
        """
        Tokenize the taxonomy tuple into a a tuple of taxon integer IDs

        Args:
            taxons (Tuple[str, ...]): The taxonomy tuple to tokenize.

        Returns:
            np.ndarray[np.int32]: The tokenized taxonomy.
        """
        result = np.empty(len(taxons), np.int32)
        for i in range(len(taxons)):
            result[i] = self.taxons_to_id_map[taxons[:i+1]]
        return result

    def detokenize_label(self, taxon_tokens: npt.NDArray[np.int32]) -> str:
        """
        Detokenize the taxonomy tokens into a taxonomy label.

        Args:
            taxon_tokens (npt.NDArray[np.int64]): The taxonomy tokens.

        Returns:
            str: The detokenized taxonomy label.
        """
        return join_taxonomy(self.detokenize_taxons(taxon_tokens), depth=self.depth)

    def detokenize_taxons(self, taxon_tokens: npt.NDArray[np.int32]) -> Tuple[str, ...]:
        """
        Detokenize the taxonomy tokens into a taxonomy tuple.

        Args:
            taxon_tokens (npt.NDArray[np.int64]): The taxonomy tokens.

        Returns:
            Tuple[str, ...]: The detokenized taxonomy tuple.
        """
        i = len(taxon_tokens) - 1
        return self.id_to_taxons_map[i][taxon_tokens[i]]

    def _mark_dirty(self):
        """
        Require a rebuild ofthe identifier maps.
        """
        self._id_to_taxons_map = None
        self._taxons_to_id_map = None

    def update(self, other_taxonomy_id_tree: "TaxonomyIdTree"):
        """
        Merge another tree into this tree.
        """
        def merge(a, b, depth: int = 0):
            if depth >= self.depth:
                return
            for taxon, children in b.items():
                if taxon not in a:
                    a[taxon] = {}
                merge(a[taxon], children, depth+1)
        merge(self.tree, other_taxonomy_id_tree.tree)

    def build(self):
        """
        Build the identifier maps given the current tree state.
        """
        self._id_to_taxons_map = tuple([] for _ in range(self.depth))
        self._taxons_to_id_map = {}
        sort_dict(self.tree)
        stack = [((), self.tree)]
        while len(stack) > 0:
            taxons, head = stack.pop()
            depth = len(taxons)
            s = []
            for taxon, children in head.items():
                next_taxons = taxons + (taxon,)
                identifier = len(self._id_to_taxons_map[depth])
                self._id_to_taxons_map[depth].append(next_taxons)
                self._taxons_to_id_map[next_taxons] = identifier
                sort_dict(children)
                s.append((next_taxons, children))
            stack += reversed(s)

    def serialize(self) -> bytes:
        """
        Serialize the tree into bytes.
        """
        return bytes(json.dumps({
            "depth": self.depth,
            "pad": self.pad,
            "tree": self.tree
        }).encode())

    def copy(self) -> "TaxonomyIdTree":
        """
        Create a copy of the tree.
        """
        return self.deserialize(self.serialize())

    @property
    def id_to_taxons_map(self) -> Tuple[List[Tuple[str, ...]], ...]:
        if self._id_to_taxons_map is None:
            self.build()
        return self._id_to_taxons_map

    @property
    def taxons_to_id_map(self) -> Dict[Tuple[str, ...], int]:
        if self._taxons_to_id_map is None:
            self.build()
        return self._taxons_to_id_map

    def __eq__(self, other: "TaxonomyIdTree"):
        return self.id_to_taxons_map == other.id_to_taxons_map

    def __iter__(self) -> Iterable[str]:
        return iter(self.id_to_taxons_map[-1])

    def __str__(self):
        return f"TaxonomyIdTree(depth={self.depth}, num_labels={len(self.id_to_taxons_map[-1])})"

    def __repr__(self):
        return str(self)

def entries(
    taxonomy: Union[io.TextIOBase, Iterable[TaxonomyEntry], str, Path],
    has_header: bool = False
) -> Iterable[TaxonomyEntry]:
    """
    Create an Iterable over a taxonomy file or iterable of taxonomy entries.
    """
    if isinstance(taxonomy, (str, Path)):
        with open_file(taxonomy, 'r') as buffer:
            yield from read(buffer, has_header=has_header)
    elif isinstance(taxonomy, io.TextIOBase):
        yield from read(taxonomy, has_header=has_header)
    elif has_header:
        yield from taxonomy[1:]
    else:
        yield from taxonomy


def read(buffer: io.TextIOBase, has_header: bool) -> Generator[TaxonomyEntry, None, None]:
    """
    Read taxonomies from a tab-separated file (TSV)
    """
    iterator = iter(buffer)
    if has_header:
        next(iterator)
    for line in iterator:
        identifier, taxonomy = line.rstrip().split('\t')
        yield TaxonomyEntry(identifier, taxonomy)


def write(buffer: io.TextIOBase, entries: Iterable[TaxonomyEntry]):
    """
    Write taxonomy entries to a tab-separate file (TSV)
    """
    for entry in entries:
        buffer.write(f"{entry.identifier}\t{entry.label}\n")
