import bisect
from dataclasses import dataclass, field, replace
import enum
from functools import singledispatchmethod
import io
import json
import numpy as np
import numpy.typing as npt
from pathlib import Path
import re
from tqdm import tqdm
from typing import Dict, Generator, Iterable, Iterator, List, Literal, Optional, overload, Tuple, TypeVar, Union

from .db import DbFactory, DbWrapper
from .fasta import FastaDb, FastaEntry
from .utils import open_file, sort_dict

RANKS = ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
RANK_PREFIXES = ''.join(rank[0] for rank in RANKS).lower()

# Utility Functions --------------------------------------------------------------------------------

def is_taxonomy(taxonomy: str) -> bool:
    """
    Check if a string is a valid taxonomy label.
    """
    return bool(re.match(r"^\w__[^;]+(;\s*\w__[^;]*)*;?$", taxonomy))


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


@overload
def taxonomy_parent(taxonomy: str) -> str: ...
@overload
def taxonomy_parent(taxonomy: Tuple[str, ...]) -> Tuple[str, ...]: ...
def taxonomy_parent(taxonomy: Union[str, Tuple[str, ...]]) -> Union[str, Tuple[str, ...]]:
    """
    Get the parent taxonomy of a taxonomy label
    """
    if isinstance(taxonomy, str):
        taxons = split_taxonomy(taxonomy, keep_empty=True)
    else:
        taxons = taxonomy
    # find the last non-empty taxon
    for i in range(len(taxons) - 1, -1, -1):
        if len(taxons[i]) > 0:
            taxons = taxons[:i] + ("",)*(len(taxons) - i)
            break
    if isinstance(taxonomy, str):
        return join_taxonomy(taxons)
    return taxons

# Taxonomy DB --------------------------------------------------------------------------------------

class ITaxonomyEntry:
    sequence_id: str
    label: str

@dataclass(frozen=True, order=True)
class TaxonomyEntry(ITaxonomyEntry):
    sequence_id: str
    label: str

    @property
    def taxons(self):
        return split_taxonomy(self.label, keep_empty=True)

    @property
    def depth(self):
        return len(self.taxons)

    def trim(self, depth: int):
        return replace(self, label=join_taxonomy(self.taxons, depth=depth))

    def __len__(self) -> int:
        return len(self.label)

    def __str__(self) -> str:
        return "\t".join([self.sequence_id, self.label])


TaxonomyDict = Dict[str, "TaxonomyDict"]
class TaxonomyTreeFactory:
    def __init__(self, depth: int = 7):
        self.depth = depth
        self._tree: TaxonomyDict = {}

    def add_taxons(self, taxons: Tuple[str, ...]):
        taxons = (taxons + ('',)*(self.depth - len(taxons)))[:self.depth]
        tree = self._tree
        for taxon in taxons:
            if taxon not in tree:
                tree[taxon] = {}
            tree = tree[taxon]

    def add_label(self, label: str):
        self.add_taxons(split_taxonomy(label, keep_empty=True))

    def add_entry(self, entry: TaxonomyEntry):
        self.add_label(entry.label)

    def add_entries(self, entries: Iterable[TaxonomyEntry]):
        for entry in entries:
            self.add_entry(entry)

    def _sort_tree_dict(self, tree: TaxonomyDict):
        sort_dict(tree)
        for value in tree.values():
            self._sort_tree_dict(value)

    def build(self) -> "TaxonomyTree":
        self._sort_tree_dict(self._tree)
        return TaxonomyTree(self.depth, self._tree)


class TaxonomyTree:

    @dataclass(frozen=True)
    class Taxon:
        taxon_label: str = field(compare=False, hash=False)
        rank: int = field(compare=True, hash=False)
        taxon_id: int = field(compare=True)
        taxonomy_id: int = field(compare=False, hash=False)
        parent: "TaxonomyTree.Taxon" = field(compare=False, hash=False)
        children: Dict[int, "TaxonomyTree.Taxon"] = field(default_factory=dict, compare=False, hash=False)
        child_ids: Dict[str, int] = field(default_factory=dict, compare=False, hash=False)

        def __init__(self, taxon_label: str, taxon_id: int = -1, taxonomy_id: int = -1, parent: Optional["TaxonomyTree.Taxon"] = None):
            object.__setattr__(self, "taxon_label", taxon_label)
            object.__setattr__(self, "rank", parent.rank+1 if parent is not None else -1)
            object.__setattr__(self, "taxon_id", taxon_id)
            object.__setattr__(self, "taxonomy_id", taxonomy_id)
            object.__setattr__(self, "parent", parent) # None is allowed
            object.__setattr__(self, "children", {})
            object.__setattr__(self, "child_ids", {})
            if self.parent is not None:
                self.parent.children[taxon_id] = self
                self.parent.child_ids[taxon_label] = self.taxon_id

        def add_child(self, taxon_label: str, taxon_id: int, taxonomy_id: int) -> "TaxonomyTree.Taxon":
            assert taxon_label not in self.child_ids, f"Taxon {repr(taxon_label)} already exists"
            self.children[taxon_id] = TaxonomyTree.Taxon(taxon_label, taxon_id, taxonomy_id, self)
            self.child_ids[taxon_label] = taxon_id
            return self.children[taxon_id]

        @property
        def taxonomy_label(self) -> str:
            return join_taxonomy(self.taxons)

        @property
        def taxons(self) -> Tuple[str, ...]:
            head = self
            taxons: Tuple[str, ...] = ()
            while head.rank != -1:
                taxons = (head.taxon_label,) + taxons
                head = head.parent
            return taxons

        @property
        def taxon_ids(self) -> Tuple[int, ...]:
            head = self
            taxon_ids: Tuple[int, ...] = ()
            while head.rank != -1:
                taxon_ids = (head.taxon_id,) + taxon_ids
                head = head.parent
            return taxon_ids

        @property
        def taxonomy_ids(self) -> Tuple[int, ...]:
            head = self
            taxonomy_ids: Tuple[int, ...] = ()
            while head.rank != -1:
                taxonomy_ids = (head.taxonomy_id,) + taxonomy_ids
                head = head.parent
            return taxonomy_ids

        def truncate(self, rank: int) -> "TaxonomyTree.Taxon":
            assert rank >= 0 and rank <= self.rank, "Invalid rank"
            head = self
            while head.rank > rank:
                head = head.parent
            return head

        def __contains__(self, label_or_taxon_id: Union[str, int]) -> bool:
            if isinstance(label_or_taxon_id, str):
                return label_or_taxon_id in self.child_ids
            return label_or_taxon_id in self.children

        def __getitem__(self, label_or_taxon_id: Union[str, int]) -> "TaxonomyTree.Taxon":
            if isinstance(label_or_taxon_id, str):
                label_or_taxon_id = self.child_ids[label_or_taxon_id]
            return self.children[label_or_taxon_id]

        def __iter__(self) -> Iterator["TaxonomyTree.Taxon"]:
            return iter(self.children.values())

        def __len__(self) -> int:
            return len(self.children)

        def __repr__(self) -> str:
            params = [
                "taxon_label=" + repr(self.taxon_label),
            ]
            if self.rank != -1:
                params += [
                    "taxon_id=" + repr(self.taxon_id),
                    "taxonomy_id=" + repr(self.taxonomy_id),
                    "taxonomy_label=" + repr(self.taxonomy_label),
                ]
            return f"Taxon({', '.join(params)})"

    @classmethod
    def deserialize(cls, taxonomy_tree_bytes: bytes) -> "TaxonomyTree":
        return cls(**json.loads(taxonomy_tree_bytes))

    def __init__(self, depth: int, tree: TaxonomyDict):
        self.depth = depth
        self.id_to_taxon_map, self.taxon_to_id_map = self._build_taxon_id_maps(tree)
        self.tree, self.taxonomy_id_map = self._build_tree_and_taxonomy_map(tree)

    def _build_taxon_id_maps(
        self,
        tree: TaxonomyDict
    ) -> Tuple[Tuple[List[str], ...], Tuple[Dict[str, int], ...]]:
        taxon_id_to_taxon_map: Tuple[List[str], ...] = tuple([] for _ in range(self.depth))
        stack: List[Tuple[int, TaxonomyDict]] = [(0, tree)]
        while len(stack) > 0:
            rank, tree = stack.pop()
            for taxon in tree:
                index = bisect.bisect_left(taxon_id_to_taxon_map[rank], taxon)
                if index >= len(taxon_id_to_taxon_map[rank]) or taxon_id_to_taxon_map[rank][index] != taxon:
                    taxon_id_to_taxon_map[rank].insert(index, taxon)
                if len(tree[taxon]) > 0:
                    stack.append((rank + 1, tree[taxon]))
        taxon_to_taxon_id_map = tuple({t: i for i, t in enumerate(g)} for g in taxon_id_to_taxon_map)
        return taxon_id_to_taxon_map, taxon_to_taxon_id_map

    def _build_tree_and_taxonomy_map(self, tree: TaxonomyDict) -> Tuple[Taxon, Tuple[List[Taxon], ...]]:
        taxonomy_id_to_taxon_map = tuple([] for _ in range(self.depth))
        root = TaxonomyTree.Taxon("Root")
        stack: List[Tuple[TaxonomyTree.Taxon, TaxonomyDict]] = [(root, tree)]
        while len(stack) > 0:
            parent, head = stack.pop()
            s = []
            for label in head:
                taxon_id = self.taxon_to_id_map[parent.rank + 1][label]
                taxonomy_id = len(taxonomy_id_to_taxon_map[parent.rank+1])
                taxon = parent.add_child(label, taxon_id, taxonomy_id)
                taxonomy_id_to_taxon_map[parent.rank+1].append(taxon)
                assert taxonomy_id_to_taxon_map[parent.rank+1][taxonomy_id] == taxon
                if len(head[taxon.taxon_label]) > 0:
                    s.append((taxon, head[taxon.taxon_label]))
            stack += reversed(s)
        return root, taxonomy_id_to_taxon_map

    def reduce_entry(self, label: TaxonomyEntry) -> TaxonomyEntry:
        return replace(label, label=self.reduce_label(label.label))

    def reduce_label(self, label: str) -> str:
        return join_taxonomy(self.reduce_taxons(split_taxonomy(label, keep_empty=True)))

    def reduce_taxons(self, taxons: Tuple[str, ...], pad: bool = True) -> Tuple[str, ...]:
        count = 0
        head = self.tree
        for taxon in taxons:
            if taxon not in head.child_ids:
                break
            count += 1
            head = head[taxon]
        if pad:
            return taxons[:count] + ("",)*(len(taxons) - count)
        return taxons[:count]

    def reduce_taxonomy(
        self,
        taxonomy: Union[TaxonomyEntry, str, int, Tuple[str, ...], Tuple[int, ...]]
    ) -> Union["TaxonomyTree.Taxon", None]:
        if isinstance(taxonomy, int):
            return self.taxonomy_id_map[-1][taxonomy]
        if isinstance(taxonomy, TaxonomyEntry):
            taxonomy = taxonomy.label
        if isinstance(taxonomy, str):
            taxonomy = split_taxonomy(taxonomy, keep_empty=True)
        head = self.tree
        for taxon in taxonomy:
            if taxon not in head:
                break
            head = head[taxon]
        if head.rank == -1:
            return None
        return head

    def has_taxonomy(
        self,
        taxonomy: Union[TaxonomyEntry, str, int, Tuple[str, ...], Tuple[int, ...]]
    ) -> bool:
        if isinstance(taxonomy, int):
            return taxonomy < len(self.taxonomy_id_map[self.depth - 1])
        if isinstance(taxonomy, TaxonomyEntry):
            taxonomy = taxonomy.label
        if isinstance(taxonomy, str):
            taxonomy = split_taxonomy(taxonomy, keep_empty=True)
        head = self.tree
        try:
            for taxon in taxonomy:
                head = head[taxon]
            if head.rank == -1 > 0:
                raise KeyError()
        except KeyError:
            return False
        return True

    def taxonomy(
        self,
        taxonomy: Union[TaxonomyEntry, str, int, Tuple[str, ...], Tuple[int, ...]]
    ) -> "TaxonomyTree.Taxon":
        if isinstance(taxonomy, int):
            return self.taxonomy_id_map[-1][taxonomy]
        if isinstance(taxonomy, TaxonomyEntry):
            taxonomy = taxonomy.label
        if isinstance(taxonomy, str):
            taxonomy = split_taxonomy(taxonomy, keep_empty=True)
        head = self.tree
        try:
            for taxon in taxonomy:
                head = head[taxon]
            if head.rank == -1 > 0:
                raise KeyError()
        except KeyError:
            raise KeyError(f"Taxonomy {taxonomy} not found")
        return head

    def serialize(self) -> bytes:
        tree: TaxonomyDict = {}
        stack: List[Tuple[TaxonomyTree.Taxon, TaxonomyDict]] = [(self.tree, tree)]
        while len(stack) > 0:
            head, tree_head = stack.pop()
            for child in head:
                tree_head[child.taxon_label] = {}
                stack.append((child, tree_head[child.taxon_label]))
        return json.dumps(dict(depth=self.depth, tree=tree)).encode()

    def __contains__(
        self,
        taxonomy: Union[TaxonomyEntry, str, int, Tuple[str, ...], Tuple[int, ...]]
    ) -> bool:
        return self.has_taxonomy(taxonomy)

    def __getitem__(
        self,
        taxonomy: Union[TaxonomyEntry, str, int, Tuple[str, ...], Tuple[int, ...]]
    ) -> "TaxonomyTree.Taxon":
        return self.taxonomy(taxonomy)

    def __len__(self) -> int:
        return len(self.taxonomy_id_map[-1])

    def __iter__(self) -> Iterator["TaxonomyTree.Taxon"]:
        return iter(self.taxonomy_id_map[-1])

    def __eq__(self, other: "TaxonomyTree"):
        return self.serialize() == other.serialize()

    def sample(self, shape: Union[int, Tuple[int, ...]], rng: np.random.Generator) -> np.ndarray:
        result = np.empty(np.product(shape), dtype=object)
        result[:] = list(map(
            lambda i: self.taxonomy_id_map[-1][i],
            rng.choice(len(self), len(result), replace=True)))
        return result.reshape(shape)


@dataclass(frozen=True)
class TaxonomyDbEntry(ITaxonomyEntry):
    db: "TaxonomyDb"
    sequence_index: int
    label_id: int

    @property
    def fasta_entry(self) -> FastaEntry:
        assert self.db.fasta_db is not None, "FASTA DB is not available."
        return self.db.fasta_db[self.sequence_index]

    @property
    def label(self):
        return self.taxonomy.taxonomy_label

    @property
    def sequence_id(self) -> str:
        return self.db.sequence_index_to_id(self.sequence_index)

    @property
    def taxonomy(self) -> "TaxonomyTree.Taxon":
        return self.db.tree.taxonomy_id_map[-1][self.label_id]

    def __repr__(self) -> str:
        return f"TaxonomyDbEntry(sequence_id={repr(self.sequence_id)}, label={repr(self.label)})"


class TaxonomyDbFactory(DbFactory):
    def __init__(self, path: Union[str, Path], fasta_db: FastaDb, depth: int = 7):
        super().__init__(path)
        self.depth = depth
        self.fasta_db = fasta_db
        self.sequences: Dict[str, List[int]] = {}
        self.num_sequences = 0

    def write_sequence(self, sequence_id: str, label: str):
        sequence_index = self.fasta_db.sequence_id_to_index(sequence_id)
        if label not in self.sequences:
            self.sequences[label] = []
        bisect.insort(self.sequences[label], sequence_index)
        self.num_sequences += 1

    def write_entry(self, entry: TaxonomyEntry):
        self.write_sequence(entry.sequence_id, entry.label)

    def write_entries(self, entries: Iterable[TaxonomyEntry]):
        for entry in entries:
            self.write_entry(entry)

    def before_close(self):
        tree = TaxonomyTreeFactory(self.depth)
        for label in self.sequences:
            tree.add_label(label)
        tree = tree.build()
        self.write("fasta_uuid", self.fasta_db.uuid.bytes)
        self.write("tree", tree.serialize())
        self.write("num_sequences", np.int32(self.num_sequences).tobytes())
        self.write("num_labels", np.int32(len(self.sequences)).tobytes())
        for label, sequence_indices in tqdm(self.sequences.items(), desc="Writing labels disk..."):
            taxonomy_id = tree.taxonomy(label).taxonomy_id
            self.write(f"sequences_{taxonomy_id}", np.array(sequence_indices, dtype=np.int32).tobytes())
            for sequence_index in sequence_indices:
                sequence_id = self.fasta_db.index_to_sequence_id(sequence_index)
                self.write(str(sequence_index), np.int32(taxonomy_id).tobytes())
                self.write(f"sequence_index_{sequence_index}", sequence_id.encode())
                self.write(f"sequence_{sequence_id}", np.int32(sequence_index).tobytes())
        return super().before_close()


class TaxonomyDb(DbWrapper):

    # sequence_index -> label_id
    # sequence_index_x -> sequence_id
    # sequence_x -> sequence_index

    class InMemory(enum.Flag):
        Nothing = 0
        SequencesWithTaxonomy = enum.auto()
        SequenceLabels = enum.auto()
        SequenceIdMaps = enum.auto()
        All = SequencesWithTaxonomy | SequenceLabels | SequenceIdMaps

    _sequences_with_label: List[npt.NDArray[np.int32]]|None = None # label_id -> sequence_indices
    _sequence_labels: Dict[int, int]|None = None                   # sequence_index -> label_id
    _sequence_id_to_index: Dict[str, int]|None = None              # sequence_id -> sequence_index
    _sequence_index_to_id: Dict[int, str]|None = None              # sequence_index -> sequence_id

    def __init__(
        self,
        path: Union[str, Path],
        fasta_db: Optional[FastaDb] = None,
        in_memory: InMemory = InMemory.Nothing
    ):
        super().__init__(path)
        self.fasta_db = fasta_db
        if self.fasta_db is not None:
            assert self.fasta_db.uuid.bytes == self.db["fasta_uuid"], "FASTA DB UUID does not match"
        self.tree = TaxonomyTree.deserialize(self.db["tree"])
        self.num_sequences = np.frombuffer(self.db["num_sequences"], dtype=np.int32)[0]
        self.num_labels = np.frombuffer(self.db["num_labels"], dtype=np.int32)[0]

        if TaxonomyDb.InMemory.SequencesWithTaxonomy in in_memory:
            self._sequences_with_label = []
            for i in range(self.num_labels):
                self._sequences_with_label.append(np.frombuffer(self.db[f"sequences_{i}"], dtype=np.int32))

        if TaxonomyDb.InMemory.SequenceLabels in in_memory:
            self._sequence_labels = {}
            for i in range(self.num_sequences):
                self._sequence_labels[i] = np.frombuffer(self.db[str(i)], dtype=np.int32)[0]

        if TaxonomyDb.InMemory.SequenceIdMaps in in_memory:
            self._sequence_id_to_index = {}
            self._sequence_index_to_id = {}
            for i in range(self.num_sequences):
                sequence_id = self.db[f"sequence_index_{i}"].decode()
                self._sequence_id_to_index[sequence_id] = i
                self._sequence_index_to_id[i] = sequence_id

    def count(
        self,
        taxonomy: Union[TaxonomyEntry, str, int, Tuple[str, ...], Tuple[int, ...]]
    ) -> int:
        return sum(1 for _ in self.sequences_with_taxonomy(taxonomy))

    def has_taxonomy(
        self,
        taxonomy: Union[TaxonomyEntry, str, int, Tuple[str, ...], Tuple[int, ...]]
    ) -> bool:
        return self.tree.has_taxonomy(taxonomy)

    def labels(self, depth: int = -1) -> Iterator[TaxonomyTree.Taxon]:
        return iter(self.tree.taxonomy_id_map[depth])

    def sequence_index_to_id(self, sequence_index: int) -> str:
        if self._sequence_index_to_id is not None:
            return self._sequence_index_to_id[sequence_index]
        return self.db[f"sequence_index_{sequence_index}"].decode()

    def sequence_id_to_index(self, sequence_id: str) -> int:
        if self._sequence_id_to_index is not None:
            return self._sequence_id_to_index[sequence_id]
        return np.frombuffer(self.db[f"sequence_{sequence_id}"], dtype=np.int32)[0]

    @singledispatchmethod
    def __contains__(self, sequence_index: int):
        return sequence_index < self.num_sequences

    @__contains__.register
    def _(self, sequence_id: str):
        if self._sequence_id_to_index is not None:
            return sequence_id in self._sequence_id_to_index
        return f"sequence_{sequence_id}" in self.db

    def sequence_indices_with_taxonomy_id(self, taxonomy_id: int) -> npt.NDArray[np.int32]:
        if self._sequences_with_label is not None:
            return self._sequences_with_label[taxonomy_id]
        return np.frombuffer(self.db[f"sequences_{taxonomy_id}"], dtype=np.int32)

    @singledispatchmethod
    def sequences_with_taxonomy(self, taxonomy_id: int) -> Generator[TaxonomyDbEntry, None, None]:
        for i in self.sequence_indices_with_taxonomy_id(taxonomy_id):
            yield TaxonomyDbEntry(self, i, taxonomy_id)

    @sequences_with_taxonomy.register
    def _(self, taxonomy: TaxonomyTree.Taxon) -> Generator[TaxonomyDbEntry, None, None]:
        return self.sequences_with_taxonomy(taxonomy.taxonomy_id)

    @sequences_with_taxonomy.register
    def _(self, taxonomy_label: str) -> Generator[TaxonomyDbEntry, None, None]:
        return self.sequences_with_taxonomy(self.tree.taxonomy(taxonomy_label))

    @singledispatchmethod
    def __getitem__(self, sequence_index: int) -> TaxonomyDbEntry:
        if self._sequence_labels is not None:
            taxonomy_id = self._sequence_labels[sequence_index]
        else:
            taxonomy_id = np.frombuffer(self.db[f"{sequence_index}"], dtype=np.int32)[0]
        return TaxonomyDbEntry(self, sequence_index, taxonomy_id)

    @__getitem__.register
    def _(self, sequence_id: str) -> TaxonomyDbEntry:
        return self[self.sequence_id_to_index(sequence_id)]

    def __len__(self) -> int:
        return self.num_sequences

    def __iter__(self) -> Iterator[TaxonomyDbEntry]:
        for i in range(self.num_sequences):
            yield self[i]

    def sample(self, shape: Union[int, Tuple[int, ...]], rng: np.random.Generator) -> np.ndarray:
        """
        Sample entries weighted by the number of sequences in each taxonomy.
        """
        result = np.empty(np.product(shape), dtype=object)
        result[:] = list(map(
            lambda i: TaxonomyDbEntry(self, rng.choice(self.sequence_indices_with_taxonomy_id(i)), i), # type: ignore
            rng.choice(self.num_labels, len(result), replace=True)))
        return result.reshape(shape)

T = TypeVar("T", bound=ITaxonomyEntry)
def entries(
    taxonomy: Union[io.TextIOBase, Iterable[T], str, Path],
    header: bool|Literal["auto"] = "auto"
) -> Generator[TaxonomyEntry, None, None]:
    """
    Create an Iterable over a taxonomy file or iterable of taxonomy entries.
    """
    if isinstance(taxonomy, (str, Path)):
        with open_file(taxonomy, 'r') as buffer:
            yield from read(buffer, header=header)
    elif isinstance(taxonomy, io.TextIOBase):
        yield from read(taxonomy, header=header)
    else:
        yield from map(lambda entry: TaxonomyEntry(entry.sequence_id, entry.label), taxonomy)


def read(
    buffer: io.TextIOBase,
    header: bool|Literal["auto"] = "auto"
) -> Generator[TaxonomyEntry, None, None]:
    """
    Read taxonomies from a tab-separated file (TSV)
    """
    iterator = iter(buffer)
    if header == "auto":
        sequence_id, taxonomy = next(iterator).strip().split('\t')
        if is_taxonomy(taxonomy):
            yield TaxonomyEntry(sequence_id, taxonomy)
    elif header:
        next(iterator)
    for line in iterator:
        sequence_id, taxonomy = line.rstrip().split('\t')
        yield TaxonomyEntry(sequence_id, taxonomy)


def write(buffer: io.TextIOBase, entries: Iterable[ITaxonomyEntry]):
    """
    Write taxonomy entries to a tab-separate file (TSV)
    """
    for entry in entries:
        buffer.write(f"{entry.sequence_id}\t{entry.label}\n")

