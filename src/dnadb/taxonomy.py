from functools import cache, cached_property
import json
from lmdbm import Lmdb
import numpy as np
from pathlib import Path
import re
from typing import Generator, Iterable, TextIO, TypedDict

from .db import DbFactory
from .utils import open_file

TAXON_PREFIXES = "kpcofgs"

class TaxonHierarchyJson(TypedDict):
    taxon_values: list[list[str]]
    children: dict[str, list[str]]
    parents: dict[str, str]

# Utility Functions --------------------------------------------------------------------------------

def split_taxonomy(taxonomy: str, max_depth: int = 7) -> tuple[str, ...]:
    """
    Split taxonomy label into a tuple
    """
    return tuple(re.findall(r"\w__([^;]+)", taxonomy))[:max_depth]


def join_taxonomy(taxonomy: tuple[str]|list[str], depth: int = 7) -> str:
    """
    Merge a taxonomy tuple into a string format
    """
    assert depth >= 1 and depth <= 7
    taxonomy = taxonomy[:depth] # Trim to depth
    taxonomy = tuple(taxonomy) + ("",) * (depth - len(taxonomy))
    return "; ".join([f"{TAXON_PREFIXES[i]}__{taxon}" for i, taxon in enumerate(taxonomy)])


def unique_labels(entries: Iterable["TaxonomyEntry"]) -> Generator[str, None, None]:
    """
    Iterate over all unique taxonomy labels.
    """
    m: dict[str, TaxonomyEntry] = {}
    for entry in entries:
        if entry.label in m:
            continue
        m[entry.label] = entry
        yield entry.label


def unique_taxons(entries: Iterable["TaxonomyEntry"], depth: int = 7) -> list[set[str]]:
    """
    Pull each taxon as a set
    """
    taxon_sets: list[set[str]] = [set() for _ in range(depth)]
    for entry in entries:
        for taxon_set, taxon in zip(taxon_sets, entry.taxons(depth)):
            taxon_set.add(taxon)
    return taxon_sets

# Taxonomy TSV Utilities ---------------------------------------------------------------------------

class TaxonomyHierarchy:
    @classmethod
    def deserialize(cls, taxonomy_hierarchy: bytes) -> "TaxonomyHierarchy":
        hierarchy_json: TaxonHierarchyJson = json.loads(taxonomy_hierarchy.decode())
        hierarchy = TaxonomyHierarchy(len(hierarchy_json["taxon_values"]))
        hierarchy.taxon_values = hierarchy_json["taxon_values"]
        hierarchy.children.update({t: set(c) for t, c in hierarchy_json["children"].items()})
        hierarchy.parents = hierarchy_json["parents"]
        return hierarchy

    @classmethod
    def from_entries(cls, entries: Iterable["TaxonomyEntry"], depth: int = 7):
        hierarchy = TaxonomyHierarchy(depth)
        for label in unique_labels(entries):
            hierarchy.add_taxons(split_taxonomy(label))
        return hierarchy

    @classmethod
    def from_labels(cls, labels: Iterable[str], depth: int = 7):
        hierarchy = TaxonomyHierarchy(depth)
        for label in set(labels):
            hierarchy.add_taxons(split_taxonomy(label))
        return hierarchy

    @classmethod
    def merge(cls, hierarchies: Iterable["TaxonomyHierarchy"]) -> "TaxonomyHierarchy":
        hierarchy_list = list(hierarchies)
        depth = min(hierarchy.depth for hierarchy in hierarchy_list)
        if any(hierarchy.depth > depth for hierarchy in hierarchy_list):
            print(
                "Warning: Merging taxonomy hierarchies with different depths.",
                f"Using depth: {depth}."
            )
        merged_hierarchy = TaxonomyHierarchy(depth)
        for hierarchy in hierarchy_list:
            for i in range(depth):
                for taxon in hierarchy.taxon_values[i]:
                    if taxon not in merged_hierarchy.children:
                        merged_hierarchy.taxon_values[i].append(taxon)
                        merged_hierarchy.children[taxon] = set()
                    merged_hierarchy.children[taxon].update(hierarchy.children[taxon])
        return merged_hierarchy

    def __init__(self, depth: int = 7):
        self.taxon_values: list[list[str]] = [[] for _ in range(depth)]
        self.children: dict[str, set[str]] = {}
        self.parents: dict[str, str] = {}

    def add_entry(self, entry: "TaxonomyEntry"):
        self.add_taxons(entry.taxons(self.depth))

    def add_taxons(self, taxons: tuple[str, ...]):
        parent = ""
        for taxon_level, taxon in zip(self.taxon_values, taxons):
            if taxon not in self.children:
                taxon_level.append(taxon)
                self.children[taxon] = set()
            if parent != "":
                self.children[parent].add(taxon)
                self.parents[taxon] = parent
            parent = taxon

    def is_valid(self, taxons: str|tuple[str, ...]) -> bool:
        taxon_list = split_taxonomy(taxons) if isinstance(taxons, str) else taxons
        for taxon in taxon_list:
            if taxon not in self.children:
                return False
        return True

    def reduce_entry(self, entry: "TaxonomyEntry") -> "TaxonomyEntry":
        """
        Reduce the provided taxonomy entry to a valid taxonomy in this hierarchy.
        """
        reduced_label = join_taxonomy(self.reduce_taxons(entry.taxons()))
        return TaxonomyEntry(entry.identifier, reduced_label)

    def reduce_taxons(self, taxons: tuple[str]) -> tuple[str]:
        """
        Reduce the provided taxons to a valid taxonomy in this hierarchy.
        """
        reduced_taxons: list[str] = [""]*len(taxons)
        for i, taxon in enumerate(taxons):
            if taxon not in self.children:
                break
            reduced_taxons[i] = taxon
        return tuple(reduced_taxons)

    def taxonomy(self, taxon: str) -> str:
        return join_taxonomy(self.taxons(taxon), len(self.taxon_values))

    def taxons(self, taxon: str) -> tuple[str, ...]:
        values = [taxon]
        while taxon in self.parents:
            taxon = self.parents[taxon]
            values.append(taxon)
        padding = ('0',)*(len(self.taxon_values) - len(values))
        return tuple(reversed(values)) + padding

    def serialize(self) -> bytes:
        json_hierarchy: TaxonHierarchyJson = {
            "taxon_values": self.taxon_values,
            "children": {taxon: list(children) for taxon, children in self.children.items()},
            "parents": self.parents
        }
        return json.dumps(json_hierarchy, separators=(',', ':')).encode()

    def __getitem__(self, key: str) -> set[str]:
        return self.children[key]

    @property
    def depth(self):
        return len(self.taxon_values)


class TaxonomyEntry:

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

    def __init__(self, identifier, label):
        self.identifier = identifier
        self.label = label

    def taxons(self, depth: int = 7) -> tuple[str, ...]:
        return split_taxonomy(self.label, depth)

    def serialize(self) -> bytes:
        return str(self).encode()

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f"{self.identifier}\t{self.label}"


class TaxonomyDbFactory(DbFactory):
    """
    A factory for creating LMDB-backed databases of FASTA entries.
    """
    def __init__(self, path: str|Path, chunk_size: int = 10000):
        super().__init__(path, chunk_size)
        self.hierarchy = TaxonomyHierarchy()
        self.num_entries = np.int32(0)

    def write_entry(self, entry: TaxonomyEntry):
        """
        Create a new FASTA LMDB database from a FASTA file.
        """
        self.hierarchy.add_entry(entry)
        self.write(entry.identifier, entry.serialize())
        self.num_entries += 1

    def write_entries(self, entries: Iterable[TaxonomyEntry]):
        for entry in entries:
            self.write_entry(entry)

    def close(self):
        self.write("hierarchy", self.hierarchy.serialize())
        self.write("length", self.num_entries.tobytes())
        super().close()


class TaxonomyDb:
    def __init__(self, taxonomy_db_path: str|Path):
        self.path = Path(taxonomy_db_path).absolute
        self.db = Lmdb.open(str(taxonomy_db_path))

    @cache
    def __len__(self):
        return np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]

    @cached_property
    def hierarchy(self):
        return TaxonomyHierarchy.deserialize(self.db["hierarchy"])

    def __getitem__(self, sequence_id: str) -> TaxonomyEntry:
        return TaxonomyEntry.deserialize(self.db[sequence_id])


def entries(taxonomy: TextIO|Iterable[TaxonomyEntry]|str|Path) -> Iterable[TaxonomyEntry]:
    """
    Create an iterator over a taxonomy file or iterable of taxonomy entries.
    """
    if isinstance(taxonomy, (str, Path)):
        with open_file(taxonomy, 'r') as buffer:
            yield from read(buffer)
    elif isinstance(taxonomy, TextIO):
        yield from read(taxonomy)
    else:
        yield from taxonomy


def read(buffer: TextIO) -> Generator[TaxonomyEntry, None, None]:
    """
    Read taxonomies from a tab-separated file (TSV)
    """
    for line in buffer:
        identifier, taxonomy = line.rstrip().split('\t')
        yield TaxonomyEntry(identifier, taxonomy)


def write(buffer: TextIO, entries: Iterable[TaxonomyEntry]):
    """
    Write taxonomy entries to a tab-separate file (TSV)
    """
    for entry in entries:
        buffer.write(f"{entry.identifier}\t{entry.label}\n")
