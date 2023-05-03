from functools import cached_property
import io
import json
from lmdbm import Lmdb
import numpy as np
from pathlib import Path
import re
from typing import Dict, Generator, Iterable, List, Optional, Tuple, TypedDict, Set, Union

from .db import DbFactory
from .utils import open_file

TAXON_LEVEL_NAMES = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAXON_PREFIXES = ''.join(name[0] for name in TAXON_LEVEL_NAMES).lower()

class TaxonHierarchyJson(TypedDict):
    max_depth: int
    parent_map: List[Dict[str, str]]

# Utility Functions --------------------------------------------------------------------------------

def split_taxonomy(taxonomy: str, max_depth: int = 7) -> Tuple[str, ...]:
    """
    Split taxonomy label into a tuple
    """
    return tuple(re.findall(r"\w__([^;]*)", taxonomy))[:max_depth]


def join_taxonomy(taxonomy: Union[Tuple[str], List[str]], depth: int = 7) -> str:
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


def unique_taxons(entries: Iterable["TaxonomyEntry"], depth: int = 7) -> List[Set[str]]:
    """
    Pull each taxon as a set
    """
    taxon_sets: List[set[str]] = [set() for _ in range(depth)]
    for entry in entries:
        for taxon_set, taxon in zip(taxon_sets, entry.taxons(depth)):
            taxon_set.add(taxon)
    return taxon_sets

# Taxonomy TSV Utilities ---------------------------------------------------------------------------

class Taxon:
    def __init__(self, name: str, parent: Optional["Taxon"] = None):
        self.depth = (parent.depth + 1) if parent is not None else 0
        self.name = name
        self.parent = parent
        self.children: Set[Taxon] = set([])
        if self.parent:
            self.parent.children.add(self)

    def resolve(self, depth: Optional[int] = None) -> Tuple[str, ...]:
        taxon = self
        taxon_list: List[str] = []
        while taxon is not None:
            taxon_list.append(taxon.name)
            taxon = taxon.parent
        taxons: Tuple[str, ...] = tuple(reversed(taxon_list))
        if depth is not None:
            taxons = taxons[:depth] + ('',)*max(depth - len(taxons), 0)
        return taxons

    def __hash__(self):
        return hash((self.depth, self.name.casefold()))

    def __lt__(self, other: "Taxon"):
        if self.depth < other.depth:
            return True
        return self.depth == other.depth and self.name.casefold() < other.name.casefold()

    def __gt__(self, other: "Taxon"):
        if self.depth > other.depth:
            return True
        return self.depth == other.depth and self.name.casefold() > other.name.casefold()

    def __le__(self, other: "Taxon"):
        if self.depth < other.depth:
            return True
        return self.depth == other.depth and self.name.casefold() <= other.name.casefold()

    def __ge__(self, other: "Taxon"):
        if self.depth > other.depth:
            return True
        return self.depth == other.depth and self.name.casefold() >= other.name.casefold()

    def __eq__(self, other: "Taxon"):
        return self.depth == other.depth and self.name.casefold() == other.name.casefold()

    def __str__(self):
        return join_taxonomy(self.resolve(), self.depth + 1)

    def __repr__(self):
        return str(self)


class TaxonomyIdentifierMap:
    def __init__(self, taxons: Iterable[Iterable[str]]):
        self.taxon_identifiers: List[Dict[str, int]] = []
        self.taxon_labels: List[Dict[int, str]] = []
        for i, taxon_group in enumerate(taxons):
            self.taxon_identifiers.append({v: i for i, v in enumerate(sorted(taxon_group))})
            self.taxon_identifiers[-1][""] = -1
            self.taxon_labels.append({v: k for k, v in self.taxon_identifiers[i].items()})

    def encode(self, taxons: Tuple[str, ...]) -> Tuple[int, ...]:
        """
        Encode the taxon strings into integer identifiers.
        """
        return tuple(self.taxon_identifiers[i][taxon] for i, taxon in enumerate(taxons))

    def decode(self, taxons: Tuple[int, ...]) -> Tuple[str, ...]:
        """
        Decode the taxon integer identifiers into their corresponding string labels.
        """
        return tuple(self.taxon_labels[i][taxon] for i, taxon in enumerate(taxons))

    def encode_label(self, taxonomy: str) -> Tuple[int, ...]:
        """
        Encode a taxonomy label into integer identifiers.
        """
        return self.encode(split_taxonomy(taxonomy, len(self.taxon_identifiers)))

    def decode_label(self, taxonomy: Tuple[int, ...]) -> str:
        """
        Decode a taxonomy label from integer identifiers into their corresponding string labels.
        """
        return join_taxonomy(self.decode(taxonomy), len(self.taxon_identifiers))

    def encode_entry(self, taxonomy_entry: "TaxonomyEntry") -> Tuple[int, ...]:
        """
        Encode a taxonomy entry into integer identifiers.
        """
        return self.encode(taxonomy_entry.taxons(len(self.taxon_identifiers)))


class TaxonomyHierarchy:
    @classmethod
    def deserialize(cls, taxonomy_hierarchy: bytes) -> "TaxonomyHierarchy":
        return cls.from_json(json.loads(taxonomy_hierarchy.decode()))

    @classmethod
    def from_json(cls, hierarchy_json: TaxonHierarchyJson) -> "TaxonomyHierarchy":
        hierarchy = TaxonomyHierarchy(hierarchy_json["max_depth"])
        hierarchy.depth = len(hierarchy_json["parent_map"])
        for i in range(len(hierarchy.taxon_maps)):
            taxon_map = hierarchy.taxon_maps[i]
            parent_map = hierarchy_json["parent_map"][i]
            for child, parent_name in parent_map.items():
                parent = hierarchy.taxon_maps[i-1][parent_name] if i > 0 else None
                taxon_map[child] = Taxon(child, parent)
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
        max_depth = min(hierarchy.max_depth for hierarchy in hierarchy_list)
        if any(hierarchy.depth > max_depth for hierarchy in hierarchy_list):
            print(
                "Warning: Merging taxonomy hierarchies with different maximum depths.",
                f"Using depth: {max_depth}."
            )
        merged_hierarchy = TaxonomyHierarchy(max_depth)
        # Largest depth smaller than the maximum depth, or the maximum depth if no such depth exists
        merged_hierarchy.depth = min(max(h.depth for h in hierarchy_list), max_depth)
        for other_hierarchy in hierarchy_list:
            for i in range(min(other_hierarchy.depth, merged_hierarchy.depth)):
                to_map = merged_hierarchy.taxon_maps[i]
                from_map = other_hierarchy.taxon_maps[i]
                for taxon in from_map.values():
                    if taxon.name in to_map:
                        continue
                    if i > 0 and taxon.parent is not None:
                        parent = merged_hierarchy.taxon_maps[i-1][taxon.parent.name]
                    else:
                        parent = None
                    to_map[taxon.name] = Taxon(taxon.name, parent)
        return merged_hierarchy

    def __init__(self, max_depth: int = 7):
        self.max_depth = max_depth
        self.taxon_maps: List[Dict[str, Taxon]] = [] # taxon -> parent

    def add_entry(self, entry: "TaxonomyEntry"):
        self.add_taxons(entry.taxons(self.max_depth))

    def add_taxonomy(self, entry: str):
        self.add_taxons(split_taxonomy(entry))

    def add_taxons(self, taxons: Tuple[str, ...]):
        parent: Optional[Taxon] = None
        taxons = tuple(taxon for taxon in taxons[:self.max_depth] if taxon != "")
        if len(taxons) > self.depth:
            self.depth = len(taxons)
        for taxon_map, taxon_name in zip(self.taxon_maps, taxons):
            if taxon_name not in taxon_map:
                taxon_map[taxon_name] = Taxon(taxon_name, parent)
            parent = taxon_map[taxon_name]

    def is_valid(self, taxons: Union[str, Tuple[str, ...]]) -> bool:
        taxon_list = split_taxonomy(taxons) if isinstance(taxons, str) else taxons
        for taxon_map, taxon in zip(self.taxon_maps, taxon_list):
            if len(taxon) == 0:
                break
            if taxon not in taxon_map:
                return False
        return True

    def reduce_entry(self, entry: "TaxonomyEntry") -> "TaxonomyEntry":
        """
        Reduce the provided taxonomy entry to a valid taxonomy in this hierarchy.
        """
        reduced_label = join_taxonomy(self.reduce_taxons(entry.taxons()))
        return TaxonomyEntry(entry.identifier, reduced_label)

    def reduce_taxons(self, taxons: Tuple[str]) -> Tuple[str]:
        """
        Reduce the provided taxons to a valid taxonomy in this hierarchy.
        """
        reduced_taxons: List[str] = [""]*len(taxons)
        for i, (taxon_map, taxon) in enumerate(zip(self.taxon_maps, taxons)):
            if taxon not in taxon_map:
                break
            reduced_taxons[i] = taxon
        return tuple(reduced_taxons)

    # def taxonomy(self, depth: int, taxon: str) -> str:
    #     return join_taxonomy(self.taxons(depth, taxon), len(self.taxon_maps))

    # def taxons(self, depth: int, taxon: str) -> Tuple[str, ...]:
    #     return self.taxon_maps[depth][taxon].resolve(len(self.taxon_maps))

    # def taxon_identifiers(self, )

    def leaves(self):
        def find_leaves(node: Taxon):
            result = []
            for child in node.children:
                if len(child.children) == 0:
                    result.append(child)
                    continue
                result += find_leaves(child)
            return result
        return [leaf for taxon in self.taxon_maps[0].values() for leaf in find_leaves(taxon)]

    def serialize(self) -> bytes:
        json_hierarchy = self.to_json()
        return json.dumps(json_hierarchy, separators=(',', ':')).encode()

    def to_json(self) -> TaxonHierarchyJson:
        return {
            "max_depth": self.max_depth,
            "parent_map": [
                {taxon.name: taxon.parent.name if taxon.parent else "" for taxon in m.values()}
                for m in self.taxon_maps
            ]
        }

    @property
    def depth(self):
        return len(self.taxon_maps)

    @depth.setter
    def depth(self, depth: int):
        assert depth >= 1
        if depth >= self.depth:
            for i in range(self.depth, depth):
                self.taxon_maps.append({})
            return
        for taxon in self.taxon_maps[depth - 1].values():
            taxon.children.clear()
        self.taxon_maps = self.taxon_maps[:depth]

    @cached_property
    def identifier_map(self) -> TaxonomyIdentifierMap:
        return TaxonomyIdentifierMap(((taxon for taxon in m) for m in self.taxon_maps))


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

    def __init__(self, identifier: str, label: str):
        self.identifier = identifier
        self.label = label

    def taxons(self, depth: int = 7) -> Tuple[str, ...]:
        return split_taxonomy(self.label, depth)

    def serialize(self) -> bytes:
        return str(self).encode()

    def __eq__(self, other: "TaxonomyEntry"):
        return self.identifier == other.identifier \
            and self.label == other.label

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f"{self.identifier}\t{self.label}"


class TaxonomyDbFactory(DbFactory):
    """
    A factory for creating LMDB-backed databases of FASTA entries.
    """
    def __init__(self, path: Union[str, Path], max_depth: int = 7, chunk_size: int = 10000):
        super().__init__(path, chunk_size)
        self.hierarchy = TaxonomyHierarchy(max_depth)
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

    def before_close(self):
        self.write("hierarchy", self.hierarchy.serialize())
        self.write("length", self.num_entries.tobytes())
        super().before_close()


class TaxonomyDb:
    def __init__(self, taxonomy_db_path: Union[str, Path]):
        self.path = Path(taxonomy_db_path).absolute()
        self.db = Lmdb.open(str(self.path), lock=False)
        self.length = np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]

    @cached_property
    def hierarchy(self):
        return TaxonomyHierarchy.deserialize(self.db["hierarchy"])

    def __len__(self):
        return self.length

    def __contains__(self, sequence_id: str) -> bool:
        return sequence_id in self.db

    def __iter__(self):
        for key in self.db.keys():
            if key in ("length", "hierarchy"):
                continue
            yield self[key]

    def __getitem__(self, sequence_id: str) -> TaxonomyEntry:
        return TaxonomyEntry.deserialize(self.db[sequence_id])


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
