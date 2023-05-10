from functools import cached_property
import io
import json
from lmdbm import Lmdb
import numpy as np
from pathlib import Path
import re
from typing import Dict, Generator, Iterable, List, Optional, Tuple, TypedDict, Set, Union

from .db import DbFactory, DbWrapper
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

# Taxonomy TSV Utilities ---------------------------------------------------------------------------

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

    [label index to fasta ids]
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
