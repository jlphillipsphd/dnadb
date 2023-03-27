from functools import cache
from lmdbm import Lmdb
import numpy as np
from pathlib import Path
import re
from typing import Generator, Iterable, TextIO

from .utils import open_file

TAXON_PREFIXES = "kpcofgs"

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


def unique_labels(entries: Iterable["TaxonomyEntry"]) -> Generator["TaxonomyEntry", None, None]:
    """
    Iterate over all unique taxonomy labels.
    """
    m: dict[str, TaxonomyEntry] = {}
    for entry in entries:
        if entry.label in m:
            continue
        m[entry.label] = entry
        yield entry


def unique_taxons(entries: Iterable["TaxonomyEntry"], depth=7) -> list[set[str]]:
    """
    Pull each taxon as a set
    """
    taxon_sets: list[set[str]] = [set() for _ in range(depth)]
    for entry in entries:
        for taxon_set, taxon in zip(taxon_sets, entry.taxons(depth)):
            taxon_set.add(taxon)
    return taxon_sets

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


class TaxonomyDb:
    @classmethod
    def create(cls, taxonomy_path: str|Path, taxonomy_db_path: str|Path, chunk_size=10000) -> None:
        db = Lmdb.open(str(taxonomy_db_path), 'n')
        chunk: dict[str, bytes] = {}
        i: int = 0
        for i, entry in enumerate(entries(taxonomy_path)):
            chunk[entry.identifier] = entry.serialize()
            if i > 0 and i % chunk_size == 0:
                db.update(chunk)
                chunk.clear()
        db.update(chunk)
        db["length"] = np.int32(i + 1).tobytes()
        db.close()

    def __init__(self, taxonomy_db_path: str|Path):
        super().__init__()
        self.db = Lmdb.open(str(taxonomy_db_path))

    @cache
    def __len__(self):
        return np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]

    def __getitem__(self, sequence_index: int) -> TaxonomyEntry:
        return TaxonomyEntry.deserialize(self.db[str(sequence_index)])


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
