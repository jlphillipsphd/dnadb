from functools import cache, singledispatchmethod
from lmdbm import Lmdb
import numpy as np
from pathlib import Path
from typing import Generator, Iterable, TextIO

from .taxonomy import TaxonomyEntry
from .utils import open_file

class FastaEntry:

    @classmethod
    def deserialize(cls, entry: bytes) -> "FastaEntry":
        """
        Deserialize a FASTA entry from a byte string
        """
        return cls.from_str(entry.decode())

    @classmethod
    def from_str(cls, entry: str) -> "FastaEntry":
        """
        Create a FASTA entry from a string
        """
        header, *sequence_parts = entry.split('\n')
        header_line = header[1:].rstrip().split(maxsplit=1)
        identifier = header_line[0]
        extra = header_line[1] if len(header_line) > 1 else ""
        sequence = "".join(sequence_parts)
        return cls(identifier, sequence, extra)

    """
    A container class to represent a FASTA entry
    """
    def __init__(self, identifier: str, sequence: str, extra: str = ""):
        self.identifier = identifier
        self.sequence = sequence
        self.extra = extra

    def serialize(self) -> bytes:
        return str(self).encode()

    def __str__(self):
        header_line = f"{self.identifier} {self.extra}".rstrip()
        return f">{header_line}\n{self.sequence}"

    def __repr__(self):
        return str(self)


class FastaDb:
    """
    An LMDB-backed database of FASTA entries
    """
    @classmethod
    def create(cls, fasta_entries: Iterable[FastaEntry], fasta_db_path: str|Path, chunk_size=10000):
        """
        Create a new FASTA LMDB database from a FASTA file
        """
        db = Lmdb.open(str(fasta_db_path), 'n')
        chunk: dict[str, bytes] = {}
        i: int = 0
        for i, entry in enumerate(fasta_entries):
            chunk[f"id_{entry.identifier}"] = np.int32(i).tobytes()
            chunk[str(i)] = entry.serialize()
            if i > 0 and i % chunk_size == 0:
                db.update(chunk)
                chunk.clear()
        db.update(chunk)
        db["length"] = np.int32(i + 1).tobytes()
        db.close()

    def __init__(self, fasta_db_path: str|Path):
        super().__init__()
        self.path = Path(fasta_db_path).absolute
        self.db = Lmdb.open(str(fasta_db_path))

    @cache
    def __len__(self):
        return np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]

    @singledispatchmethod
    def __getitem__(self, sequence_index: int) -> FastaEntry:
        return FastaEntry.deserialize(self.db[str(sequence_index)])

    @__getitem__.register
    def _(self, sequence_id: str) -> FastaEntry:
        index = np.frombuffer(self.db[f"id_{sequence_id}"], dtype=np.int32, count=1)[0]
        return self[index]


def entries(sequences: TextIO|Iterable[FastaEntry]|str|Path) -> Iterable[FastaEntry]:
    """
    Create an iterator over a FASTA file or iterable of FASTA entries.
    """
    if isinstance(sequences, (str, Path)):
        with open_file(sequences, 'r') as buffer:
            yield from read(buffer)
    elif isinstance(sequences, TextIO):
        yield from read(sequences)
    else:
        yield from sequences


def entries_with_taxonomy(
    sequences: Iterable[FastaEntry],
    taxonomies: Iterable[TaxonomyEntry],
) -> Generator[tuple[FastaEntry, TaxonomyEntry], None, None]:
    """
    Efficiently iterate over a FASTA file with a corresponding taxonomy file
    """
    labels = {}
    taxonomy_iterator = iter(taxonomies)
    taxonomy: TaxonomyEntry
    for sequence in sequences:
        while sequence.identifier not in labels:
            taxonomy = next(taxonomy_iterator)
            labels[taxonomy.identifier] = taxonomy
        taxonomy = labels[sequence.identifier]
        del labels[sequence.identifier]
        yield sequence, taxonomy


def read(buffer: TextIO) -> Generator[FastaEntry, None, None]:
    """
    Read entries from a FASTA file buffer.
    """
    entry_str = buffer.readline()
    for line in buffer:
        if line.startswith('>'):
            yield FastaEntry.from_str(entry_str)
            entry_str = ""
        entry_str += line
    if len(entry_str) > 0:
        yield FastaEntry.from_str(entry_str)


def write(buffer: TextIO, entries: Iterable[FastaEntry]) -> int:
    """
    Write entries to a FASTA file.
    """
    bytes_written = 0
    for entry in entries:
        bytes_written += buffer.write(str(entry) + '\n')
    return bytes_written
