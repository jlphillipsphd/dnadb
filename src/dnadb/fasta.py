from dataclasses import dataclass
from functools import singledispatchmethod
import io
import numpy as np
from pathlib import Path
from typing import Generator, Iterable, Tuple, Union

from .db import DbFactory, DbWrapper
from .taxonomy import TaxonomyEntry
from .utils import open_file

_int_t = Union[int, np.int32, np.int64]

@dataclass(frozen=True, order=True)
class FastaEntry:
    """
    A container class to represent a FASTA entry
    """
    __slots__ = ("identifier", "sequence", "extra")

    identifier: str
    sequence: str
    extra: str

    @classmethod
    def deserialize(cls, entry: bytes) -> "FastaEntry":
        """
        Deserialize a FASTA entry from a byte string
        """
        return cls(*entry.decode().split('\x00'))

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

    def serialize(self) -> bytes:
        return "\x00".join((self.identifier, self.sequence, self.extra)).encode()

    def __str__(self):
        header_line = f"{self.identifier} {self.extra}".rstrip()
        return f">{header_line}\n{self.sequence}"


class FastaDbFactory(DbFactory):
    """
    A factory for creating LMDB-backed databases of FASTA entries.
    """
    __slots__ = ("num_entries", "has_ambiguous_bases")

    def __init__(self, path: Union[str, Path], chunk_size: int = 10000):
        super().__init__(path, chunk_size)
        self.num_entries = np.int32(0)
        self.has_ambiguous_bases = False

    def write_entry(self, entry: FastaEntry):
        """
        Create a new FASTA LMDB database from a FASTA file.
        """
        self.write(f"id_{entry.identifier}", np.int32(self.num_entries).tobytes())
        self.write(str(self.num_entries), entry.serialize())
        self.num_entries += 1

    def write_entries(self, entries: Iterable[FastaEntry]):
        for entry in entries:
            self.write_entry(entry)

    def before_close(self):
        self.write("length", self.num_entries.tobytes())
        super().before_close()


class FastaDb(DbWrapper):
    """
    An LMDB-backed database of FASTA entries.
    """
    def __init__(self, fasta_db_path: Union[str, Path]):
        super().__init__(fasta_db_path)
        self.length = np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]

    def __len__(self):
        return self.length

    @singledispatchmethod
    def __contains__(self, sequence_index: _int_t) -> bool:
        return str(sequence_index) in self.db

    @__contains__.register
    def _(self, sequence_id: str) -> bool:
        return f"id_{sequence_id}" in self.db

    @__contains__.register
    def _(self, entry: FastaEntry) -> bool:
        return entry.identifier in self

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    @singledispatchmethod
    def __getitem__(self, sequence_index: _int_t) -> FastaEntry:
        return FastaEntry.deserialize(self.db[str(sequence_index)])

    @__getitem__.register
    def _(self, sequence_id: str) -> FastaEntry:
        index = np.frombuffer(self.db[f"id_{sequence_id}"], dtype=np.int32, count=1)[0]
        return self[index]


def entries(
    sequences: Union[io.TextIOBase, Iterable[FastaEntry], str, Path]
) -> Iterable[FastaEntry]:
    """
    Create an iterator over a FASTA file or iterable of FASTA entries.
    """
    if isinstance(sequences, (str, Path)):
        with open_file(sequences, 'r') as buffer:
            yield from read(buffer)
    elif isinstance(sequences, io.TextIOBase):
        yield from read(sequences)
    else:
        yield from sequences


def entries_with_taxonomy(
    sequences: Iterable[FastaEntry],
    taxonomies: Iterable[TaxonomyEntry],
) -> Generator[Tuple[FastaEntry, TaxonomyEntry], None, None]:
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


def read(buffer: io.TextIOBase) -> Generator[FastaEntry, None, None]:
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


def write(buffer: io.TextIOBase, entries: Iterable[FastaEntry]) -> int:
    """
    Write entries to a FASTA file.
    """
    bytes_written = 0
    for entry in entries:
        bytes_written += buffer.write(str(entry) + '\n')
    return bytes_written
