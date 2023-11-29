import bisect
from dataclasses import dataclass, field
import enum
from functools import singledispatchmethod
import io
import numpy as np
import numpy.typing as npt
from pathlib import Path
from typing import Dict, Generator, List, Iterable, Tuple, Union

from .db import DbFactory, DbWrapper
from .dna import AbstractSequenceWrapper
from .sample import ISample
from .utils import open_file

class IFastaEntry(AbstractSequenceWrapper):
    identifier: str
    extra: str

@dataclass(frozen=True, order=True)
class FastaEntry(IFastaEntry):
    """
    A container class to represent a FASTA entry
    """
    # __slots__ = ("identifier", "extra")

    identifier: str
    sequence: str
    extra: str = field(default="")

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
    __slots__ = ("num_entries")

    def __init__(self, path: Union[str, Path], chunk_size: int = 10000):
        super().__init__(path, chunk_size)
        self.num_entries = np.int32(0)

    def write_entry(self, entry: FastaEntry):
        """
        Create a new FASTA LMDB database from a FASTA file.
        """
        index = self.num_entries
        self.write(f">{entry.identifier}", np.int32(index).tobytes())
        self.write(f"{index}_id", entry.identifier.encode())
        self.write(f"{index}", entry.sequence.encode())
        self.write(f"{index}_extra", entry.extra.encode())
        self.num_entries += 1

    def write_entries(self, entries: Iterable[FastaEntry]):
        for entry in entries:
            self.write_entry(entry)

    def before_close(self):
        self.write("length", self.num_entries.tobytes())
        super().before_close()


class FastaDbEntry(IFastaEntry):

    __slots__ = ("fasta_db", "index")

    fasta_db: "FastaDb"
    index: int

    def __init__(self, fasta_db: "FastaDb", index: int):
        self.fasta_db = fasta_db
        self.index = index

    @property
    def extra(self) -> str:
        return self.fasta_db.extra(self.index)

    @property
    def identifier(self) -> str:
        return self.fasta_db.sequence_id(self.index)

    @property
    def sequence(self) -> str:
        return self.fasta_db.sequence(self.index)

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return f">{self.identifier} {self.extra}\n{self.sequence}"

    def __repr__(self):
        return f"FastaDbEntry(identifier={repr(self.identifier)}, sequence={repr(self.sequence)}, extra={repr(self.extra)})"


class FastaDb(ISample[FastaEntry], DbWrapper):
    """
    An LMDB-backed database of FASTA entries.
    """
    # __slots__ = ("length", "contains_sequence_id", "sequence_id_to_index", "entry", "_id_map", "_sequences")

    length: int

    _index_to_sequence_id: List[str]|None = None
    _sequence_id_to_index: Dict[str, int]|None = None
    _sequences: List[str]|None = None

    class InMemory(enum.Flag):
        Nothing = 0
        IdMaps = enum.auto()
        Sequences = enum.auto()
        Extra = enum.auto()
        All = IdMaps | Sequences | Extra

    def __init__(
        self,
        fasta_db_path: Union[str, Path],
        in_memory: InMemory = InMemory.Nothing,
    ):
        super().__init__(fasta_db_path)
        self.length = np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]

        if FastaDb.InMemory.IdMaps in in_memory:
            self._index_to_sequence_id = [self.db[f"{i}_id"].decode() for i in range(self.length)]
            self._sequence_id_to_index = {identifier: i for i, identifier in enumerate(self._index_to_sequence_id)}

        if FastaDb.InMemory.Sequences in in_memory:
            self._sequences = [self.db[f"{i}"].decode() for i in range(self.length)]

        if FastaDb.InMemory.Extra in in_memory:
            self._extra = [self.db[f"{i}_extra"].decode() for i in range(self.length)]

    def __len__(self):
        return self.length

    def contains_index(self, sequence_index: int) -> bool:
        return sequence_index < self.length

    def contains_sequence_id(self, sequence_id: str) -> bool:
        if self._sequence_id_to_index is not None:
            return sequence_id in self._sequence_id_to_index
        return f">{sequence_id}" in self.db

    def entry(self, sequence_index: int) -> FastaDbEntry:
        return FastaDbEntry(self, sequence_index)

    def extra(self, sequence_index: int) -> str:
        if self._extra is not None:
            return self._extra[sequence_index]
        return self.db[f"{sequence_index}_extra"].decode()

    def sequence(self, sequence_index: int) -> str:
        if self._sequences is not None:
            return self._sequences[sequence_index]
        return self.db[f"{sequence_index}"].decode()

    def sequence_id(self, sequence_index: int) -> str:
        if self._index_to_sequence_id is not None:
            return self._index_to_sequence_id[sequence_index]
        return self.db[f"{sequence_index}_id"].decode()

    def sequence_id_to_index(self, sequence_id: str) -> int:
        if self._sequence_id_to_index is not None:
            return self._sequence_id_to_index[sequence_id]
        return np.frombuffer(self.db[f">{sequence_id}"], dtype=np.int32, count=1)[0]

    @singledispatchmethod
    def __contains__(self, sequence_index: int) -> bool:
        return self.contains_index(sequence_index)

    @__contains__.register
    def _(self, sequence_id: str) -> bool:
        return self.contains_sequence_id(sequence_id)

    @__contains__.register
    def _(self, entry: IFastaEntry) -> bool:
        return self.contains_sequence_id(entry.identifier)

    def __iter__(self):
        for i in range(len(self)):
            yield self.entry(i)

    @singledispatchmethod
    def __getitem__(self, sequence_index: int) -> FastaDbEntry:
        return self.entry(sequence_index)

    @__getitem__.register
    def _(self, sequence_id: str) -> FastaDbEntry:
        return self.entry(self.sequence_id_to_index(sequence_id))

    def mappings(
        self,
        fasta_mapping_db_path: Union[str, Path],
        load_into_memory: bool = False
    ) -> "Tuple[FastaMappingEntry, ...]":
        return FastaMappingDb(fasta_mapping_db_path, self, load_into_memory).entries

    def sample(self, shape: Union[int, Tuple[int, ...]], rng: np.random.Generator) -> np.ndarray:
        """
        Sample sequences from the FASTA database.
        """
        result = np.empty(np.product(shape), dtype=object)
        result[:] = list(map(self.entry, rng.choice(self.length, len(result), replace=True)))
        return result.reshape(shape)


class FastaMappingEntryFactory:
    def __init__(self, name: str, fasta_mapping_db_factory: "FastaMappingDbFactory"):
        self.name = name
        self.fasta_mapping_db_factory = fasta_mapping_db_factory
        self.sequence_indices: list[int] = []

    def write_entry(self, entry: FastaEntry, abundance: int = 1):
        index = self.fasta_mapping_db_factory.fasta_db.sequence_id_to_index(entry.identifier)
        insert_index = bisect.bisect_right(self.sequence_indices, index)
        for i in range(abundance):
            self.sequence_indices.insert(insert_index+i, index)

    def write_entries(self, entries: Iterable[FastaEntry]):
        for entry in entries:
            self.write_entry(entry)

    def __enter__(self) -> "FastaMappingEntryFactory":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.fasta_mapping_db_factory.write_entry(self)


class FastaMappingEntry(ISample[FastaEntry]):

    _sequence_indices: npt.NDArray[np.uint32]|None

    def __init__(
        self,
        index: int,
        fasta_mapping_db: "FastaMappingDb",
        load_into_memory: bool = False,
    ):
        self.index = index
        self.fasta_mapping_db = fasta_mapping_db
        self.name = self.fasta_mapping_db.db[f"{index}_name"].decode()
        self.length = np.frombuffer(self.fasta_mapping_db.db[f"{index}_length"], dtype=np.int32, count=1)[0]

        if load_into_memory:
            self._sequence_indices = sequence_indices = np.empty(self.length, dtype=np.uint32)
            for i in range(self.length):
                self._sequence_indices[i] = np.frombuffer(self.fasta_mapping_db.db[f"{index}_{i}"], dtype=np.uint32, count=1)[0]
            self.sequence_index = lambda mapped_sequence_index: sequence_indices[mapped_sequence_index]
        else:
            self._sequence_indices = None
            self.sequence_index = lambda mapped_sequence_index: np.frombuffer(self.fasta_mapping_db.db[f"{index}_{mapped_sequence_index}"], dtype=np.uint32, count=1)[0]

    def entry(self, mapped_sequence_index: int) -> FastaDbEntry:
        return self.fasta_mapping_db.fasta_db.entry(self.sequence_index(mapped_sequence_index))

    @singledispatchmethod
    def __contains__(self, sequence_index: int) -> bool:
        if self._sequence_indices is not None:
            return sequence_index == np.searchsorted(self._sequence_indices, sequence_index)
        low = 0
        high = self.length - 1
        while low <= high:
            mid = (low + high) // 2
            mid_val = self.sequence_index(mid)
            if mid_val < sequence_index:
                low = mid + 1
            elif mid_val > sequence_index:
                high = mid - 1
            else:
                return True
        return False

    @__contains__.register
    def _(self, sequence_id: str) -> bool:
        sequence_index = self.fasta_mapping_db.fasta_db.sequence_id_to_index(sequence_id)
        return sequence_index in self

    @__contains__.register
    def _(self, entry: FastaEntry) -> bool:
        return entry.identifier in self

    @singledispatchmethod
    def __getitem__(self, mapped_sequence_index: int) -> FastaDbEntry:
        return self.entry(mapped_sequence_index)

    @__getitem__.register
    def _(self, sequence_id: str) -> FastaDbEntry:
        return self.fasta_mapping_db.fasta_db[sequence_id]

    def __iter__(self) -> Generator[FastaDbEntry, None, None]:
        for i in range(self.length):
            yield self.entry(i)

    def __len__(self) -> int:
        return self.length

    def sample(self, shape: Union[int, Tuple[int, ...]], rng: np.random.Generator) -> np.ndarray:
        """
        Sample sequences from the FASTA database.
        """
        result = np.empty(np.product(shape), dtype=object)
        result[:] = list(map(self.entry, rng.choice(self.length, len(result), replace=True)))
        return result.reshape(shape)


class FastaMappingDbFactory(DbFactory):
    def __init__(self, path: Union[str, Path], fasta_db: FastaDb, chunk_size: int = 10000):
        super().__init__(path, chunk_size)
        self.fasta_db = fasta_db
        self.num_entries = np.int32(0)
        self.write("fasta_db_uuid", self.fasta_db.uuid.bytes)

    def create_entry(self, name: str) -> FastaMappingEntryFactory:
        return FastaMappingEntryFactory(name, self)

    def write_entry(self, entry: FastaMappingEntryFactory):
        self.write(f"{self.num_entries}_name", entry.name.encode())
        self.write(f"{self.num_entries}_length", np.int32(len(entry.sequence_indices)).tobytes())
        for i, index in enumerate(np.array(entry.sequence_indices, dtype=np.uint32)):
            self.write(f"{self.num_entries}_{i}", index.tobytes())
        self.num_entries += 1

    def write_entries(self, entries: Iterable[FastaMappingEntryFactory]):
        for entry in entries:
            self.write_entry(entry)

    def before_close(self):
        self.write("length", self.num_entries.tobytes())
        super().before_close()


class FastaMappingDb(DbWrapper):
    def __init__(
        self,
        path: Union[str, Path],
        fasta_db: FastaDb,
        load_into_memory: bool = False
    ):
        super().__init__(path)
        assert fasta_db.uuid.bytes == self.db["fasta_db_uuid"], "This FASTA Mapping was not created with the given FASTA DB."
        self.fasta_db = fasta_db
        self.length = np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]
        self.entries = tuple(FastaMappingEntry(i, self, load_into_memory) for i in range(self.length))

    def __getitem__(self, index: int) -> FastaMappingEntry:
        return self.entries[index]

    def __len__(self):
        return self.length

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


# def entries_with_taxonomy(
#     sequences: Iterable[FastaEntry],
#     taxonomies: Iterable[TaxonomyEntry],
# ) -> Generator[Tuple[FastaEntry, TaxonomyEntry], None, None]:
#     """
#     Efficiently iterate over a FASTA file with a corresponding taxonomy file
#     """
#     labels = {}
#     taxonomy_iterator = iter(taxonomies)
#     taxonomy: TaxonomyEntry
#     for sequence in sequences:
#         while sequence.identifier not in labels:
#             taxonomy = next(taxonomy_iterator)
#             labels[taxonomy.identifier] = taxonomy
#         taxonomy = labels[sequence.identifier]
#         del labels[sequence.identifier]
#         yield sequence, taxonomy


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
