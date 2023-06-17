import abc
from dataclasses import dataclass, field
import numpy as np
import numpy.typing as npt
from pathlib import Path
from typing import Generator, Iterable, Optional, Tuple, Union

from .db import DbFactory, DbWrapper
from .dna import AbstractSequenceWrapper
from .fasta import FastaDb, FastaEntry, FastaIndexDb
from .fastq import FastqDb, FastqEntry
from .types import int_t

@dataclass(frozen=True, order=True)
class SampleMappingEntry:
    name: str
    indices: npt.NDArray[np.int64] = field(compare=False, hash=False)
    abundances: npt.NDArray[np.int64] = field(compare=False, hash=False)
    fasta_index_db: FastaIndexDb = field(compare=False, hash=False, repr=False)

    @classmethod
    def deserialize(cls, entry: bytes, fasta_index_db: FastaIndexDb):
        name, values = entry.split(b'\x00', maxsplit=1)
        values = np.frombuffer(values, dtype=np.int64)
        name = name.decode()
        indices = values[:len(values)//2]
        abundances = values[len(values)//2:]
        return cls(name, indices, abundances, fasta_index_db)

    def contains_fasta_id(self, fasta_id: str) -> bool:
        index = self.fasta_index_db.fasta_id_to_index(fasta_id)
        return self.contains_index(index)

    def contains_key(self, key: str) -> bool:
        index = self.fasta_index_db.key_to_index(key)
        return self.contains_index(index)

    def contains_index(self, index: int_t) -> bool:
        return index in self.indices

    def serialize(self) -> bytes:
        return self.name.encode() + b'\x00' + self.indices.tobytes() + self.abundances.tobytes()

    @property
    def total_abundance(self) -> int:
        return int(np.sum(self.abundances))

    def __getitem__(self, index: int_t) -> str:
        return self.fasta_index_db.index_to_fasta_id(self.indices[index])

    def __len__(self):
        return len(self.indices)

    def __iter__(self) -> Generator[Tuple[str, int], None, None]:
        for index, abundance in zip(self.indices, self.abundances):
            yield (self.fasta_index_db.index_to_fasta_id(index), abundance)

    def __eq__(self, other: 'SampleMappingEntry'):
        return self.name == other.name \
            and np.array_equal(self.indices, other.indices) \
            and np.array_equal(self.abundances, other.abundances)


class SampleMappingEntryFactory:
    def __init__(self, name: str, fasta_index_db: FastaIndexDb):
        self.name = name
        self.fasta_index_db = fasta_index_db
        self.abundance_map: dict[int, int] = {}

    def add_fasta_id(self, fasta_id: str, abundance: int = 1) -> 'SampleMappingEntryFactory':
        index = self.fasta_index_db.fasta_id_to_index(fasta_id)
        if index not in self.abundance_map:
            self.abundance_map[index] = 0
        self.abundance_map[index] += abundance
        return self

    def add_entry(self, entry: FastaEntry, abundance: int = 1) -> 'SampleMappingEntryFactory':
        return self.add_fasta_id(entry.identifier, abundance)

    def add_entries(self, entries: Iterable[FastaEntry]) -> 'SampleMappingEntryFactory':
        for entry in entries:
            self.add_entry(entry)
        return self

    def build(self) -> SampleMappingEntry:
        indices = np.array(sorted(self.abundance_map.keys()))
        abundances = np.array([self.abundance_map[i] for i in indices])
        return SampleMappingEntry(self.name, indices, abundances, self.fasta_index_db)


class SampleMappingDbFactory(DbFactory):
    """

    """
    def __init__(self, path: Union[str, Path], chunk_size: int_t = 10000):
        super().__init__(path, chunk_size)
        self.num_entries = np.int32(0)

    def write_entry(self, entry: SampleMappingEntry):
        """
        Create a new FASTA LMDB database from a FASTA file.
        """
        self.write(f"name_{entry.name}", np.int32(self.num_entries).tobytes())
        self.write(str(self.num_entries), entry.serialize())
        self.num_entries += 1

    def write_entries(self, entries: Iterable[SampleMappingEntry]):
        for entry in entries:
            self.write_entry(entry)

    def before_close(self):
        self.write("length", self.num_entries.tobytes())
        super().before_close()


class SampleMappingDb(DbWrapper):
    def __init__(self, path: Union[str, Path], fasta_index_db: FastaIndexDb):
        super().__init__(path)
        self.fasta_index_db = fasta_index_db
        self.length = np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]

    def __getitem__(self, index_or_name: Union[int, str]) -> SampleMappingEntry:
        if isinstance(index_or_name, str):
            index_or_name = np.frombuffer(
                self.db[f"name_{index_or_name}"],
                dtype=np.int32,
                count=1)[0]
        return SampleMappingEntry.deserialize(self.db[str(index_or_name)], self.fasta_index_db)

    def __len__(self):
        return self.length

    def __iter__(self) -> Generator[SampleMappingEntry, None, None]:
        for i in range(len(self)):
            yield self[i]


class SampleInterface(abc.ABC):
    def __init__(self, name: str):
        self.name = name

    @abc.abstractmethod
    def sample(
        self,
        n: int,
        rng: np.random.Generator = np.random.default_rng()
    ) -> Generator[AbstractSequenceWrapper, None, None]:
        raise NotImplementedError()

    @abc.abstractmethod
    def __contains__(self, key) -> bool:
        raise NotImplementedError()

    @abc.abstractmethod
    def __len__(self) -> bool:
        raise NotImplementedError()

    @abc.abstractmethod
    def __getitem__(self, key):
        raise NotImplementedError()


class FastaSample(SampleInterface):
    def __init__(self, fasta_db: FastaDb, name: Optional[str] = None):
        super().__init__(name or fasta_db.path.name)
        self.fasta_db = fasta_db

    def sample(
        self,
        n: int_t,
        replace: bool = True,
        rng: np.random.Generator = np.random.default_rng()
    ) -> Generator[FastaEntry, None, None]:
        indices = rng.choice(len(self.fasta_db), size=n, replace=replace)
        indices, counts = np.unique(indices, return_counts=True)
        for index, count in zip(indices, counts):
            fasta_entry = self.fasta_db[index]
            for _ in range(count):
                yield fasta_entry

    def __contains__(self, index_or_id: Union[int, str]) -> bool:
        return index_or_id in self.fasta_db

    def __len__(self) -> int:
        return len(self.fasta_db)

    def __getitem__(self, index_or_id: Union[int, str]) -> FastaEntry:
        return self.fasta_db[index_or_id]

    def __iter__(self) -> Generator[FastaEntry, None, None]:
        return iter(self.fasta_db)


class FastqSample(SampleInterface):
    def __init__(self, fastq_db: FastqDb, name: Optional[str] = None):
        super().__init__(name or fastq_db.path.name)
        self.fastq_db = fastq_db

    def sample(
        self,
        n: int_t,
        replace: bool = True,
        rng: np.random.Generator = np.random.default_rng()
    ) -> Generator[FastqEntry, None, None]:
        indices = rng.choice(len(self.fastq_db), size=n, replace=replace)
        indices, counts = np.unique(indices, return_counts=True)
        for index, count in zip(indices, counts):
            fastq_entry = self.fastq_db[index]
            for _ in range(count):
                yield fastq_entry

    def __contains__(self, index: int_t) -> bool:
        return index in self.fastq_db

    def __len__(self) -> int:
        return len(self.fastq_db)

    def __getitem__(self, index: int_t) -> FastqEntry:
        return self.fastq_db[index]

    def __iter__(self):
        return iter(self.fastq_db)


class DemultiplexedFastaSample(FastaSample):
    def __init__(self, fasta_db: FastaDb, sample_mapping: SampleMappingEntry):
        super().__init__(fasta_db, sample_mapping.name)
        self.sample_mapping = sample_mapping
        self.cumulative_abundance = np.cumsum(self.sample_mapping.abundances)

    def sample(
        self,
        n: int_t,
        replace: bool = True,
        rng: np.random.Generator = np.random.default_rng()
    ) -> Generator[FastaEntry, None, None]:
        """
        Sample random entries from the FASTA database according to the abundance distribution.
        """
        p = self.sample_mapping.abundances / self.sample_mapping.abundances.sum()
        indices = rng.choice(len(self.sample_mapping), size=n, replace=replace, p=p)
        indices, counts = np.unique(indices, return_counts=True)
        for index, count in zip(indices, counts):
            fasta_entry = self.fasta_db[self.sample_mapping[index]]
            for _ in range(count):
                yield fasta_entry

    @property
    def abundances(self):
        return self.sample_mapping.abundances

    @property
    def total_abundance(self):
        return self.sample_mapping.total_abundance

    def __getitem__(self, index_or_id: Union[int, str]) -> FastaEntry:
        if isinstance(index_or_id, str):
            if not self.sample_mapping.contains_fasta_id(index_or_id):
                raise KeyError(f"FASTA ID {index_or_id} not found in sample mapping.")
            return self.fasta_db[index_or_id]
        index = np.searchsorted(self.cumulative_abundance, index_or_id+1, side='left')
        return self.fasta_db[self.sample_mapping[index]]

    def __len__(self):
        return self.cumulative_abundance[-1]

    def __iter__(self) -> Generator[FastaEntry, None, None]:
        for fasta_id, abundance in self.sample_mapping:
            entry = self.fasta_db[fasta_id]
            for _ in range(abundance):
                yield entry

# Utility functions for loading/wrapping FASTA/FASTQ databases with sample interfaces.

def load_fasta(
    fasta_db_or_path: Union[FastaDb, Union[str, Path]],
    name: Optional[str] = None
) -> FastaSample:
    if not isinstance(fasta_db_or_path, FastaDb):
        fasta_db_or_path = FastaDb(fasta_db_or_path)
    return FastaSample(fasta_db_or_path, name)


def load_fastq(
    fastq_db_or_path: Union[FastqDb, Union[str, Path]],
    name: Optional[str] = None
) -> FastqSample:
    if not isinstance(fastq_db_or_path, FastqDb):
        fastq_db_or_path = FastqDb(fastq_db_or_path)
    return FastqSample(fastq_db_or_path, name)


def load_multiplexed_fasta(
    fasta_db_or_path: Union[FastaDb, Union[str, Path]],
    sample_mapping_db_or_path: Union[SampleMappingDb, Union[str, Path]],
    fasta_index_db_or_path: Optional[Union[FastaIndexDb, Union[str, Path]]] = None
) -> Tuple[DemultiplexedFastaSample, ...]:
    if not isinstance(fasta_db_or_path, FastaDb):
        fasta_db_or_path = FastaDb(fasta_db_or_path)
    if not isinstance(sample_mapping_db_or_path, SampleMappingDb):
        if fasta_index_db_or_path is None:
            assert str(sample_mapping_db_or_path).endswith(".mapping.db"), \
                "Unable to automatically resolve the index path."
            fasta_index_db_or_path = str(sample_mapping_db_or_path)[:-11] + ".index.db"
        if not isinstance(fasta_index_db_or_path, FastaIndexDb):
            fasta_index_db_or_path = FastaIndexDb(fasta_index_db_or_path)
        sample_mapping_db_or_path = SampleMappingDb(
            sample_mapping_db_or_path,
            fasta_index_db_or_path)
    sample_entries = sorted(iter(sample_mapping_db_or_path))
    return tuple(DemultiplexedFastaSample(fasta_db_or_path, entry) for entry in sample_entries)
