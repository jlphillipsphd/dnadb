from dataclasses import dataclass
import io
import numpy as np
import numpy.typing as npt
from pathlib import Path
from typing import Generator, Iterable, Tuple, Union

from .db import DbFactory, DbWrapper
from .dna import AbstractSequenceWrapper
from .sample import ISample
from .types import int_t
from .utils import open_file

def phred_encode(probabilities: npt.ArrayLike, encoding: int_t = 33) -> str:
    scores = (-10 * np.log10(np.array(probabilities))).astype(int)
    return ''.join((chr(score + encoding)) for score in scores)


def phred_decode(qualities: str, encoding: int_t = 33) -> npt.NDArray[np.float64]:
    scores = np.array([(ord(token) - encoding) for token in qualities])
    return 10**(scores / -10)


@dataclass(frozen=True, order=True)
class FastqHeader:
    """
    A class representation of the header of a FASTQ entry.
    """
    __slots__ = (
        "instrument",
        "run_number",
        "flowcell_id",
        "lane",
        "tile",
        "pos",
        "read_type",
        "is_filtered",
        "control_number",
        "sequence_index"
    )

    instrument: str
    run_number: int_t
    flowcell_id: str
    lane: int_t
    tile: int_t
    pos: Tuple[int, int]
    read_type: int_t
    is_filtered: bool
    control_number: int_t
    sequence_index: str

    @classmethod
    def deserialize(cls, sequence_id: bytes):
        return cls.from_str(sequence_id.decode())

    @classmethod
    def from_str(cls, sequence_id: str) -> "FastqHeader":
        # Split up the sequence ID information
        left, right = sequence_id.strip()[1:].split(' ')
        left = left.split(':')
        right = right.split(':')
        return cls(
            instrument=left[0],
            run_number=int(left[1]),
            flowcell_id=left[2],
            lane=int(left[3]),
            tile=int(left[4]),
            pos=(int(left[5]), int(left[6])),
            read_type=int(right[0]),
            is_filtered=right[1] == 'Y',
            control_number=int(right[2]),
            sequence_index=right[3]
        )

    # Serialize a FastqHeader object to a byte string
    def serialize(self) -> bytes:
        return str(self).encode()

    def __str__(self):
        sequence_id = '@'
        sequence_id += ':'.join(map(str, [
            self.instrument,
            self.run_number,
            self.flowcell_id,
            self.lane,
            self.tile,
            *self.pos
        ]))
        sequence_id += ' '
        sequence_id += ':'.join(map(str, [
            self.read_type,
            'Y' if self.is_filtered else 'N',
            self.control_number,
            self.sequence_index
        ]))
        return sequence_id


@dataclass(frozen=True, order=True)
class FastqEntry(AbstractSequenceWrapper):
    """
    A class representation of a FASTQ entry containing the sequnce identifier, sequence, and quality
    scores.
    """
    __slots__ = ("header_str", "quality_scores")

    header_str: str
    quality_scores: str

    @classmethod
    def deserialize(cls, entry: bytes) -> "FastqEntry":
        return cls(*entry.decode().split('\x00'))

    @classmethod
    def from_str(cls, entry: str) -> "FastqEntry":
        header, sequence, _, quality_scores= entry.rstrip().split('\n')
        return cls(header, sequence, quality_scores)

    def __init__(self, header: str, sequence: str, quality_scores: str):
        object.__setattr__(self, "header_str", header)
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "quality_scores", quality_scores)

    def serialize(self) -> bytes:
        return '\x00'.join((self.header_str, self.sequence, self.quality_scores)).encode()

    @property
    def header(self):
        return FastqHeader.from_str(self.header_str)

    def __str__(self):
        return f"{self.header_str}\n{self.sequence}\n+\n{self.quality_scores}"


class FastqDbFactory(DbFactory):
    """
    A factory for creating LMDB-backed databases of FASTA entries.
    """
    def __init__(self, path: Union[str, Path], chunk_size: int_t = 10000):
        super().__init__(path, chunk_size)
        self.num_entries = np.int32(0)

    def write_entry(self, entry: FastqEntry):
        """
        Create a new FASTA LMDB database from a FASTA file.
        """
        self.write(str(self.num_entries), entry.serialize())
        self.num_entries += 1

    def write_entries(self, entries: Iterable[FastqEntry]):
        for entry in entries:
            self.write_entry(entry)

    def before_close(self):
        self.write("length", self.num_entries.tobytes())
        super().before_close()


class FastqDb(ISample[FastqEntry], DbWrapper):
    def __init__(self, fastq_db_path: Union[str, Path]):
        super().__init__(fastq_db_path)
        self.length = np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]

    def __len__(self):
        return self.length

    def __contains__(self, sequence_index: int_t) -> bool:
        return sequence_index < self.length

    def __iter__(self) -> Generator[FastqEntry, None, None]:
        for i in range(len(self)):
            yield self[i]

    def __getitem__(self, sequence_index: int_t) -> FastqEntry:
        return FastqEntry.deserialize(self.db[str(sequence_index)])

    def sample(self, shape: Union[int, Tuple[int, ...]], rng: np.random.Generator) -> np.ndarray:
        """
        Sample sequences from the FASTA database.
        """
        result = np.empty(np.product(shape), dtype=object)
        result[:] = list(map(self.__getitem__, rng.choice(self.length, len(result), replace=True)))
        return result.reshape(shape)


def entries(
    sequences: Union[io.TextIOBase, Iterable[FastqEntry], str, Path]
) -> Iterable[FastqEntry]:
    """
    Create an iterator over a FASTQ file or iterable of FASTQ entries.
    """
    if isinstance(sequences, (str, Path)):
        with open_file(sequences, 'r') as buffer:
            yield from read(buffer)
    elif isinstance(sequences, io.TextIOBase):
        yield from read(sequences)
    else:
        yield from sequences


def read(buffer: io.TextIOBase) -> Generator[FastqEntry, None, None]:
    """
    Read entries from a FASTQ file buffer.
    """
    line = buffer.readline()
    while len(line) > 0:
        entry_str = "".join((line, *(buffer.readline() for _ in range(3))))
        yield FastqEntry.from_str(entry_str)
        line = buffer.readline()


def write(buffer: io.TextIOBase, entries: Iterable[FastqEntry]) -> int:
    """
    Write entries to a FASTQ file.
    """
    bytes_written = 0
    for entry in entries:
        bytes_written += buffer.write(str(entry) + '\n')
    return bytes_written
