from lmdbm import Lmdb
import numpy as np
import numpy.typing as npt
from pathlib import Path
from typing import Generator, Iterable, TextIO, Tuple, Union

from .db import DbFactory
from .utils import open_file

def phred_encode(probabilities: npt.ArrayLike, encoding: int = 33) -> str:
    scores = (-10 * np.log10(np.array(probabilities))).astype(int)
    return ''.join((chr(score + encoding)) for score in scores)


def phred_decode(qualities: str, encoding: int = 33) -> npt.NDArray[np.float64]:
    scores = np.array([(ord(token) - encoding) for token in qualities])
    return 10**(scores / -10)


class FastqHeader:
    """
    A class representation of the sequence identifier of a FASTQ entry.
    """
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
            pos=tuple(map(int, left[5:])),
            read_type=int(right[0]),
            is_filtered=right[1] == 'Y',
            control_number=int(right[2]),
            sequence_index=right[3]
        )

    def __init__(
        self,
        instrument: str,
        run_number: int,
        flowcell_id: str,
        lane: int,
        tile: int,
        pos: Tuple[int, int],
        read_type: int,
        is_filtered: bool,
        control_number: int,
        sequence_index: str
    ):
        self.instrument = instrument
        self.run_number = run_number
        self.flowcell_id = flowcell_id
        self.lane = lane
        self.tile = tile
        self.pos = pos
        self.read_type = read_type
        self.is_filtered = is_filtered
        self.control_number = control_number
        self.sequence_index = sequence_index

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

    def __repr__(self):
        return str(self)


class FastqEntry:
    """
    A class representation of a FASTQ entry containing the sequnce identifier, sequence, and quality
    scores.
    """
    @classmethod
    def deserialize(cls, entry: bytes) -> "FastqEntry":
        return cls.from_str(entry.decode())

    @classmethod
    def from_str(cls, entry: str) -> "FastqEntry":
        header, sequence, _, quality_scores= entry.rstrip().split('\n')
        return FastqEntry(FastqHeader.from_str(header), sequence, quality_scores)

    def __init__(
        self,
        header: FastqHeader,
        sequence: str,
        quality_scores: str,
    ):
        self.header = header
        self.sequence = sequence.strip()
        self.quality_scores = quality_scores.strip()

    def serialize(self) -> bytes:
        return str(self).encode()

    def __str__(self):
        return f"{self.header}\n{self.sequence}\n+\n{self.quality_scores}"

    def __repr__(self):
        return "FastqEntry:\n" + '\n'.join(f"  {s}" for s in str(self).split('\n'))


class FastqDbFactory(DbFactory):
    """
    A factory for creating LMDB-backed databases of FASTA entries.
    """
    def __init__(self, path: Union[str, Path], chunk_size: int = 10000):
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

    def close(self):
        self.write("length", self.num_entries.tobytes())
        super().close()


class FastqDb:
    def __init__(self, fastq_db_path: Union[str, Path]):
        self.path = Path(fastq_db_path).absolute
        self.db = Lmdb.open(str(fastq_db_path), lock=False)
        self.length = np.frombuffer(self.db["length"], dtype=np.int32, count=1)[0]

    def __len__(self):
        return self.length

    def __getitem__(self, sequence_index: int) -> FastqEntry:
        return FastqEntry.deserialize(self.db[str(sequence_index)])


def entries(sequences: Union[TextIO, Iterable[FastqEntry], str, Path]) -> Iterable[FastqEntry]:
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


def read(buffer: TextIO) -> Generator[FastqEntry, None, None]:
    """
    Read entries from a FASTQ file buffer.
    """
    line = buffer.readline()
    while len(line) > 0:
        entry_str = "".join((line, *(buffer.readline() for _ in range(3))))
        yield FastqEntry.from_str(entry_str)
        line = buffer.readline()


def write(buffer: TextIO, entries: Iterable[FastqEntry]) -> int:
    """
    Write entries to a FASTQ file.
    """
    bytes_written = 0
    for entry in entries:
        bytes_written += buffer.write(str(entry) + '\n')
    return bytes_written
