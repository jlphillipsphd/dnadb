from lmdbm import Lmdb
from pathlib import Path
from typing import TypeVar, Union

from .types import int_t

T = TypeVar("T")

class DbFactory:
    """
    A factory for creating LMDB-backed databases of FASTA entries.
    """
    def __init__(self, path: Union[str, Path], chunk_size: int_t = 10000):
        self.path = Path(path)
        if self.path.suffix != ".db":
            self.path = Path(str(self.path) + ".db")
        self.db = Lmdb.open(str(self.path), "n", lock=True)
        self.buffer: dict[Union[str, bytes], bytes] = {}
        self.chunk_size = chunk_size
        self.is_closed = False

    def flush(self):
        self.db.update(self.buffer)
        self.buffer.clear()

    def contains(self, key: Union[str, bytes]) -> bool:
        return key in self.buffer or key in self.db

    def read(self, key: Union[str, bytes]) -> bytes:
        return self.buffer[key] if key in self.buffer else self.db[key]

    def append(self, key: Union[str, bytes], value: bytes):
        self.write(key, self.read(key) + value)

    def write(self, key: Union[str, bytes], value: bytes):
        self.buffer[key] = value
        if len(self.buffer) >= self.chunk_size:
            self.flush()

    def before_close(self):
        self.flush()

    def close(self):
        if self.is_closed:
            return
        self.before_close()
        self.db.close()
        self.is_closed = True

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()


class DbWrapper:
    __slots__ = ("__path", "__db", "__is_closed")

    def __init__(self, path: Union[str, Path]):
        self.__path = Path(path).absolute()
        self.__db = Lmdb.open(str(path), lock=False)
        self.__is_closed = False

    def close(self):
        if self.__is_closed:
            return
        self.__is_closed = True
        return self.__db.close()

    @property
    def db(self) -> Lmdb:
        return self.__db

    @property
    def path(self) -> Path:
        return self.__path

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()
