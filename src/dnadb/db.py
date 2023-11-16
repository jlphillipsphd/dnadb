from lmdbm import Lmdb
from pathlib import Path
from typing import TypeVar, Union
import uuid

from .types import int_t

T = TypeVar("T")

class DbFactory:
    """
    A wrapper around lmdbm.Lmdb that allows for buffered writes.
    """
    def __init__(self, path: Union[str, Path], chunk_size: int_t = 10000):
        self.path = Path(path)
        if self.path.suffix != ".db":
            self.path = Path(str(self.path) + ".db")
        self.db = Lmdb.open(str(self.path), "n", lock=True)
        self.buffer: dict[Union[str, bytes], bytes] = {}
        self.chunk_size = chunk_size
        self.is_closed = False
        self.uuid = uuid.uuid4()
        self.write("uuid", self.uuid.bytes)

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
        self.is_closed = True
        self.before_close()
        self.db.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()


class DbWrapper:
    __slots__ = ("_path", "_db", "_is_closed", "_uuid")

    def __init__(self, path: Union[str, Path]):
        self._path = Path(path).absolute()
        self._db = Lmdb.open(str(path), lock=False)
        self._is_closed = False
        try:
            self._uuid = uuid.UUID(bytes=self.db["uuid"])
        except KeyError:
            raise Exception(f"Database at {self.path} does not have a UUID.")

    def close(self):
        if self._is_closed:
            return
        self._is_closed = True
        return self._db.close()

    @property
    def db(self) -> Lmdb:
        return self._db

    @property
    def path(self) -> Path:
        return self._path

    @property
    def uuid(self) -> uuid.UUID:
        return self._uuid

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()
