from lmdbm import Lmdb
from pathlib import Path
from typing import Union

class DbFactory:
    """
    A factory for creating LMDB-backed databases of FASTA entries.
    """
    def __init__(self, path: Union[str, Path], chunk_size: int = 10000):
        self.path = Path(path)
        self.db = Lmdb.open(str(path), "n", lock=True)
        self.buffer: dict[Union[str, bytes], bytes] = {}
        self.chunk_size = chunk_size

    def flush(self):
        self.db.update(self.buffer)
        self.buffer.clear()

    def write(self, key: Union[str, bytes], value: bytes):
        self.buffer[key] = value
        if len(self.buffer) >= self.chunk_size:
            self.flush()

    def close(self):
        self.flush()
        self.db.close()

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()
