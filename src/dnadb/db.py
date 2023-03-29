from lmdbm import Lmdb
from pathlib import Path

class DbFactory:
    """
    A factory for creating LMDB-backed databases of FASTA entries.
    """
    def __init__(self, path: str|Path, chunk_size: int = 10000):
        self.path = Path(path)
        self.db = Lmdb.open(str(path), "n")
        self.buffer: dict[str|bytes, bytes] = {}
        self.chunk_size = chunk_size

    def flush(self):
        self.db.update(self.buffer)
        self.buffer.clear()

    def write(self, key: str|bytes, value: bytes):
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
