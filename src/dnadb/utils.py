import io
from pathlib import Path
import requests
from skbio.io.util import open as skbio_open
import subprocess
from tqdm.auto import tqdm
from typing import cast, Union


def download(url: str, destination: Union[str, Path], chunk_size: int = 1024):
    """
    Download a file from the internet to the provided destination.
    """
    # https://gist.github.com/yanqd0/c13ed29e29432e3cf3e7c38467f42f51
    response = requests.get(url, stream=True)
    if response.status_code != 200:
        raise Exception(f"Could not download: {url}\n{response.content}")
    total = int(response.headers.get('content-length', 0))
    with open(str(destination), 'wb') as file, tqdm(
        desc=f"Downloading: {destination}",
        total=total,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
        leave=False
    ) as bar:
        for data in response.iter_content(chunk_size=chunk_size):
            size = file.write(data)
            bar.update(size)


def compress(path: Union[str, Path]):
    """
    Compress the given file with gzip.
    """
    subprocess.run(["gzip", "-f", str(path)])


def decompress(path: Union[str, Path]):
    """
    Decompress the given file with gunzip.
    """
    subprocess.run(["gunzip", "-f", str(path)])


def open_file(path: Union[str, Path], mode: str = "r") -> io.TextIOWrapper:
    """
    Open a file wihtout worrying about compression
    """
    return cast(io.TextIOWrapper, skbio_open(str(path), mode))
