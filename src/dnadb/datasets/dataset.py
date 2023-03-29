import abc
import os
from pathlib import Path
from tqdm.auto import tqdm
from typing import Generator, List, Optional, Tuple, Union

from .. import fasta, fastq

from .. import taxonomy
from ..utils import open_file, download

# Utility Functions --------------------------------------------------------------------------------

def __has_methods(subclass: object, methods: List[str]):
    """
    Determine if a class has all the methods in the list.
    """
    for method in methods:
        if not (hasattr(subclass, method) and callable(getattr(subclass, method))):
            return False
    return True

# Interface Classes --------------------------------------------------------------------------------

class IPathProvider(metaclass=abc.ABCMeta):
    @classmethod
    def __subclasshook__(cls, subclass):
        return __has_methods(subclass, ["path"]) or NotImplemented

    @property
    @abc.abstractmethod
    def path(self) -> Path:
        raise NotImplementedError()

# Dataset Classes ----------------------------------------------------------------------------------

class Dataset(abc.ABC, IPathProvider):

    NAME: str

    def __init__(self, path: Optional[Union[str, Path]] = None):
        if path is None:
            self.__path = Path(os.environ["DATASETS_PATH"]) / self.name
        else:
            self.__path = Path(path)

    @property
    @abc.abstractmethod
    def name(self) -> str:
        raise NotImplementedError()

    @property
    def path(self) -> Path:
        return self.__path


class VersionedDataset(Dataset):
    def __init__(self, version: str, path: Optional[Union[str, Path]] = None):
        super().__init__(path=path)
        self.version = version

    @property
    def version(self) -> str:
        return self.__version

    @version.setter
    def version(self, version: str):
        self.__version = version.replace('_', '.')

# Dataset Mixins -----------------------------------------------------------------------------------

# For more about mixin architectures: https://stackoverflow.com/a/50465583
class RemoteDatasetMixin(IPathProvider, metaclass=abc.ABCMeta):
    @classmethod
    def __subclasshook__(cls, subclass):
        return __has_methods(subclass, ["fetch", "files", "on_dataset_downloaded"]) \
            or NotImplemented

    def __init__(self, *args, force_download=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.path.mkdir(exist_ok=True)
        self.__download_files(force_download)

    def __download_files(self, force_download: bool = True):
        """
        Fetch the dataset from the remote source.
        """
        to_download = []
        for filename in self.remote_files:
            if not (self.path / filename).exists() or force_download:
                to_download.append(filename)
        if len(to_download) == 0:
            return
        base_url = self.url.rstrip('/')
        for filename in tqdm(tuple(map(str, to_download)), desc="Downloading dataset files"):
            filename = filename.lstrip('/')
            url = f"{base_url}/{filename}"
            path = self.path / filename
            path.parent.mkdir(parents=True, exist_ok=True)
            download(url, path)
        self.on_dataset_downloaded(force_download)

    @property
    @abc.abstractmethod
    def remote_files(self) -> List[Path]:
        """
        The list of remote files for this dataset.
        """
        raise NotImplementedError()

    def on_dataset_downloaded(self, force_download: bool):
        """
        Invoked when the dataset has finished downloading.
        """
        pass

    @property
    @abc.abstractmethod
    def url(self) -> str:
        """
        The URL prefix to the dataset files.
        """
        raise NotImplementedError()

# Dataset Traits -----------------------------------------------------------------------------------

class InterfacesWithFasta(abc.ABC, IPathProvider):
    @classmethod
    def __subclasshook__(cls, subclass):
        return __has_methods(subclass, ["sequences", "fasta_file"]) or NotImplemented

    def sequences(self) -> Generator[fasta.FastaEntry, None, None]:
        yield from fasta.read(open_file(self.path / self.fasta_file))

    @property
    @abc.abstractmethod
    def fasta_file(self) -> Union[str, Path]:
        raise NotImplementedError()


class InterfacesWithFastq(abc.ABC, IPathProvider):
    @classmethod
    def __subclasshook__(cls, subclass):
        return __has_methods(subclass, ["sequences", "fastq_file"]) or NotImplemented

    def sequences(self) -> Generator[fastq.FastqEntry, None, None]:
        yield from fastq.read(open_file(self.path / self.fastq_file))

    @property
    @abc.abstractmethod
    def fastq_file(self) -> Union[str, Path]:
        raise NotImplementedError()


class InterfacesWithTaxonomy(abc.ABC, IPathProvider):
    @classmethod
    def __subclasshook__(cls, subclass):
        return __has_methods(subclass, ["sequences", "taxonomy_file"]) or NotImplemented

    def taxonomies(self) -> Generator[taxonomy.TaxonomyEntry, None, None]:
        yield from taxonomy.read(open_file(self.path / self.taxonomy_file))

    @property
    @abc.abstractmethod
    def taxonomy_file(self) -> Union[str, Path]:
        raise NotImplementedError()


class InterfacesWithFastas(abc.ABC, IPathProvider):
    @classmethod
    def __subclasshook__(cls, subclass):
        return __has_methods(subclass, ["sequences", "fasta_files"]) or NotImplemented

    def sequences(self) -> Generator[Tuple[str, fasta.FastaEntry], None, None]:
        for path in self.fasta_files:
            sample_file_name = Path(path).name
            with open_file(self.path / path) as file:
                for sequence_entry in fasta.read(file):
                    yield sample_file_name, sequence_entry

    @property
    @abc.abstractmethod
    def fasta_files(self) -> List[Union[str, Path]]:
        raise NotImplementedError()


class InterfacesWithFastqs(abc.ABC, IPathProvider):
    @classmethod
    def __subclasshook__(cls, subclass):
        return __has_methods(subclass, ["sequences", "fastq_files"]) or NotImplemented

    def sequences(self) -> Generator[Tuple[str, fastq.FastqEntry], None, None]:
        for path in self.fastq_files:
            sample_file_name = Path(path).name
            with open_file(self.path / path) as file:
                for sequence_entry in fastq.read(file):
                    yield sample_file_name, sequence_entry

    @property
    @abc.abstractmethod
    def fastq_files(self) -> List[Union[str, Path]]:
        raise NotImplementedError()


class InterfacesWithTaxonomies(abc.ABC, IPathProvider):
    @classmethod
    def __subclasshook__(cls, subclass):
        return __has_methods(subclass, ["taxonomies", "taxonomy_files"]) or NotImplemented

    def taxonomies(self) -> Generator[Tuple[str, taxonomy.TaxonomyEntry], None, None]:
        for path in self.taxonomy_files:
            taxonomy_file_name = Path(path).name
            with open_file(self.path / path) as file:
                for taxonomy_entry in taxonomy.read(file):
                    yield taxonomy_file_name, taxonomy_entry

    @property
    @abc.abstractmethod
    def taxonomy_files(self) -> List[Union[str, Path]]:
        raise NotImplementedError()
