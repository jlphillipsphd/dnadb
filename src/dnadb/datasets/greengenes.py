from functools import cached_property
from pathlib import Path
from typing import List, Optional, Union

from .dataset import VersionedDataset, RemoteDatasetMixin, InterfacesWithFasta, \
    InterfacesWithTaxonomy

class Greengenes(RemoteDatasetMixin, VersionedDataset, InterfacesWithFasta, InterfacesWithTaxonomy):

    NAME = "Greengenes"
    DEFAULT_VERSION = "13.5" # 13.8 is not on their website...

    # https://greengenes.secondgenome.com/
    BASE_URL = "https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database"

    def __init__(
        self,
        path: Optional[Union[str, Path]] = None,
        version: str = DEFAULT_VERSION,
        force_download: bool = False
    ):
        super().__init__(path=path, version=version, force_download=force_download)

    @property
    def name(self) -> str:
        return self.NAME

    @property
    def remote_files(self) -> List[Path]:
        return [self.fasta_file, self.taxonomy_file]

    @cached_property
    def fasta_file(self) -> Path:
        version = self.version.replace('.', '_')
        return Path(f"gg_{version}") / f"gg_{version}.fasta.gz"

    @cached_property
    def taxonomy_file(self) -> Path:
        version = self.version.replace('.', '_')
        return Path(f"gg_{version}") / f"gg_{version}_taxonomy.txt.gz"

    @cached_property
    def url(self) -> str:
        return self.BASE_URL
