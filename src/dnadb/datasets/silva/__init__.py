import enum
from functools import cached_property
from pathlib import Path
from typing import List, Optional, Union


from . import taxonomy_map as taxmap
from ..dataset import VersionedDataset, RemoteDatasetMixin, InterfacesWithFasta, \
    InterfacesWithTaxonomy
from ...dna import to_dna
from ...utils import compress, open_file
from ... import fasta

class Subunit(enum.Enum):
    LSU = enum.auto() # Large subunit (23S/28S ribosomal RNA)
    SSU = enum.auto() # Small subunit (16S/18S ribosomal RNA)

class Group(enum.Enum):
    Parc = enum.auto()
    Ref = enum.auto()
    RefNr99 = enum.auto()

class Preprocessing(enum.Enum):
    Raw = enum.auto()
    Truncated = enum.auto()
    FullAligned = enum.auto()

class Silva(RemoteDatasetMixin, VersionedDataset, InterfacesWithFasta, InterfacesWithTaxonomy):

    NAME = "Silva"
    DEFAULT_VERSION = "138.1"

    # https://www.arb-silva.de/download/arb-files/
    BASE_URL = "https://ftp.arb-silva.de"

    def __init__(
        self,
        path: Optional[Union[str, Path]] = None,
        version: str = DEFAULT_VERSION,
        subunit: Subunit = Subunit.SSU,
        group: Group = Group.Ref,
        preprocessing: Preprocessing = Preprocessing.Raw,
        force_download: bool = False
    ):
        self.subunit = subunit
        self.group = group
        self.preprocessing = preprocessing
        super().__init__(path=path, version=version, force_download=force_download)
        if not (self.path / self.fasta_file).exists() or force_download:
            self.__build_dna_fasta()
        if not (self.path / self.taxonomy_file).exists() or force_download:
            self.__build_taxonomy_map()

    def __build_dna_fasta(self):
        print("Converting FASTA to DNA...")
        fasta_path = str(self.path / self.fasta_file)[:-3] # remove .gz extension

        def map_to_dna(entry: fasta.FastaEntry) -> fasta.FastaEntry:
            entry.sequence = to_dna(entry.sequence)
            return entry
        entry_iterator = fasta.entries(self.path / self.__rna_fasta_file)
        entry_iterator = map(map_to_dna, entry_iterator)
        with open(fasta_path, 'w') as file:
            fasta.write(file, entry_iterator)

        print("Compressing FASTA...")
        compress(fasta_path)

    def __build_taxonomy_map(self):
        print("Building taxonomy map...")

        taxonomy_file_path = str(self.path / self.taxonomy_file)[:-3] # remove .gz extension
        Path(taxonomy_file_path).parent.mkdir(exist_ok=True)

        input_taxonomy = open_file(self.path / self.__taxonomy_file, 'r')
        input_taxonomy_tree = open_file(self.path / self.__taxonomy_tree_file, 'r')
        input_taxonomy_map = open_file(self.path / self.__taxonomy_map_file, 'r')
        ouput_taxonomy = open(taxonomy_file_path, 'w')

        tax_dict = taxmap.make_taxid_dict(input_taxonomy)
        sts_dict = taxmap.build_base_silva_taxonomy(input_taxonomy_tree, tax_dict)
        prop_dict = taxmap.propagate_upper_taxonomy(sts_dict, taxmap.rank_prefixes)
        taxmap_dict = taxmap.make_acc_to_species_tid_dict(input_taxonomy_map)

        taxmap.write_tax_strings(taxmap_dict, prop_dict, ouput_taxonomy, sp_label=False)
        ouput_taxonomy.close()

        print("Compressing taxonomy map...")
        compress(taxonomy_file_path)

    def __fasta_file_name(self) -> str:
        suffix = ""
        if self.preprocessing == Preprocessing.Truncated:
            suffix = "_trunc"
        elif self.preprocessing == Preprocessing.FullAligned:
            suffix = "_full_align"
        filename = "SILVA_{}_{}{}_tax_silva{}.fasta.gz".format(
            self.version,
            self.subunit.name,
            "Ref_NR99" if self.group == Group.RefNr99 else self.group.name,
            suffix)
        return filename

    def __taxonomy_file_path(self, name: str) -> Path:
        group = ""
        if self.group in (Group.Parc, Group.Ref):
            group = self.group.name.lower()
        elif self.group == Group.RefNr99:
            group = "ref_nr"
        prefix = Path(f"release_{self.version}/Exports/taxonomy")
        filename = name.format(subunit=self.subunit.name.lower(), group=group, version=self.version)
        return prefix / filename

    @property
    def name(self) -> str:
        return Silva.NAME

    @property
    def url(self) -> str:
        return self.BASE_URL

    @property
    def remote_files(self) -> List[Path]:
        return [
            self.__rna_fasta_file,
            self.__taxonomy_file,
            self.__taxonomy_tree_file,
            self.__taxonomy_map_file
        ]

    @property
    def __rna_fasta_file(self) -> Path:
        return Path(f"release_{self.version}") / "Exports" / self.__fasta_file_name()

    @property
    def __taxonomy_file(self) -> Path:
        return Path(self.__taxonomy_file_path("tax_slv_{subunit}_{version}.txt.gz"))

    @property
    def __taxonomy_tree_file(self) -> Path:
        return Path(self.__taxonomy_file_path("tax_slv_{subunit}_{version}.tre.gz"))

    @property
    def __taxonomy_map_file(self) -> Path:
        return Path(self.__taxonomy_file_path("taxmap_slv_{subunit}_{group}_{version}.txt.gz"))

    @cached_property
    def fasta_file(self) -> Path:
        return Path(f"release_{self.version}") / "cache" / self.__fasta_file_name()

    @cached_property
    def taxonomy_file(self) -> Path:
        return Path(f"release_{self.version}") / "cache" / "taxonomy.txt.gz"
