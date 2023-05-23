from dataclasses import dataclass
import itertools
import numpy as np
import numpy.typing as npt
import scipy as sp

BASES = "ACGT"
AMBIGUOUS_BASES = "MRWSYKVHDBN"
ALL_BASES = BASES + AMBIGUOUS_BASES

# Map each ambiguous base to a tuple of bases as per the IUAPC standard.
# https://droog.gs.washington.edu/parc/images/iupac.html
# https://iubmb.qmul.ac.uk/misc/naseq.html
IUPAC_MAP = {__b: __c for __b, __c in zip(ALL_BASES, [
    __c for __n in range(1, 5) for __c in itertools.combinations(BASES, __n)])}

# Hashtables ---------------------------------------------------------------------------------------

# A lookup table for converting bases to integers.
# The index is the ASCII value of the base minus the ASCII value of 'A'.
BASE_LOOKUP_TABLE = np.full(np.max(list(map(ord, ALL_BASES))) - ord('A') + 1, 255, dtype=np.uint8)
for __i, __base in enumerate(ALL_BASES):
    BASE_LOOKUP_TABLE[ord(__base) - ord('A')] = __i

# A reverse lookup table for converting integers to bases.
BASE_REVERSE_LOOKUP_TABLE = np.array(list(map(ord, ALL_BASES)), dtype=np.uint8)

# A lookup table for augmenting ambiguous bases to concrete bases.
__iupac_lcm = np.lcm.reduce([len(v) for v in IUPAC_MAP.values()])
IUPAC_AUGMENT_LOOKUP_TABLE = np.full((len(IUPAC_MAP), __iupac_lcm), 255, dtype=np.uint8)
for __base, __bases in IUPAC_MAP.items():
    __repeat = __iupac_lcm // len(__bases)
    __repeated_bases = BASE_LOOKUP_TABLE[np.array(list(map(ord, __bases*__repeat)), np.uint8) - ord('A')]
    IUPAC_AUGMENT_LOOKUP_TABLE[BASE_LOOKUP_TABLE[ord(__base) - ord('A')]] = __repeated_bases

# DNA Sequence Encoding/Decoding -------------------------------------------------------------------

def encode(ascii_bases: npt.NDArray[np.uint8]) ->  npt.NDArray[np.uint8]:
    """
    Encode the given DNA bases in ASCII form into an integer vector representation.
    """
    return BASE_LOOKUP_TABLE[ascii_bases - 65]


def decode(sequences: npt.NDArray[np.uint8]) -> npt.NDArray[np.uint8]:
    """
    Decoed the given DNA bases into ASCII an ascii vector representation.
    """
    return BASE_REVERSE_LOOKUP_TABLE[sequences]


def encode_sequence(sequence: str) -> npt.NDArray[np.uint8]:
    """
    Encode a DNA sequence into an integer vector representation.
    """
    return encode(np.frombuffer(sequence.encode(), np.uint8))


def decode_sequence(sequence: npt.NDArray[np.uint8]) -> str:
    """
    Decode a DNA sequence integer vector representation into a string of bases.
    """
    return decode(sequence).tobytes().decode()


def encode_kmers(
    sequences: npt.NDArray[np.uint8],
    kmer: int,
    ambiguous_bases: bool = False
) -> npt.NDArray[np.int64]:
    """
    Convert DNA sequences into sequences of k-mers.
    """
    slices = [slice(0, s) for s in sequences.shape[:-1]]
    edge_slices = slice((kmer - 1) // 2, (kmer - 1) // -2)
    num_bases = len(BASES + (AMBIGUOUS_BASES if ambiguous_bases else ""))
    powers = np.arange(kmer).reshape((1,)*len(slices) + (-1,))
    kernel = num_bases**powers
    return sp.ndimage.convolve(sequences, kernel)[(*slices, edge_slices)]


def decode_kmers(
    sequences: np.ndarray,
    kmer: int,
    ambiguous_bases: bool = False
) -> npt.NDArray[np.uint8]:
    """
    Decode sequence of k-mers into 1-mer DNA sequences.
    """
    slices = [slice(0, s) for s in sequences.shape[:-1]]
    edge_slice = slice(-1, sequences.shape[-1])
    num_bases = len(BASES + (AMBIGUOUS_BASES if ambiguous_bases else ""))
    powers = np.arange(kmer - 1, -1, -1)
    kernel = num_bases**powers
    edge = (sequences[(*slices, edge_slice)] % kernel[:-1]) // kernel[1:]
    return np.concatenate([sequences // kernel[0], edge], axis=-1).astype(np.uint8)


def augment_ambiguous_bases(
    sequence: str,
    rng: np.random.Generator = np.random.default_rng()
) -> str:
    """
    Replace the ambiguous bases in a DNA sequence at random with a valid concrete base.
    """
    return decode_sequence(replace_ambiguous_encoded_bases(encode_sequence(sequence), rng))


def replace_ambiguous_encoded_bases(
    encoded_sequences: npt.NDArray[np.uint8],
    rng: np.random.Generator = np.random.default_rng()
) -> npt.NDArray[np.uint8]:
    """
    Replace the ambiguous bases in an encoded DNA sequence at random with a valid concrete base.
    """
    augment_indices = rng.integers(0, 12, size=encoded_sequences.shape)
    return IUPAC_AUGMENT_LOOKUP_TABLE[encoded_sequences, augment_indices]


def to_rna(dna_sequence: str) -> str:
    """
    Convert an RNA sequence to DNA.
    """
    return dna_sequence.replace('T', 'U')


def to_dna(rna_sequence: str) -> str:
    """
    Convert a DNA sequence to RNA.
    """
    return rna_sequence.replace('U', 'T')

# Data Classes -------------------------------------------------------------------------------------

@dataclass(frozen=True)
class AbstractSequenceWrapper:
    __slots__ = ("sequence",)

    sequence: str

    def __len__(self):
        return len(self.sequence)
