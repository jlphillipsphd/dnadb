import itertools
import numpy as np
import numpy.typing as npt
import scipy as sp

BASES = "ACGT"
AMBIGUOUS_BASES = "MRWSYKVHDBN" # https://iubmb.qmul.ac.uk/misc/naseq.html
ALL_BASES = BASES + AMBIGUOUS_BASES

# Map each base to an integer
BASE_MAP = {c: i for i, c in enumerate(ALL_BASES)}

# Map each ambiguous base to a tuple of bases
AMBIGUOUS_BASE_MAP = {b: c for b, c in zip(AMBIGUOUS_BASES, [
    c for n in range(2, 5) for c in itertools.combinations(BASES, n)])}

# Map each ambiguous base to a tuple of encoded bases (integerss)
ENC_AMBIGUOUS_BASE_MAP = {BASE_MAP[b]: tuple(BASE_MAP[c] for c in cs)
    for b, cs in AMBIGUOUS_BASE_MAP.items()}

# DNA Sequence Encoding/Decoding -------------------------------------------------------------------

def encode_sequence(sequence: str) -> npt.NDArray[np.uint8]:
    """
    Encode a DNA sequence into an integer vector representation.
    """
    return np.array([BASE_MAP[base] for base in sequence], dtype=np.uint8)


def decode_sequence(sequence: npt.ArrayLike) -> str:
    """
    Decode a DNA sequence integer vector representation into a string of bases.
    """
    return ''.join(ALL_BASES[base] for base in np.asanyarray(sequence))


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


def replace_ambiguous_base(base: np.uint8, rng: np.random.Generator = np.random.default_rng()):
    """
    Replace an encoded base with a random (likely) base in its place.
    """
    return rng.choice(ENC_AMBIGUOUS_BASE_MAP[base])


def replace_ambiguous_bases(
    sequences: npt.NDArray[np.uint8],
    rng: np.random.Generator = np.random.default_rng()
) -> npt.NDArray[np.uint8]:
    """
    Replace the ambiguous bases with a random (likely) base in its place given encoded sequences.
    """
    sequences = sequences.copy()
    indices = np.where(sequences >= 4)
    sequences[indices] = [replace_ambiguous_base(b, rng) for b in sequences[indices]]
    return sequences


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
