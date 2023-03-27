import itertools
import numpy as np
import numpy.typing as npt

BASES = "ACGT"
INCOMPLETE_BASES = "MRWSYKVHDBN" # https://iubmb.qmul.ac.uk/misc/naseq.html
ALL_BASES = BASES + INCOMPLETE_BASES

BASE_MAP = {c: i for i, c in enumerate(ALL_BASES)}
INCOMPLETE_BASE_MAP = {b: c for b, c in zip(INCOMPLETE_BASES, [
    c for n in range(2, 5) for c in itertools.combinations(BASES, n)])}
ENC_INCOMPLETE_BASE_MAP = {BASE_MAP[b]: tuple(BASE_MAP[c] for c in cs)
    for b, cs in INCOMPLETE_BASE_MAP.items()}

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
