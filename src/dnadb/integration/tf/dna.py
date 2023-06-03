import numpy as np
import numpy.typing as npt
import tensorflow as tf

from ... import dna

BASE_LOOKUP_TABLE = tf.constant(dna.BASE_LOOKUP_TABLE, dtype=tf.int32)
BASE_REVERSE_LOOKUP_TABLE = tf.constant(dna.BASE_REVERSE_LOOKUP_TABLE, dtype=tf.int32)
IUPAC_AUGMENT_LOOKUP_TABLE = tf.constant(dna.IUPAC_AUGMENT_LOOKUP_TABLE, dtype=tf.int32)

# DNA Sequence Encoding/Decoding -------------------------------------------------------------------

def encode(sequences: str|npt.NDArray[np.str_]|tf.Tensor) -> tf.Tensor:
    """
    Encode a DNA sequence into an integer vector representation.
    """
    ascii = tf.cast(tf.io.decode_raw(sequences, tf.uint8), tf.int32)
    return tf.gather(BASE_LOOKUP_TABLE, ascii - 65)


def decode(sequences: tf.Tensor) -> str:
    """
    Decode a DNA sequence integer vector representation into a string of bases.
    """
    ascii = tf.gather(BASE_REVERSE_LOOKUP_TABLE, sequences)
    return tf.strings.unicode_encode(ascii, output_encoding="UTF-8")


def encode_kmers(
    sequences: tf.Tensor,
    kmer: int|tf.Tensor,
    ambiguous_bases: bool = False
) -> tf.Tensor:
    """
    Convert DNA sequences into sequences of k-mers.
    """
    original_shape = tf.shape(sequences)
    sequence_length = tf.shape(sequences)[-1]
    sequences = tf.reshape(sequences, (-1, sequence_length, 1))
    num_bases = len(dna.BASES + (dna.AMBIGUOUS_BASES if ambiguous_bases else ""))
    kernel = tf.reshape(tf.pow(num_bases, tf.range(kmer-1, -1, -1, dtype=tf.int32)), (kmer, 1, 1))
    result = tf.nn.convolution(sequences, kernel, padding="VALID")
    return tf.reshape(
        result,
        tf.concat((original_shape[:-1], (sequence_length - kmer + 1,)), axis=0))


def decode_kmers(
    kmer_sequences: tf.Tensor,
    kmer: int|tf.Tensor,
    ambiguous_bases: bool = False
) -> tf.Tensor:
    """
    Decode sequence of k-mers into 1-mer DNA sequences.
    """
    num_bases = len(dna.BASES + (dna.AMBIGUOUS_BASES if ambiguous_bases else ""))
    powers = tf.range(kmer - 1, -1, -1)
    kernel = num_bases**powers
    return tf.concat([
        kmer_sequences // kernel[0],
        tf.repeat(kmer_sequences[:,-1:], kmer - 1, axis=-1) % kernel[:-1] // kernel[1:]
    ], axis=-1)


# def augment_ambiguous_bases(
#     sequence: str,
#     rng: np.random.Generator = np.random.default_rng()
# ) -> str:
#     """
#     Replace the ambiguous bases in a DNA sequence at random with a valid concrete base.
#     """
#     return decode_sequence(replace_ambiguous_encoded_bases(encode_sequence(sequence), rng))


# def replace_ambiguous_encoded_bases(
#     encoded_sequences: npt.NDArray[np.uint8],
#     rng: np.random.Generator = np.random.default_rng()
# ) -> npt.NDArray[np.uint8]:
#     """
#     Replace the ambiguous bases in an encoded DNA sequence at random with a valid concrete base.
#     """
#     augment_indices = rng.integers(0, 12, size=encoded_sequences.shape)
#     return IUPAC_AUGMENT_LOOKUP_TABLE[encoded_sequences, augment_indices]


# def to_rna(dna_sequence: str) -> str:
#     """
#     Convert an RNA sequence to DNA.
#     """
#     return dna_sequence.replace('T', 'U')


# def to_dna(rna_sequence: str) -> str:
#     """
#     Convert a DNA sequence to RNA.
#     """
#     return rna_sequence.replace('U', 'T')
