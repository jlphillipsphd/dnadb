import numpy as np
import numpy.typing as npt
import tensorflow as tf
from typing import Optional

from ... import dna

BASE_LOOKUP_TABLE = tf.constant(dna.BASE_LOOKUP_TABLE, dtype=tf.int32)
BASE_REVERSE_LOOKUP_TABLE = tf.constant(dna.BASE_REVERSE_LOOKUP_TABLE, dtype=tf.int32)
IUPAC_AUGMENT_LOOKUP_TABLE = tf.constant(dna.IUPAC_AUGMENT_LOOKUP_TABLE, dtype=tf.int32)

# DNA Sequence Encoding/Decoding -------------------------------------------------------------------

def encode(
    sequences: str|bytes|npt.NDArray[np.str_]|tf.Tensor,
    kmer: Optional[int|tf.Tensor] = 1,
    ambiguous_bases: bool = False,
    augment_ambiguous_bases: bool = True
) -> tf.Tensor:
    """
    Encode the given DNA sequences.
    """
    result = _encode_bases(sequences)
    if ambiguous_bases and augment_ambiguous_bases:
        result = augment_ambiguous_bases(result)
    if kmer > 1:
        result = encode_kmers(result, kmer, ambiguous_bases)
    return result


def decode(
    sequences: str|bytes|npt.NDArray[np.str_]|tf.Tensor,
    kmer: Optional[int|tf.Tensor] = 1,
    ambiguous_bases: bool = False
) -> tf.Tensor:
    """
    Decode the given DNA sequences.
    """
    if kmer > 1:
        sequences = decode_kmers(ambiguous_bases)
    return _decode_bases(sequences)

def _encode_bases(sequences: str|bytes|npt.NDArray[np.str_]|tf.Tensor) -> tf.Tensor:
    """
    Encode a DNA sequence into an integer vector representation.
    """
    ascii = tf.cast(tf.io.decode_raw(sequences, tf.uint8), tf.int32)
    return tf.gather(BASE_LOOKUP_TABLE, ascii - 65)


def _decode_bases(sequences: tf.Tensor) -> tf.Tensor:
    """
    Decode a DNA sequence integer vector representation into a string of bases.
    """
    ascii = tf.gather(BASE_REVERSE_LOOKUP_TABLE, sequences)
    return tf.strings.unicode_encode(ascii, output_encoding="UTF-8")


def _encode_kmers(
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


def _decode_kmers(
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


def _augment_ambiguous_bases(sequences: tf.Tensor):
    sequences = tf.expand_dims(sequences, -1)
    return tf.gather_nd(
        IUPAC_AUGMENT_LOOKUP_TABLE,
        tf.concat((
            sequences,
            tf.random.uniform(tf.shape(sequences), minval=0, maxval=len(dna.ALL_BASES), dtype=tf.int32)),
            axis=-1))