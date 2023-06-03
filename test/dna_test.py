import numpy as np
import unittest

from dnadb import dna

class TestEncodeDecode(unittest.TestCase):
    def test_encode(self):
        np.testing.assert_array_equal(dna.encode_sequence('ACGT'), [0, 1, 2, 3])

    def test_decode(self):
        np.testing.assert_array_equal(dna.decode_sequence(np.array([0, 1, 2, 3])), 'ACGT')

class TestEncodeDecodeKmers(unittest.TestCase):
    def test_encode_kmers(self):
        sequences = np.array([
            dna.encode_sequence("AAACT"),
            dna.encode_sequence("ACTGC"),
            dna.encode_sequence("TTTAA")
        ])
        result = [[0, 1, 7],
                  [7, 30, 57],
                  [63, 60, 48]]
        np.testing.assert_array_equal(dna.encode_kmers(sequences, 3), result)

    def test_decode_kmers(self):
        sequences = np.array([
            dna.encode_sequence("AAACT"),
            dna.encode_sequence("ACTGC"),
            dna.encode_sequence("TTTAA")
        ])
        kmers = dna.encode_kmers(sequences, 3)
        np.testing.assert_array_equal(dna.decode_kmers(kmers, 3), sequences)


if __name__ == "__main__":
    unittest.main()
