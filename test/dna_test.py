import numpy as np
import unittest

from dnadb import dna

class TestEncode(unittest.TestCase):
    def test_encode(self):
        np.testing.assert_array_equal(dna.encode_sequence('ACGT'), [0, 1, 2, 3])

    def test_decode(self):
        np.testing.assert_array_equal(dna.decode_sequence(np.array([0, 1, 2, 3])), 'ACGT')
