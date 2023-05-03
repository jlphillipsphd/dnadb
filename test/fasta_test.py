import io
import sys
import unittest

sys.path.append("./src")

from dnadb import fasta

FASTA_SAMPLE = """\
>12345
TACGTAGGGTGCAAGCGTTAAACGGAATTACTGGGCGTAAAGCGTGCGAAGGCGGTTTTATAAGTCTGTAGTGAAAGCACCGGGCTCAACCTGGGAAATGCGAACGAGACTGCAAGGCTTAAATATGGCAGAGGTGGGTAGAATTACACGT
>12346
AACGTAGGTACCGAGCGTTATCCGGATTTACTGGGCGTAAAGCGTGTTCAGGCGGCCTGGCAAGTCGGGCATGAAATCTCTCGGCTCAACCGAGAGAGGCTGTCCGATACTGCTGGGCTTGAGGACGGTAGAGGGTGGTGGAATTCCGCGT
"""

class TestFastaEntry(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))

    def test_length(self):
        self.assertEqual(len(self.fasta_entries), 2)

    def test_entry_sequence_id(self):
        self.assertEqual(str(self.fasta_entries[0].identifier), self.fasta_lines[0][1:].rstrip())
        self.assertEqual(str(self.fasta_entries[1].identifier), self.fasta_lines[2][1:].rstrip())

    def test_entry_sequence(self):
        self.assertEqual(self.fasta_entries[0].sequence, self.fasta_lines[1])
        self.assertEqual(self.fasta_entries[1].sequence, self.fasta_lines[3])

    def test_write(self):
        fasta_file = io.StringIO()
        fasta.write(fasta_file, self.fasta_entries)
        self.assertEqual(fasta_file.getvalue().rstrip(), FASTA_SAMPLE.rstrip())


class TestFastaDb(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create DB
        factory = fasta.FastaDbFactory("/tmp/test.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        # Open DB for testing
        self.db = fasta.FastaDb("/tmp/test.db")

    def test_length(self):
        self.assertEqual(len(self.db), 2)

    def test_has_sequence_index(self):
        self.assertIn(0, self.db)

    def test_has_sequence_id(self):
        self.assertIn(self.fasta_entries[0].identifier, self.db)

    def test_has_sequence_entry(self):
        self.assertIn(self.fasta_entries[0], self.db)

    def test_iter_sequence_entries(self):
        for entry in self.db:
            self.assertIn(entry, self.fasta_entries)

    def test_get_sequence(self):
        self.assertEqual(self.db[0], self.fasta_entries[0])
        self.assertEqual(self.db[1], self.fasta_entries[1])

    def test_get_sequence_by_id(self):
        self.assertEqual(self.db["12345"], self.fasta_entries[0])
        self.assertEqual(self.db["12346"], self.fasta_entries[1])


if __name__ == "__main__":
    unittest.main()
