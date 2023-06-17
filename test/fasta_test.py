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
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create DB
        factory = fasta.FastaDbFactory("/tmp/test.fasta.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        # Open DB for testing
        self.db = fasta.FastaDb(factory.path)

    def tearDown(self):
        """
        Remove the taxonomy database.
        """
        self.db.close()
        for f in self.db.path.iterdir():
            f.unlink()
        self.db.path.rmdir()

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


class TestFastaIndexDb(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create FASTA DB
        factory = fasta.FastaDbFactory("/tmp/test.fasta.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        # Create FASTA Index DB
        factory = fasta.FastaIndexDbFactory("/tmp/test.index.db")
        for i, entry in enumerate(self.fasta_entries):
            factory.write_entry(entry.identifier, f"otu_{i}")
        factory.close()
        # Open DB for testing
        self.index_db = fasta.FastaIndexDb(factory.path)

    def tearDown(self) -> None:
        for db in (self.index_db,):
            db.close()
            for f in db.path.iterdir():
                f.unlink()
            db.path.rmdir()

    def test_length(self):
        self.assertEqual(len(self.index_db), 2)

    def test_contains_fasta_id(self):
        self.assertTrue(self.index_db.contains_fasta_id("12345"))
        self.assertTrue(self.index_db.contains_fasta_id("12346"))

    def test_contains_index(self):
        self.assertTrue(self.index_db.contains_index(0))
        self.assertTrue(self.index_db.contains_index(1))

    def test_contains_key(self):
        self.assertTrue(self.index_db.contains_key("otu_0"))
        self.assertTrue(self.index_db.contains_key("otu_1"))

    def test_fasta_id_to_index(self):
        self.assertEqual(self.index_db.fasta_id_to_index("12345"), 0)
        self.assertEqual(self.index_db.fasta_id_to_index("12346"), 1)

    def test_fasta_id_to_key(self):
        self.assertEqual(self.index_db.fasta_id_to_key("12345"), "otu_0")
        self.assertEqual(self.index_db.fasta_id_to_key("12346"), "otu_1")

    def test_index_to_fasta_id(self):
        self.assertEqual(self.index_db.index_to_fasta_id(0), "12345")
        self.assertEqual(self.index_db.index_to_fasta_id(1), "12346")

    def test_index_to_key(self):
        self.assertEqual(self.index_db.index_to_key(0), "otu_0")
        self.assertEqual(self.index_db.index_to_key(1), "otu_1")

    def test_key_to_fasta_id(self):
        self.assertEqual(self.index_db.key_to_fasta_id("otu_0"), "12345")
        self.assertEqual(self.index_db.key_to_fasta_id("otu_1"), "12346")

    def test_key_to_index(self):
        self.assertEqual(self.index_db.key_to_index("otu_0"), 0)
        self.assertEqual(self.index_db.key_to_index("otu_1"), 1)


if __name__ == "__main__":
    unittest.main()
