import io
import numpy as np
import sys
from typing import cast
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

    def test_loaded_id_map(self):
        self.assertIsNone(self.db._id_map)

    def test_loaded_sequences(self):
        self.assertIsNone(self.db._sequences)

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


class TestFastaDbLoadedIntoMemory(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create DB
        factory = fasta.FastaDbFactory("/tmp/test.fasta.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        # Open DB for testing
        self.db = fasta.FastaDb(factory.path, load_id_map_into_memory=True, load_sequences_into_memory=True)

    def tearDown(self):
        """
        Remove the taxonomy database.
        """
        self.db.close()
        for f in self.db.path.iterdir():
            f.unlink()
        self.db.path.rmdir()

    def test_loaded_id_map(self):
        self.assertIsNotNone(self.db._id_map)
        assert self.db._id_map is not None # needed to remove warning below for some reason...
        self.assertDictEqual(self.db._id_map, {self.fasta_entries[0].identifier: 0, self.fasta_entries[1].identifier: 1})

    def test_loaded_sequences(self):
        self.assertIsNotNone(self.db._sequences)
        self.assertTrue(np.all(self.db._sequences == self.fasta_entries))

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


class TestFastaMappingEntryFactory(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create DB
        factory = fasta.FastaDbFactory("/tmp/test.fasta.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        # Open DB for testing
        self.db = fasta.FastaDb(factory.path, load_id_map_into_memory=True, load_sequences_into_memory=True)
        self.factory = fasta.FastaMappingDbFactory("/tmp/test.fasta.mapping.db", self.db)
        self.mapping_entry_factory = fasta.FastaMappingEntryFactory("Test", self.factory)

    def tearDown(self) -> None:
        self.db.close()
        for f in self.db.path.iterdir():
            f.unlink()
        self.db.path.rmdir()
        self.factory.close()
        for f in self.factory.path.iterdir():
            f.unlink()
        self.factory.path.rmdir()

    def test_write_entries_are_sorted(self):
        self.mapping_entry_factory.write_entries([self.fasta_entries[1], self.fasta_entries[0]])
        self.assertEqual(self.mapping_entry_factory.sequence_indices[0], 0)
        self.assertEqual(self.mapping_entry_factory.sequence_indices[1], 1)


class TestFastaMappingDb(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create DB
        factory = fasta.FastaDbFactory("/tmp/test.fasta.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        # Open DB for testing
        self.db = fasta.FastaDb(factory.path, load_id_map_into_memory=True, load_sequences_into_memory=True)
        # Create mapping DB
        factory = fasta.FastaMappingDbFactory("/tmp/test.fasta.mapping.db", self.db)
        with factory.create_entry("Test 1") as entry:
            entry.write_entries(self.fasta_entries)
        with factory.create_entry("Test 2") as entry:
            entry.write_entry(self.fasta_entries[1])
        factory.close()
        self.mapping_db = fasta.FastaMappingDb(factory.path, self.db)

    def test_length(self):
        self.assertEqual(len(self.mapping_db), 2)

    def test_mapping_lengths(self):
        self.assertEqual(len(self.mapping_db[0]), 2)
        self.assertEqual(len(self.mapping_db[1]), 1)

    def test_mapping_names(self):
        self.assertEqual(self.mapping_db[0].name, "Test 1")
        self.assertEqual(self.mapping_db[1].name, "Test 2")

    def test_get_sequences(self):
        self.assertEqual(self.mapping_db[0][0], self.fasta_entries[0])
        self.assertEqual(self.mapping_db[0][1], self.fasta_entries[1])

    def test_contains_fasta_entry(self):
        self.assertIn(self.fasta_entries[0], self.mapping_db[0])
        self.assertIn(self.fasta_entries[1], self.mapping_db[0])
        self.assertIn(self.fasta_entries[1], self.mapping_db[1])
        self.assertNotIn(self.fasta_entries[0], self.mapping_db[1])

    def test_contains_fasta_entry_by_id(self):
        self.assertIn(self.fasta_entries[0].identifier, self.mapping_db[0])
        self.assertIn(self.fasta_entries[1].identifier, self.mapping_db[0])
        self.assertIn(self.fasta_entries[1].identifier, self.mapping_db[1])
        self.assertNotIn(self.fasta_entries[0].identifier, self.mapping_db[1])

    def test_contains_fasta_entry_by_index(self):
        self.assertIn(0, self.mapping_db[0])
        self.assertIn(1, self.mapping_db[0])
        self.assertIn(1, self.mapping_db[1])
        self.assertNotIn(0, self.mapping_db[1])

    def test_iter_sequences(self):
        self.assertEqual(list(self.mapping_db[0]), self.fasta_entries)
        self.assertEqual(list(self.mapping_db[1]), [self.fasta_entries[1]])

    def test_get_sequences_by_id(self):
        self.assertEqual(self.mapping_db[0]["12345"], self.fasta_entries[0])
        self.assertEqual(self.mapping_db[0]["12346"], self.fasta_entries[1])
        self.assertEqual(self.mapping_db[1]["12346"], self.fasta_entries[1])

    def test_get_sequences_by_index(self):
        self.assertEqual(self.mapping_db[0][0], self.fasta_entries[0])
        self.assertEqual(self.mapping_db[0][1], self.fasta_entries[1])
        self.assertEqual(self.mapping_db[1][0], self.fasta_entries[1])

if __name__ == "__main__":
    unittest.main()
