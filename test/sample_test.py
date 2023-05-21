import io
import numpy as np
import sys
import unittest
from unittest.mock import MagicMock
from typing import List

from .mocking import MockableRng

sys.path.append("./src")

from dnadb import fasta, fastq, sample

FASTA_SAMPLE = """\
>aaa
TACGTAGGGTGCAAGCGTTAAACGGAATTACTGGGCGTAAAGCGTGCGAAGGCGGTTTTATAAGTCTGTAGTGAAAGCACCGGGCTCAACCTGGGAAATGCGAACGAGACTGCAAGGCTTAAATATGGCAGAGGTGGGTAGAATTACACGT
>bbb
AACGTAGGTACCGAGCGTTATCCGGATTTACTGGGCGTAAAGCGTGTTCAGGCGGCCTGGCAAGTCGGGCATGAAATCTCTCGGCTCAACCGAGAGAGGCTGTCCGATACTGCTGGGCTTGAGGACGGTAGAGGGTGGTGGAATTCCGCGT
>ccc
CGAGCGTTATCCCGGCCTGGCAAGTCGGGCATGAAATCTCTCGGCTCAACCGAGAGAGGCTGTCCGATACTGCTGGGCTTGAGGACGGTAGAGGGTGGTGGAATTCCGCGTGGATTTAACGTAGGTACACTGGGCGTAAAGCGTGTTCAGG
>ddd
TCCGGCGTAAAGCGTGTTCAGGCGGCCTGGCAAGTCGGGCATGAAATCTCTCGGCTCAACCGAGAGAGGCTGGGGCTTGAGGACGGTAGAGGGTGGTGGAATTCCGCGTGATACTGCTAACGTAGGTACCGAGCGTTATCCGGATTTACTG
>eee
GGACGGTAGAGGGTGGTGGAATTCCGCGTAACGTAGGTACCGAGCGTTATCCGGATTTACTGGGCGTAAAGCGTGTTCAGGCGGCCTGGCAAGTCGGGCATGAAATCTCTCGGCTCAACCGAGAGAGGCTGTCCGATACTGCTGGGCTTGA
"""

FASTQ_SAMPLE = """\
@MN00371:50:000H2735W:1:11102:4617:1056 1:N:0:TAGGTGGAAT
GACGAAGGATCCAAAAAAAATCCACATACATTAGATTAANACGGCCAGTAGCTGGCCAAAAAACACACACAAANAAAATCTGAGCTAAAACCACAAATTACAAATGATACTGTAGCCCATCAGAATAAATCCCCATGCCGCAAAAATACAT
+
/F=/6AF/FAAA=//F///FFFF/F/FA//F///F/AFF#/AA=/A//AFA/AF//F//////FF/FF/=/F6#6F///AAA///////F/=FFF//F/=//A/F/FF/=A/FFFA=/FFFFFF/F/AAF/AFF/F/AF/F/FFF///=//
@MN00371:50:000H2735W:1:11102:11664:1070 1:N:0:TAGGTGGAAT
CACGTAGGGTGCGAGCGTTTTCCGGAATTACTGGGCGTAAAGCGCGCGCAGGCGGCTTCGCGCGCCCGCCGTGAAAGCCCCCGGCTTAACCGGGGAGAGTCGGTGGGGACGGCGGAGCTTGAGGGCGGGAGAGGTCGGTGGAATTCCCGGT
+
F/FAFFFFFFFFFFFFFFF=FFFFFFFFFFFFFFFFFFFFFFFFAFFFF=FFFFFFFFFFFFA/AFAFFFFFFFFFAAF=FFFFFFFAFFFFFF/FF6FFFFFFF/FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAFFFF/FFFFF
@MN00371:50:000H2735W:1:11102:24608:1073 1:N:0:TAGGTGGAAT
TCCAGACCAAGTCTCTGCTACGTATTAGATACCCGTGTAGTCCAGACCAAGTCTCTGCTACCGTATAGGTGGAATATATCGTATGCCGTCAAGTGGTTGAAAAAAAAAACGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+
AFAFFFFFFFFFAFFFFFFA=FFFFFFFFFFFFFFFFF/FFF/FFF/FF/FF/F/=/FFF/FFFFFFAAFFFF/F/////A//F///=F/////////AFF/AFF///////AFAFFFFFF/FAFFFF/FAFFFFFFAFFFFFFFFFFFFF
@MN00371:50:000H2735W:1:11102:10397:1085 1:N:0:TAGGTGGAAT
TACAGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGGGTGCGTAGGCGGTTATATAAGTTAGATGTGAAAGCCCTAGGCTTAACCTGGGAACTGCGTTTAATACTGTATAGCTAGAGTACTGAAGAGGTTAGTGGAATTTCCAGT
+
FAAFFFFF=/AF/FFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAFFFF=FFFF/FFFF/FFFF=FAAFF/FF/FF/FFF/FFF/A//FFFFFAFFFF/FFFFFFF/F/F/FFFFFAFF/FFFFFFFFFFFFFFFFF/FFFFFFF
@MN00371:50:000H2735W:1:11102:19984:1087 1:N:0:TAGGTGGAAT
TACGGAGGGTGCTAGCGTTGTTCGGAATCACTGGGCGTAAAGGGCGTGTAGGCGGTCTGTTAAGTCTGATGTGAAATCCCCTCGCTCAACGAGGGAACTGCGTCGGATACTGGCAGACTTGAGTACCGGAGAGGAAGGTGGAATTCCCGGT
+
FFA=FFFFAFFFFFFFFFFFFFAFF/FFFAFFFFFFFFFFFFFFFFFFFAFFFFFAFFFFFFFFFFFFFFFFFFAFFFFFFFFFFFFFFFF6FFFFFAFAFFFFFFFFFFFFFFF/FFFFFFFFFFAFF/FFFFFFFFFFFFFFFFAFAFF
"""

def create_samples(fasta_entries: List[fasta.FastaEntry], index_db: fasta.FastaIndexDb):
    sample_a = sample.SampleMappingEntryFactory("Sample A", index_db) \
        .add_entry(fasta_entries[1]) \
        .add_entry(fasta_entries[2], abundance=5) \
        .add_entry(fasta_entries[0], abundance=3) \
        .add_entry(fasta_entries[0], abundance=7) \
        .build()
    sample_b = sample.SampleMappingEntryFactory("Sample B", index_db) \
        .add_entry(fasta_entries[3], abundance=10) \
        .add_entry(fasta_entries[4], abundance=20) \
        .build()
    return sample_a, sample_b


class TestSampleMappingEntry(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create index DB
        factory = fasta.FastaIndexDbFactory("/tmp/test.index.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        self.index_db = fasta.FastaIndexDb("/tmp/test.index.db")
        self.sample_entries = create_samples(self.fasta_entries, self.index_db)

    def tearDown(self) -> None:
        self.index_db.close()
        for f in self.index_db.path.iterdir():
            f.unlink()
        self.index_db.path.rmdir()

    def test_length(self):
        self.assertEqual(len(self.sample_entries[0]), 3, self.sample_entries[0])
        self.assertEqual(len(self.sample_entries[1]), 2)

    def test_abundance(self):
        self.assertTrue(np.array_equal(self.sample_entries[0].abundances, np.array([10, 1, 5])))
        self.assertTrue(np.array_equal(self.sample_entries[1].abundances, np.array([10, 20])))

    def test_total_abundance(self):
        self.assertEqual(self.sample_entries[0].total_abundance, 16)
        self.assertEqual(self.sample_entries[1].total_abundance, 30)

    def test_sorted_indices(self):
        self.assertTrue(np.array_equal(np.sort(self.sample_entries[0].indices), np.array([0, 1, 2])))
        self.assertTrue(np.array_equal(np.sort(self.sample_entries[1].indices), np.array([3, 4])))

    def test_contains_fasta_id(self):
        # Sample A
        self.assertTrue(self.sample_entries[0].contains_fasta_id("aaa"))
        self.assertTrue(self.sample_entries[0].contains_fasta_id("bbb"))
        self.assertTrue(self.sample_entries[0].contains_fasta_id("ccc"))
        self.assertFalse(self.sample_entries[0].contains_fasta_id("ddd"))
        self.assertFalse(self.sample_entries[0].contains_fasta_id("eee"))
        # Sample B
        self.assertFalse(self.sample_entries[1].contains_fasta_id("aaa"))
        self.assertFalse(self.sample_entries[1].contains_fasta_id("bbb"))
        self.assertFalse(self.sample_entries[1].contains_fasta_id("ccc"))
        self.assertTrue(self.sample_entries[1].contains_fasta_id("ddd"))
        self.assertTrue(self.sample_entries[1].contains_fasta_id("eee"))


class TestSampleMappingDb(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create DB
        factory = fasta.FastaDbFactory("/tmp/test.fasta.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        self.db = fasta.FastaDb(factory.path)
        # Create index DB
        factory = fasta.FastaIndexDbFactory("/tmp/test.index.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        self.index_db = fasta.FastaIndexDb(factory.path)
        # Create the sample mapping DB
        self.sample_entries = create_samples(self.fasta_entries, self.index_db)
        factory = sample.SampleMappingDbFactory("/tmp/test.mapping.db")
        factory.write_entries(self.sample_entries)
        factory.close()
        # Open DB for testing
        self.sample_db = sample.SampleMappingDb(factory.path, self.index_db)

    def tearDown(self) -> None:
        for db in (self.db, self.index_db, self.sample_db):
            db.close()
            for f in db.path.iterdir():
                f.unlink()
            db.path.rmdir()

    def test_length(self):
        self.assertEqual(len(self.sample_db), 2)

    def test_getitem(self):
        self.assertEqual(self.sample_db[0], self.sample_entries[0])
        self.assertEqual(self.sample_db[1], self.sample_entries[1])
        self.assertEqual(self.sample_db["Sample A"], self.sample_entries[0])
        self.assertEqual(self.sample_db["Sample B"], self.sample_entries[1])


class TestFastaSample(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create DB
        factory = fasta.FastaDbFactory("/tmp/test.fasta.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        self.db = fasta.FastaDb(factory.path)
        self.sample = sample.load_fasta(self.db, "test")

    def tearDown(self) -> None:
        self.db.close()
        for f in self.db.path.iterdir():
            f.unlink()
        self.db.path.rmdir()

    def test_length(self):
        self.assertEqual(len(self.sample), 5)

    def test_name(self):
        self.assertEqual(self.sample.name, "test")

    def test_getitem(self):
        self.assertEqual(self.sample[0], self.fasta_entries[0])
        self.assertEqual(self.sample[1], self.fasta_entries[1])
        self.assertEqual(self.sample[2], self.fasta_entries[2])
        self.assertEqual(self.sample[3], self.fasta_entries[3])
        self.assertEqual(self.sample[4], self.fasta_entries[4])
        self.assertEqual(self.sample["aaa"], self.fasta_entries[0])
        self.assertEqual(self.sample["bbb"], self.fasta_entries[1])
        self.assertEqual(self.sample["ccc"], self.fasta_entries[2])
        self.assertEqual(self.sample["ddd"], self.fasta_entries[3])
        self.assertEqual(self.sample["eee"], self.fasta_entries[4])

    def test_iter(self):
        self.assertEqual(list(self.sample), self.fasta_entries)

    def test_sample(self):
        # Create the mocked RNG
        n = 5
        indices = [4, 2, 0, 3, 0]
        rng = MockableRng()
        rng.choice = MagicMock(return_value=indices)
        # Sample
        samples = list(self.sample.sample(n, rng=rng))
        rng.choice.assert_called_with(5, size=n, replace=True)
        self.assertEqual(sorted(samples), sorted([self.fasta_entries[i] for i in indices]))


class TestFastqSample(unittest.TestCase):
    def setUp(self):
        fastq_file = io.StringIO(FASTQ_SAMPLE)
        self.fastq_lines = FASTQ_SAMPLE.split('\n')
        self.fastq_entries = list(fastq.read(fastq_file))
        # Create DB
        factory = fastq.FastqDbFactory("/tmp/test.fastq.db")
        factory.write_entries(self.fastq_entries)
        factory.close()
        self.db = fastq.FastqDb(factory.path)
        self.sample = sample.load_fastq(self.db, "test")

    def tearDown(self) -> None:
        self.db.close()
        for f in self.db.path.iterdir():
            f.unlink()
        self.db.path.rmdir()

    def test_length(self):
        self.assertEqual(len(self.sample), 5)

    def test_name(self):
        self.assertEqual(self.sample.name, "test")

    def test_getitem(self):
        self.assertEqual(self.sample[0], self.fastq_entries[0])
        self.assertEqual(self.sample[1], self.fastq_entries[1])
        self.assertEqual(self.sample[2], self.fastq_entries[2])
        self.assertEqual(self.sample[3], self.fastq_entries[3])
        self.assertEqual(self.sample[4], self.fastq_entries[4])

    def test_iter(self):
        self.assertEqual(list(self.sample), self.fastq_entries)

    def test_sample(self):
        # Create the mocked RNG
        n = 5
        indices = [4, 2, 0, 3, 0]
        rng = MockableRng()
        rng.choice = MagicMock(return_value=indices)
        # Sample
        samples = list(self.sample.sample(n, rng=rng))
        rng.choice.assert_called_with(5, size=n, replace=True)
        self.assertEqual(sorted(samples), sorted([self.fastq_entries[i] for i in indices]))


class TestDemultiplexedFastaSample(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create DB
        factory = fasta.FastaDbFactory("/tmp/test.fasta.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        self.db = fasta.FastaDb(factory.path)
        # Create index DB
        factory = fasta.FastaIndexDbFactory("/tmp/test.index.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        self.index_db = fasta.FastaIndexDb(factory.path)
        # Create the sample mapping DB
        self.sample_entries = create_samples(self.fasta_entries, self.index_db)
        factory = sample.SampleMappingDbFactory("/tmp/test.mapping.db")
        factory.write_entries(self.sample_entries)
        factory.close()
        # Open DB for testing
        self.sample_db = sample.SampleMappingDb(factory.path, self.index_db)
        self.samples = sample.load_multiplexed_fasta(self.db, self.sample_db)

    def tearDown(self) -> None:
        for db in (self.db, self.index_db, self.sample_db):
            db.close()
            for f in db.path.iterdir():
                f.unlink()
            db.path.rmdir()

    def test_length(self):
        self.assertEqual(len(self.samples[0]), 1+5+3+7)
        self.assertEqual(len(self.samples[1]), 10+20)

    def test_name(self):
        self.assertEqual(self.samples[0].name, "Sample A")
        self.assertEqual(self.samples[1].name, "Sample B")

    def test_getitem(self):
        # Sample A
        self.assertEqual(self.samples[0][0], self.fasta_entries[0])
        self.assertEqual(self.samples[0][9], self.fasta_entries[0])
        self.assertEqual(self.samples[0][10], self.fasta_entries[1])
        self.assertEqual(self.samples[0][11], self.fasta_entries[2])
        self.assertEqual(self.samples[0][15], self.fasta_entries[2])
        self.assertEqual(self.samples[0]["ccc"], self.fasta_entries[2])
        self.assertRaises(KeyError, lambda: self.samples[0]["ddd"])
        self.assertRaises(KeyError, lambda: self.samples[0]["eee"])
        # Sample B
        self.assertEqual(self.samples[1][0], self.fasta_entries[3])
        self.assertEqual(self.samples[1][9], self.fasta_entries[3])
        self.assertEqual(self.samples[1][10], self.fasta_entries[4])
        self.assertEqual(self.samples[1][29], self.fasta_entries[4])
        self.assertRaises(KeyError, lambda: self.samples[1]["aaa"])
        self.assertRaises(KeyError, lambda: self.samples[1]["bbb"])
        self.assertRaises(KeyError, lambda: self.samples[1]["ccc"])

    def test_sample(self):
        # Create the mocked RNG
        n = 5
        indices = [0, 1, 1, 0, 0]
        rng = MockableRng()
        rng.choice = MagicMock(return_value=indices)
        # Sample
        samples = list(self.samples[0].sample(n, rng=rng))
        self.assertEqual(rng.choice.call_count, 1)
        self.assertEqual(rng.choice.call_args.args[0], 3) # n
        self.assertEqual(rng.choice.call_args.kwargs["size"], n) # size
        self.assertEqual(rng.choice.call_args.kwargs["replace"], True) # replace
        np.testing.assert_array_equal(rng.choice.call_args.kwargs["p"], self.samples[0].abundances / self.samples[0].total_abundance) # p
        self.assertEqual(sorted(samples), sorted([self.fasta_entries[i] for i in indices]))


class TestLoadFastaSample(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create DB
        factory = fasta.FastaDbFactory("/tmp/test.fasta.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        self.db = fasta.FastaDb(factory.path)
        # Create index DB
        factory = fasta.FastaIndexDbFactory("/tmp/test.index.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        self.index_db = fasta.FastaIndexDb(factory.path)
        # Create the sample mapping DB
        self.sample_entries = create_samples(self.fasta_entries, self.index_db)
        factory = sample.SampleMappingDbFactory("/tmp/test.mapping.db")
        factory.write_entries(self.sample_entries)
        factory.close()
        # Open DB for testing
        self.sample_db = sample.SampleMappingDb(factory.path, self.index_db)
        self.samples = sample.load_multiplexed_fasta(self.db, self.sample_db)

    def tearDown(self) -> None:
        for db in (self.db, self.index_db, self.sample_db):
            db.close()
            for f in db.path.iterdir():
                f.unlink()
            db.path.rmdir()

    def test_load_from_path(self):
        fasta_samples = sample.load_multiplexed_fasta(self.db.path, self.sample_db.path)
        self.assertEqual(len(fasta_samples), 2)
        self.assertEqual(fasta_samples[0].name, "Sample A")
        self.assertEqual(fasta_samples[1].name, "Sample B")

    def test_load_from_path_with_index(self):
        fasta_samples = sample.load_multiplexed_fasta(self.db.path, self.sample_db.path, self.index_db.path)
        self.assertEqual(len(fasta_samples), 2)
        self.assertEqual(fasta_samples[0].name, "Sample A")
        self.assertEqual(fasta_samples[1].name, "Sample B")

    def test_load_from_db(self):
        fasta_samples = sample.load_multiplexed_fasta(self.db, self.sample_db)
        self.assertEqual(len(fasta_samples), 2)
        self.assertEqual(fasta_samples[0].name, "Sample A")
        self.assertEqual(fasta_samples[1].name, "Sample B")

    def test_load_from_db_with_index(self):
        fasta_samples = sample.load_multiplexed_fasta(self.db, self.sample_db, self.index_db)
        self.assertEqual(len(fasta_samples), 2)
        self.assertEqual(fasta_samples[0].name, "Sample A")
        self.assertEqual(fasta_samples[1].name, "Sample B")


class TestLoadDemultiplexedFastaSample(unittest.TestCase):
    def setUp(self):
        fasta_file = io.StringIO(FASTA_SAMPLE)
        self.fasta_lines = FASTA_SAMPLE.split('\n')
        self.fasta_entries = list(fasta.read(fasta_file))
        # Create DB
        factory = fasta.FastaDbFactory("/tmp/test.fasta.db")
        factory.write_entries(self.fasta_entries)
        factory.close()
        self.db = fasta.FastaDb(factory.path)

    def tearDown(self) -> None:
        self.db.close()
        for f in self.db.path.iterdir():
            f.unlink()
        self.db.path.rmdir()

    def test_load_from_path(self):
        fasta_sample = sample.load_fasta(self.db.path, "test")
        self.assertEqual(fasta_sample.name, "test")

    def test_load_from_db(self):
        fasta_sample = sample.load_fasta(self.db, "test")
        self.assertEqual(fasta_sample.name, "test")

if __name__ == "__main__":
    unittest.main()
