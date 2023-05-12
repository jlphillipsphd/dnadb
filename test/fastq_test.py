import io
import sys
import unittest

sys.path.append("./src")

from dnadb import fastq

FASTQ_SAMPLE = """\
@MN00371:50:000H2735W:1:11102:23071:1116 1:N:0:CACTTGTGTC
TACGTAGGGTGCAAGCGTTAAACGGAATTACTGGGCGTAAAGCGTGCGAAGGCGGTTTTATAAGTCTGTAGTGAAAGCACCGGGCTCAACCTGGGAAATGCGAACGAGACTGCAAGGCTTAAATATGGCAGAGGTGGGTAGAATTACACGT
+
AFA/FFF=F/FFFFFFFFAAF6FFFFFFFF/FFFF/AFAF/A/F/=/FAFFFF/F/F/A/AFFFFFAFF/FFFAF//F/AF/AFFFFFFFFFA//FF/F/F//A//F/FFFA/FF////FAFF/////F//FFF//F/AF//F/A/F///A
@MN00371:50:000H2735W:1:11102:16811:1134 1:N:0:CACTTGTGTG
AACGTAGGTACCGAGCGTTATCCGGATTTACTGGGCGTAAAGCGTGTTCAGGCGGCCTGGCAAGTCGGGCATGAAATCTCTCGGCTCAACCGAGAGAGGCTGTCCGATACTGCTGGGCTTGAGGACGGTAGAGGGTGGTGGAATTCCGCGT
+
FFAAFF=FFFFAF/FFFFFFFF=FFFFFFFFFAF//FFAFFFFFFAFFFFFFFFFFFF/FFFA/AFFAFFFAAFAFAFFF=F/F//FFFFFFFAFF/FFFFFAFFF=FFFAFFFFFAA/A/FAFFF/AF/FAFFFFFFFA/AFFFFFFFF/
"""

class TestPhredEncoding(unittest.TestCase):
    def test_phred_33_encode(self):
        self.assertEqual(fastq.phred_encode([1.0]), "!")
        self.assertEqual(fastq.phred_encode([10**(40 / -10)]), "I")

    def test_phred_33_decode(self):
        self.assertEqual(fastq.phred_decode("!", 33), 1.0)
        self.assertEqual(fastq.phred_decode("I", 33), 10**(40 / -10))


class TestFastqHeader(unittest.TestCase):
    def setUp(self):
        self.header_str: str = "@MN00371:50:000H2735W:1:11102:23071:1116 1:N:0:CACTTGTGTC\n"
        self.header: fastq.FastqHeader = fastq.FastqHeader.from_str(self.header_str)

    def test_instrument(self):
        self.assertEqual(self.header.instrument, "MN00371")

    def test_run_number(self):
        self.assertEqual(self.header.run_number, 50)

    def test_flowcell_id(self):
        self.assertEqual(self.header.flowcell_id, "000H2735W")

    def test_lane(self):
        self.assertEqual(self.header.lane, 1)

    def test_tile(self):
        self.assertEqual(self.header.tile, 11102)

    def test_pos(self):
        self.assertEqual(self.header.pos, (23071, 1116))

    def test_read_type(self):
        self.assertEqual(self.header.read_type, 1)

    def test_is_filtered(self):
        self.assertIs(self.header.is_filtered, False)

    def test_control_number(self):
        self.assertEqual(self.header.control_number, 0)

    def test_sequence_index(self):
        self.assertEqual(self.header.sequence_index, "CACTTGTGTC")

    def test_encode_to_str(self):
        self.assertEqual(str(self.header), self.header_str.rstrip())


class TestFastqEntry(unittest.TestCase):
    def setUp(self):
        fastq_file = io.StringIO(FASTQ_SAMPLE)
        self.fastq_lines = FASTQ_SAMPLE.split('\n')
        self.fastq_entries = list(fastq.read(fastq_file))

    def test_length(self):
        self.assertEqual(len(self.fastq_entries), 2)

    def test_entry_header(self):
        self.assertEqual(str(self.fastq_entries[0].header), self.fastq_lines[0])
        self.assertEqual(str(self.fastq_entries[1].header), self.fastq_lines[4])

    def test_entry_sequence(self):
        self.assertEqual(self.fastq_entries[0].sequence, self.fastq_lines[1])
        self.assertEqual(self.fastq_entries[1].sequence, self.fastq_lines[5])

    def test_entry_qality_scores(self):
        self.assertEqual(self.fastq_entries[0].quality_scores, self.fastq_lines[3])
        self.assertEqual(self.fastq_entries[1].quality_scores, self.fastq_lines[7])

    def test_write(self):
        fastq_file = io.StringIO()
        fastq.write(fastq_file, self.fastq_entries)
        self.assertEqual(fastq_file.getvalue().rstrip(), FASTQ_SAMPLE.rstrip())


class TestFastqDb(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        fastq_file = io.StringIO(FASTQ_SAMPLE)
        self.fastq_lines = FASTQ_SAMPLE.split('\n')
        self.fastq_entries = list(fastq.read(fastq_file))
        # Create DB
        factory = fastq.FastqDbFactory("/tmp/test.db")
        factory.write_entries(self.fastq_entries)
        factory.close()
        # Open DB for testing
        self.db = fastq.FastqDb("/tmp/test.db")

    def test_length(self):
        self.assertEqual(len(self.db), 2)

    def test_has_sequence_index(self):
        self.assertIn(0, self.db)

    def test_iter_sequence_entries(self):
        for entry in self.db:
            self.assertIn(entry, self.fastq_entries)

    def test_get_sequence(self):
        self.assertEqual(self.db[0], self.fastq_entries[0])
        self.assertEqual(self.db[1], self.fastq_entries[1])


if __name__ == "__main__":
    unittest.main()
