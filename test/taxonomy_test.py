import io
import sys
import unittest

sys.path.append("./src")

from dnadb import taxonomy

# Sampled from SILVA 138.1
#
TAXONOMY_SAMPLE = """\
AY855839.1.1390\tk__Bacteria; p__; c__; o__; f__; g__; s__
FW343016.1.1511\tk__Bacteria; p__Firmicutes; c__; o__; f__; g__; s__
AY835431.189876.191345\tk__Bacteria; p__Cyanobacteria; c__; o__; f__; g__; s__
FW369114.1.1462\tk__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__; f__; g__; s__
FW369795.1.1413\tk__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Acetobacterales; f__; g__; s__
AY846383.1.1790\tk__Eukaryota; p__Eukaryota; c__Chlorophyceae; o__Sphaeropleales; f__Sphaeropleales; g__Monoraphidium; s__
AB001440.1.1538\tk__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__Pseudomonas; s__test_species
"""

class TestTaxonomySplits(unittest.TestCase):
    def setUp(self):
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        self.taxonomy_entries = list(taxonomy.read(taxonomy_file))

    def test_split_max_depth(self):
        self.assertEqual(
            taxonomy.split_taxonomy(self.taxonomy_entries[0].label, 7),
            ("Bacteria", "", "", "", "", "", ""))
        self.assertEqual(
            taxonomy.split_taxonomy(self.taxonomy_entries[6].label, 7),
            ("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Pseudomonadales", "Pseudomonadaceae", "Pseudomonas", "test_species"))

    def test_split_truncated_depth(self):
        self.assertEqual(
            taxonomy.split_taxonomy(self.taxonomy_entries[0].label, 4),
            ("Bacteria", "", "", ""))
        self.assertEqual(
            taxonomy.split_taxonomy(self.taxonomy_entries[6].label, 4),
            ("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Pseudomonadales"))


class TestTaxonomyReading(unittest.TestCase):
    def setUp(self):
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        self.taxonomy_entries = list(taxonomy.read(taxonomy_file))
        with open("/tmp/test.tsv", 'w') as f:
            f.write(TAXONOMY_SAMPLE)

    def test_read(self):
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        taxonomy_entries = list(taxonomy.read(taxonomy_file))
        self.assertEqual(taxonomy_entries[0].identifier, "AY855839.1.1390")
        self.assertEqual(taxonomy_entries[0].label, "k__Bacteria; p__; c__; o__; f__; g__; s__")
        self.assertEqual(taxonomy_entries[1].identifier, "FW343016.1.1511")
        self.assertEqual(taxonomy_entries[1].label, "k__Bacteria; p__Firmicutes; c__; o__; f__; g__; s__")

    def test_entries_from_path(self):
        taxonomy_entries = list(taxonomy.entries("/tmp/test.tsv"))
        for a, b in zip(taxonomy_entries, self.taxonomy_entries):
            self.assertEqual(a, b)

    def test_entries_from_io(self):
        taxonomy_entries = taxonomy.entries(io.StringIO(TAXONOMY_SAMPLE))
        for a, b in zip(taxonomy_entries, self.taxonomy_entries):
            self.assertEqual(a, b)

    def test_entries_from_entries(self):
        taxonomy_entries = taxonomy.entries(self.taxonomy_entries)
        for a, b in zip(taxonomy_entries, self.taxonomy_entries):
            self.assertEqual(a, b)


class TestTaxonomyEntry(unittest.TestCase):
    def setUp(self):
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        self.taxonomy_lines = TAXONOMY_SAMPLE.split('\n')
        self.taxonomy_entries = list(taxonomy.entries(taxonomy_file))

    def test_length(self):
        self.assertEqual(len(self.taxonomy_entries), 7)

    def test_entry_sequence_id(self):
        for line, entry in zip(self.taxonomy_lines, self.taxonomy_entries):
            self.assertEqual(entry.identifier, line.split('\t')[0])

    def test_entry_taxonomy(self):
        for line, entry in zip(self.taxonomy_lines, self.taxonomy_entries):
            self.assertEqual(entry.label, line.split('\t')[1])

    def test_entry_serialization(self):
        deserialized = taxonomy.TaxonomyEntry.deserialize(self.taxonomy_entries[0].serialize())
        self.assertEqual(deserialized, self.taxonomy_entries[0])

    def test_write(self):
        taxonomy_file = io.StringIO()
        taxonomy.write(taxonomy_file, self.taxonomy_entries)
        self.assertEqual(taxonomy_file.getvalue().rstrip(), TAXONOMY_SAMPLE.rstrip())
