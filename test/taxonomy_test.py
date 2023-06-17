import io
import os
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
FW369795.1.xxxx\tk__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Acetobacterales; f__; g__; s__
"""

class TestTaxonomySplits(unittest.TestCase):
    def setUp(self):
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        self.taxonomy_entries = list(taxonomy.read(taxonomy_file))

    def test_split(self):
        """
        Split the taxonomy into a tuple of taxons.
        """
        self.assertEqual(
            taxonomy.split_taxonomy(self.taxonomy_entries[0].label),
            ("Bacteria",))
        self.assertEqual(
            taxonomy.split_taxonomy(self.taxonomy_entries[0].label, keep_empty=True),
            ("Bacteria", "", "", "", "", "", ""))
        self.assertEqual(
            taxonomy.split_taxonomy(self.taxonomy_entries[6].label),
            ("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Pseudomonadales",
             "Pseudomonadaceae", "Pseudomonas", "test_species"))

    def test_join(self):
        """
        Join the taxons to the maximum depth.
        """
        self.assertEqual(
            taxonomy.join_taxonomy(taxonomy.split_taxonomy(self.taxonomy_entries[0].label, keep_empty=True)),
            self.taxonomy_entries[0].label)
        self.assertEqual(
            taxonomy.join_taxonomy(taxonomy.split_taxonomy(self.taxonomy_entries[6].label, keep_empty=True)),
            self.taxonomy_entries[6].label)
        self.assertEqual(
            taxonomy.join_taxonomy(["Bacteria"]), "k__Bacteria"
        )
        self.assertEqual(
            taxonomy.join_taxonomy(["Bacteria", ""]), "k__Bacteria; p__"
        )
        self.assertEqual(
            taxonomy.join_taxonomy(["Bacteria"], depth=2), "k__Bacteria; p__"
        )

class TestTaxonomyReading(unittest.TestCase):
    def setUp(self):
        """
        Write the taxonomy sample to a file for the tests to read from.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        self.taxonomy_entries = list(taxonomy.read(taxonomy_file))
        with open("/tmp/test.tsv", 'w') as f:
            f.write(TAXONOMY_SAMPLE)

    def tearDown(self):
        """
        Remove the taxonomy database.
        """
        os.remove("/tmp/test.tsv")

    def test_read(self):
        """
        Read taxonomy entries from a file-like object.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        taxonomy_entries = list(taxonomy.read(taxonomy_file))
        self.assertEqual(taxonomy_entries[0].identifier, "AY855839.1.1390")
        self.assertEqual(taxonomy_entries[0].label, "k__Bacteria; p__; c__; o__; f__; g__; s__")
        self.assertEqual(taxonomy_entries[1].identifier, "FW343016.1.1511")
        self.assertEqual(taxonomy_entries[1].label, "k__Bacteria; p__Firmicutes; c__; o__; f__; g__; s__")

    def test_entries_from_path(self):
        """
        Read taxonomy entries from a file path.
        """
        taxonomy_entries = list(taxonomy.entries("/tmp/test.tsv"))
        for a, b in zip(taxonomy_entries, self.taxonomy_entries):
            self.assertEqual(a, b)

    def test_entries_from_io(self):
        """
        Read taxonomy entries from a file-like object.
        """
        taxonomy_entries = taxonomy.entries(io.StringIO(TAXONOMY_SAMPLE))
        for a, b in zip(taxonomy_entries, self.taxonomy_entries):
            self.assertEqual(a, b)

    def test_entries_from_entries(self):
        """
        Read taxonomy entries from a list of entries.
        """
        taxonomy_entries = taxonomy.entries(self.taxonomy_entries)
        for a, b in zip(taxonomy_entries, self.taxonomy_entries):
            self.assertEqual(a, b)

    def test_write(self):
        """
        Write entries to a file-like object.
        """
        taxonomy_file = io.StringIO()
        taxonomy.write(taxonomy_file, self.taxonomy_entries)
        self.assertEqual(taxonomy_file.getvalue().rstrip(), TAXONOMY_SAMPLE.rstrip())


class TestTaxonomyEntry(unittest.TestCase):
    def setUp(self):
        """
        Creating taxonomy entries from a file-like object.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        self.taxonomy_lines = TAXONOMY_SAMPLE.split('\n')
        self.taxonomy_entries = list(taxonomy.entries(taxonomy_file))

    def test_length(self):
        """
        Check that the correct number of entries were read.
        """
        self.assertEqual(len(self.taxonomy_entries), 8)

    def test_entry_sequence_id(self):
        """
        Check the sequence ID of each entry.
        """
        for line, entry in zip(self.taxonomy_lines, self.taxonomy_entries):
            self.assertEqual(entry.identifier, line.split('\t')[0])

    def test_entry_taxonomy(self):
        """
        Check the taxonomy label of each entry.
        """
        for line, entry in zip(self.taxonomy_lines, self.taxonomy_entries):
            self.assertEqual(entry.label, line.split('\t')[1])

    def test_entry_serialization(self):
        """
        Serialize a taxonomy entry into a byte string.
        """
        self.assertIsInstance(self.taxonomy_entries[0].serialize(), bytes)

    def test_entry_deserialization(self):
        """
        Deserialize a taxonomy entry from a byte string.
        """
        deserialized = taxonomy.TaxonomyEntry.deserialize(self.taxonomy_entries[0].serialize())
        self.assertEqual(deserialized, self.taxonomy_entries[0])


class TestTaxonomyDb(unittest.TestCase):
    def setUp(self):
        """
        Create a taxonomy database from a file-like object.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        with taxonomy.TaxonomyDbFactory("/tmp/test.tax.db") as factory:
            self.taxonomy_entries = list(taxonomy.entries(taxonomy_file))
            factory.write_entries(self.taxonomy_entries)
        self.taxonomy_db = taxonomy.TaxonomyDb(factory.path)
        self.unique_labels = []
        for entry in self.taxonomy_entries:
            if entry.label not in self.unique_labels:
                self.unique_labels.append(entry.label)

    def tearDown(self):
        """
        Remove the taxonomy database.
        """
        self.taxonomy_db.close()
        for f in self.taxonomy_db.path.iterdir():
            f.unlink()
        self.taxonomy_db.path.rmdir()

    def test_length(self):
        """
        Check that the correct number of entries were inserted.
        """
        # 7 unique label entries
        self.assertEqual(len(self.taxonomy_db), 7)

    def test_contains_fasta_id(self):
        """
        Check that the taxonomy database contains a FASTA identifier.
        """
        self.assertTrue(self.taxonomy_db.contains_fasta_id(self.taxonomy_entries[0].identifier))

    def test_contains_label(self):
        """
        Check that the taxonomy database contains a label.
        """
        self.assertTrue(self.taxonomy_db.contains_label(self.taxonomy_entries[0].label))

    def test_count(self):
        """
        Check the sequence count by label.
        """
        self.assertEqual(self.taxonomy_db.count(0), 1)
        self.assertEqual(self.taxonomy_db.count(4), 2) # 2 sequences, same label

    def test_counts(self):
        """
        Check the counts of the taxonomy database.
        """
        self.assertEqual(list(self.taxonomy_db.counts()), [1, 1, 1, 1, 2, 1, 1])

    def test_fasta_id(self):
        """
        Check that the taxonomy database maps a label index to its corresponding FASTA identifier.
        """
        self.assertEqual(self.taxonomy_db.fasta_id_with_label(0, 0), self.taxonomy_entries[0].identifier)
        self.assertEqual(self.taxonomy_db.fasta_id_with_label(4, 0), self.taxonomy_entries[4].identifier)
        self.assertEqual(self.taxonomy_db.fasta_id_with_label(4, 1), self.taxonomy_entries[7].identifier)

    def test_fasta_id_to_index(self):
        """
        Check that the taxonomy database maps a FASTA identifier to its corresponding label index.
        """
        self.assertEqual(self.taxonomy_db.fasta_id_to_index(self.taxonomy_entries[0].identifier), 0)
        self.assertEqual(self.taxonomy_db.fasta_id_to_index(self.taxonomy_entries[7].identifier), 4)

    def test_fasta_id_to_label(self):
        """
        Check that the taxonomy database maps a FASTA identifier to its corresponding label.
        """
        self.assertEqual(
            self.taxonomy_db.fasta_id_to_label(self.taxonomy_entries[0].identifier),
            self.taxonomy_entries[0].label)
        self.assertEqual(
            self.taxonomy_db.fasta_id_to_label(self.taxonomy_entries[7].identifier),
            self.taxonomy_entries[4].label)

    def test_label(self):
        """
        Check that the taxonomy database maps a label index to its corresponding label.
        """
        self.assertEqual(self.taxonomy_db.label(0), self.taxonomy_entries[0].label)
        self.assertEqual(self.taxonomy_db.label(4), self.taxonomy_entries[4].label)

    def test_labels(self):
        """
        Check that the taxonomy database maps a label index to its corresponding label.
        """
        self.assertEqual(list(self.taxonomy_db.labels()), self.unique_labels)

    def test_label_to_index(self):
        """
        Check that the taxonomy database maps a label to its corresponding index.
        """
        self.assertEqual(self.taxonomy_db.label_to_index(self.taxonomy_entries[0].label), 0)
        self.assertEqual(self.taxonomy_db.label_to_index(self.taxonomy_entries[4].label), 4)
        self.assertEqual(self.taxonomy_db.label_to_index(self.taxonomy_entries[7].label), 4) # same label


class TestTaxon(unittest.TestCase):
    def setUp(self):
        self.a = taxonomy.Taxon(rank=0, name="A", parent=None)
        self.b = taxonomy.Taxon(rank=1, name="B", parent=self.a)

    def test_taxon_parent(self):
        self.assertIsNone(self.a.parent)
        self.assertIs(self.b.parent, self.a)

    def test_taxon_has_child(self):
        self.assertIn(self.b, self.a)
        self.assertIn(self.b.name, self.a)

    def test_taxon_is_child(self):
        self.assertIs(self.a[self.b.name], self.b)
        self.assertIs(self.a[self.b], self.b)


class TestTaxonomyHierarchy(unittest.TestCase):
    def setUp(self):
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        self.taxonomy_lines = TAXONOMY_SAMPLE.split('\n')
        self.taxonomy_entries = list(taxonomy.entries(taxonomy_file))
        self.hierarchy = taxonomy.TaxonomyHierarchy(depth=7)
        self.hierarchy.add_entries(self.taxonomy_entries)
        self.invalid_label = "k__Bacteria; p__Proteobacteria; c__XYZ; o__Acetobacterales; f__; g__; s__"
        self.test_label = "k__Bacteria; p__Proteobacteria; c__XYZ; o__; f__; g__; s__"

    def test_taxon_depth(self):
        self.assertEqual(self.hierarchy.depth, 7)

    def test_taxon_counts(self):
        self.assertEqual(len(self.hierarchy.taxon_to_id_map[0]), 2, "Kingdom")
        self.assertEqual(len(self.hierarchy.taxon_to_id_map[1]), 4, "Phylum")
        self.assertEqual(len(self.hierarchy.taxon_to_id_map[2]), 3, "Class")
        self.assertEqual(len(self.hierarchy.taxon_to_id_map[3]), 3, "Order")
        self.assertEqual(len(self.hierarchy.taxon_to_id_map[4]), 2, "Family")
        self.assertEqual(len(self.hierarchy.taxon_to_id_map[5]), 2, "Genus")
        self.assertEqual(len(self.hierarchy.taxon_to_id_map[6]), 1, "Species")

    def test_has_taxonomy(self):
        for entry in self.taxonomy_entries:
            self.assertTrue(self.hierarchy.has_entry(entry), "TaxonomyEntry")
            self.assertTrue(self.hierarchy.has_taxonomy(entry.label), "taxonomy string")
        self.assertFalse(self.hierarchy.has_taxonomy(self.invalid_label))

    def test_reduce_taxons(self):
        taxons = taxonomy.split_taxonomy(self.invalid_label)
        self.assertEqual(self.hierarchy.reduce_taxons(taxons), ("Bacteria", "Proteobacteria"))

    def test_reduce_taxonomy(self):
        self.assertEqual(self.hierarchy.reduce_taxonomy(self.test_label), "k__Bacteria; p__Proteobacteria; c__; o__; f__; g__; s__")
        self.assertEqual(self.hierarchy.reduce_taxonomy(self.invalid_label), "k__Bacteria; p__Proteobacteria; c__; o__; f__; g__; s__")

    def test_rebuild_taxon_id_maps_on_add(self):
        a = self.hierarchy.taxon_to_id_map
        b = self.hierarchy.id_to_taxon_map
        self.hierarchy.add_taxonomy(self.invalid_label)
        self.assertIsNot(self.hierarchy.taxon_to_id_map, a)
        self.assertIsNot(self.hierarchy.id_to_taxon_map, b)

    def test_prevent_rebuild_taxon_id_maps_on_add(self):
        a = self.hierarchy.taxon_to_id_map
        b = self.hierarchy.id_to_taxon_map
        self.hierarchy.add_entry(self.taxonomy_entries[3])
        self.assertIs(self.hierarchy.taxon_to_id_map, a)
        self.assertIs(self.hierarchy.id_to_taxon_map, b)

    def test_taxon_id_map_order(self):
        # Ensure that the taxons in each rank are grouped together.
        for rank, taxon_group in enumerate(self.hierarchy.id_to_taxon_map[:-1]):
            offset = 0
            for parent in taxon_group:
                j = 0
                for j in range(len(parent.children)):
                    self.assertIs(self.hierarchy.id_to_taxon_map[rank+1][offset+j].parent, parent)
                offset += len(parent.children)

    def test_tokenize(self):
        tokens = self.hierarchy.tokenize(self.taxonomy_entries[3].label)
        self.assertEqual(len(tokens), 3)

    def test_tokenize_with_truncated_hierarchy(self):
        truncated_hierarchy = taxonomy.TaxonomyHierarchy.merged([self.hierarchy], 2)
        tokens = truncated_hierarchy.tokenize(self.taxonomy_entries[3].label)
        self.assertEqual(len(tokens), 2)

    def test_tokenize_with_pad(self):
        tokens = self.hierarchy.tokenize(self.taxonomy_entries[3].label, pad=True)
        self.assertEqual(len(tokens), 7)
        self.assertEqual(tuple(tokens[3:]), (-1, -1, -1, -1))

    def test_tokenize_with_pad_and_include_missing(self):
        tokens = self.hierarchy.tokenize(self.taxonomy_entries[3].label, pad=True, include_missing=True)
        self.assertEqual(len(tokens), 7)
        self.assertEqual(tuple(tokens[3:]), (0, 0, 0, 0))

    def test_detokenize(self):
        for entry in self.taxonomy_entries:
            tokens = self.hierarchy.tokenize(entry.label, pad=False)
            self.assertEqual(self.hierarchy.detokenize(tokens), entry.label)

    def test_detokenize_with_pad(self):
        for entry in self.taxonomy_entries:
            tokens = self.hierarchy.tokenize(entry.label, pad=True)
            self.assertEqual(self.hierarchy.detokenize(tokens), entry.label)

    def test_detokenize_with_pad_and_include_missing(self):
        for entry in self.taxonomy_entries:
            tokens = self.hierarchy.tokenize(entry.label, pad=True, include_missing=True)
            self.assertEqual(self.hierarchy.detokenize(tokens, include_missing=True), entry.label)

    def test_serialize_deserialize(self):
        serialized = self.hierarchy.serialize()
        deserialized = taxonomy.TaxonomyHierarchy.deserialize(serialized)
        self.assertEqual(self.hierarchy, deserialized)

class TestTaxonomyHierarchyMerge(unittest.TestCase):
    def setUp(self):
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        self.taxonomy_lines = TAXONOMY_SAMPLE.split('\n')
        self.taxonomy_entries = list(taxonomy.entries(taxonomy_file))
        self.hierarchy = taxonomy.TaxonomyHierarchy(depth=7)
        self.hierarchy.add_entries(self.taxonomy_entries)
        self.a = taxonomy.TaxonomyHierarchy(depth=7)
        self.a.add_entries(self.taxonomy_entries[:4])
        self.b = taxonomy.TaxonomyHierarchy(depth=7)
        self.b.add_entries(self.taxonomy_entries[4:])

    def test_merge_hierarchy_max_depth(self):
        merged = taxonomy.TaxonomyHierarchy.merged([self.a, self.b])
        for taxon_a, taxon_b in zip(self.hierarchy, merged):
            self.assertEqual(taxon_a, taxon_b)

    def test_merge_hierarchy_truncated_depth(self):
        depth = 4
        merged = taxonomy.TaxonomyHierarchy.merged([self.a, self.b], depth=depth)
        for taxon_a, taxon_b in zip(self.hierarchy, merged):
            self.assertEqual(taxon_a, taxon_b)


if __name__ == "__main__":
    unittest.main()
