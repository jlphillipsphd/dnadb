import io
import os
import sys
import unittest

sys.path.append("./src")

from dnadb import fasta, taxonomy

# Sampled from SILVA 138.1

FASTA_SAMPLE = """\
>AY855839.1.1390
GCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGCCTTGTAGTCCGTGAGGCGGCGGACGGGTGAGTAACACGTGGGCA
>FW343016.1.1511
ACCAGCGGCGGCGTGCTTAACACATGCAAGTCGAACGGCCTTGTAGTCCGTGAGGCGGCGGACGGGTGAGTAACACGTGG
>AY835431.189876.191345
GCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGCCTTGTAGTCCGTGAGGCGGCGGACGGGTGAGTAACACGTGGGCA
>FW369114.1.1462
GGCCAACGCGTGAGTGATGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGGGGATGACGTCAAATCAT
>FW369795.1.1413
ACTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACGTCTTCGG
>AY846383.1.1790
TTGCAACGGTGGAGCATGTGGTTTAATTCGAAGCAACCGGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAG
>AB001440.1.1538
GATCCAGCCATCCGCGTGTCTTCATCGATGAAGAACGCAGCGAAATGCGTGGGGAGCAAACAGGATTAGATACCCTGGTA
>FW369795.1.xxxx
AAAGGTTCAGTGGCGCAAGGCTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCG
"""

TAXONOMY_SAMPLE = """\
AY855839.1.1390\td__Bacteria;p__;c__;o__;f__;g__;s__
FW343016.1.1511\td__Bacteria;p__Firmicutes;c__;o__;f__;g__;s__
AY835431.189876.191345\td__Bacteria;p__Cyanobacteria;c__;o__;f__;g__;s__
FW369114.1.1462\td__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__;f__;g__;s__
FW369795.1.1413\td__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Acetobacterales;f__;g__;s__
AY846383.1.1790\td__Eukaryota;p__Eukaryota;c__Chlorophyceae;o__Sphaeropleales;f__Sphaeropleales;g__Monoraphidium;s__
AB001440.1.1538\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__test_species
FW369795.1.xxxx\td__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Acetobacterales;f__;g__;s__
"""

TAXONOMY_SAMPLE_WITH_HEADER = "Sequence ID\tTaxonomy\n" + TAXONOMY_SAMPLE

class TestIsTaxonomy(unittest.TestCase):
    def test_is_taxonomy(self):
        self.assertFalse(taxonomy.is_taxonomy("Taxonomy"))
        self.assertTrue(taxonomy.is_taxonomy("d__Bacteria;p__;c__;o__;f__;g__;s__;"))
        self.assertTrue(taxonomy.is_taxonomy("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__test_species"))


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
            taxonomy.join_taxonomy(["Bacteria"]), "d__Bacteria"
        )
        self.assertEqual(
            taxonomy.join_taxonomy(["Bacteria", ""]), "d__Bacteria;p__"
        )
        self.assertEqual(
            taxonomy.join_taxonomy(["Bacteria"], depth=2), "d__Bacteria;p__"
        )

class TestTaxonomyReading(unittest.TestCase):
    def setUp(self):
        """
        Write the taxonomy sample to a file for the tests to read from.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        self.taxonomy_entries = list(taxonomy.read(taxonomy_file, header=False))
        with open("/tmp/test.tsv", 'w') as f:
            f.write(TAXONOMY_SAMPLE)

    def tearDown(self):
        """
        Remove the taxonomy database.
        """
        os.remove("/tmp/test.tsv")

    def test_read_without_header(self):
        """
        Read taxonomy entries from a file-like object.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        taxonomy_entries = list(taxonomy.read(taxonomy_file, header=False))
        self.assertEqual(taxonomy_entries[0].sequence_id, "AY855839.1.1390")
        self.assertEqual(taxonomy_entries[0].label, "d__Bacteria;p__;c__;o__;f__;g__;s__")
        self.assertEqual(taxonomy_entries[1].sequence_id, "FW343016.1.1511")
        self.assertEqual(taxonomy_entries[1].label, "d__Bacteria;p__Firmicutes;c__;o__;f__;g__;s__")

    def test_read_without_header_auto(self):
        """
        Read taxonomy entries from a file-like object.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        taxonomy_entries = list(taxonomy.read(taxonomy_file, header="auto"))
        self.assertEqual(taxonomy_entries[0].sequence_id, "AY855839.1.1390")
        self.assertEqual(taxonomy_entries[0].label, "d__Bacteria;p__;c__;o__;f__;g__;s__")
        self.assertEqual(taxonomy_entries[1].sequence_id, "FW343016.1.1511")
        self.assertEqual(taxonomy_entries[1].label, "d__Bacteria;p__Firmicutes;c__;o__;f__;g__;s__")

    def test_read_with_header(self):
        """
        Read taxonomy entries from a file-like object.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE_WITH_HEADER)
        taxonomy_entries = list(taxonomy.read(taxonomy_file, header=True))
        self.assertEqual(taxonomy_entries[0].sequence_id, "AY855839.1.1390")
        self.assertEqual(taxonomy_entries[0].label, "d__Bacteria;p__;c__;o__;f__;g__;s__")
        self.assertEqual(taxonomy_entries[1].sequence_id, "FW343016.1.1511")
        self.assertEqual(taxonomy_entries[1].label, "d__Bacteria;p__Firmicutes;c__;o__;f__;g__;s__")

    def test_read_with_header_auto(self):
        """
        Read taxonomy entries from a file-like object.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE_WITH_HEADER)
        taxonomy_entries = list(taxonomy.read(taxonomy_file, header="auto"))
        self.assertEqual(taxonomy_entries[0].sequence_id, "AY855839.1.1390")
        self.assertEqual(taxonomy_entries[0].label, "d__Bacteria;p__;c__;o__;f__;g__;s__")
        self.assertEqual(taxonomy_entries[1].sequence_id, "FW343016.1.1511")
        self.assertEqual(taxonomy_entries[1].label, "d__Bacteria;p__Firmicutes;c__;o__;f__;g__;s__")

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
        self.assertEqual(taxonomy_file.getvalue().strip(), TAXONOMY_SAMPLE.strip())


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
            self.assertEqual(entry.sequence_id, line.split('\t')[0])

    def test_entry_taxonomy(self):
        """
        Check the taxonomy label of each entry.
        """
        for line, entry in zip(self.taxonomy_lines, self.taxonomy_entries):
            self.assertEqual(entry.label, line.split('\t')[1])


class TestTaxonomyTree(unittest.TestCase):
    def setUp(self):
        """
        Create a taxonomy tree from a file-like object.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        self.taxonomy_lines = TAXONOMY_SAMPLE.split('\n')
        self.taxonomy_entries = list(taxonomy.entries(taxonomy_file))
        factory = taxonomy.TaxonomyTreeFactory(depth=7)
        factory.add_entries(self.taxonomy_entries)
        self.tree = factory.build()
        self.invalid_label = "d__Bacteria;p__Proteobacteria;c__XYZ;o__Acetobacterales;f__;g__;s__"
        self.test_label = "d__Bacteria;p__Proteobacteria;c__XYZ;o__;f__;g__;s__"

    def test_taxon_depth(self):
        self.assertEqual(self.tree.depth, 7)

    def test_taxon_counts(self):
        self.assertEqual(len(self.tree.id_to_taxon_map[0]), 2, "Kingdom")
        self.assertEqual(len(self.tree.id_to_taxon_map[1]), 5, "Phylum")
        self.assertEqual(len(self.tree.id_to_taxon_map[3]), 4, "Class")
        self.assertEqual(len(self.tree.id_to_taxon_map[4]), 3, "Family")
        self.assertEqual(len(self.tree.id_to_taxon_map[5]), 3, "Genus")
        self.assertEqual(len(self.tree.id_to_taxon_map[6]), 2, "Species")

    def test_taxonomy_counts(self):
        self.assertEqual(len(self.tree.taxonomy_id_map[0]), 2, "Kingdom")
        self.assertEqual(len(self.tree.taxonomy_id_map[1]), 5, "Phylum")
        self.assertEqual(len(self.tree.taxonomy_id_map[2]), 6, "Class")
        self.assertEqual(len(self.tree.taxonomy_id_map[3]), 7, "Order")
        self.assertEqual(len(self.tree.taxonomy_id_map[4]), 7, "Family")
        self.assertEqual(len(self.tree.taxonomy_id_map[5]), 7, "Genus")
        self.assertEqual(len(self.tree.taxonomy_id_map[6]), 7, "Species")

    def test_bidirectional_taxon_mapping(self):
        for depth in range(self.tree.depth):
            for i, taxon in enumerate(self.tree.id_to_taxon_map[depth]):
                self.assertEqual(i, self.tree.taxon_to_id_map[depth][taxon])

    def test_taxon_sorted_order(self):
        for depth in range(self.tree.depth):
            for t1, t2 in zip(self.tree.id_to_taxon_map[depth], sorted(self.tree.id_to_taxon_map[depth])):
                self.assertEqual(t1, t2)

    def test_taxonomy_sorted_order(self):
        for depth in range(self.tree.depth):
            for t1, t2 in zip(self.tree.taxonomy_id_map[depth], sorted(self.tree.taxonomy_id_map[depth], key=lambda x: x.taxonomy_id)):
                self.assertEqual(t1, t2)

    def test_taxonomy_parent_sorted_order(self):
        for depth in range(1, self.tree.depth):
            parents = [t.parent.taxonomy_id for t in self.tree.taxonomy_id_map[depth]]
            self.assertEqual(parents, sorted(parents))

    def test_has_taxonomy(self):
        for entry in self.taxonomy_entries:
            self.assertTrue(self.tree.has_taxonomy(entry), "TaxonomyEntry")
            self.assertIn(entry.label, self.tree, "taxonomy string")
        self.assertFalse(self.tree.has_taxonomy(self.invalid_label))

    def test_get_taxonomy(self):
        for entry in self.taxonomy_entries:
            self.assertEqual(self.tree.taxonomy(entry).taxonomy_label, entry.label, "TaxonomyEntry")
            self.assertEqual(self.tree[entry.label].taxonomy_label, entry.label, "taxonomy string")
            self.assertEqual(self.tree[self.tree[entry].taxonomy_id].taxonomy_label, entry.label, "taxonomy ID lookup")
        self.assertRaises(KeyError, self.tree.taxonomy, self.invalid_label)

    def test_iter(self):
        for t1, t2 in zip(self.tree, self.tree.taxonomy_id_map[self.tree.depth - 1]):
            self.assertEqual(t1, t2)

    def test_reduce_entry(self):
        self.assertEqual(self.tree.reduce_entry(self.taxonomy_entries[0]).label, "d__Bacteria;p__;c__;o__;f__;g__;s__")
        self.assertEqual(self.tree.reduce_entry(self.taxonomy_entries[6]).label, "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__test_species")

    def test_reduce_taxonomy(self):
        self.assertEqual(self.tree.reduce_label(self.test_label), "d__Bacteria;p__Proteobacteria;c__;o__;f__;g__;s__")
        self.assertEqual(self.tree.reduce_label(self.invalid_label), "d__Bacteria;p__Proteobacteria;c__;o__;f__;g__;s__")

    def test_num_taxonomies(self):
        self.assertEqual(self.tree.taxonomy_id_map[0][0].num_taxonomies, 6)
        self.assertEqual(self.tree.taxonomy_id_map[0][1].num_taxonomies, 1)

    def test_taxonomy_range(self):
        # A better test case could be written for this...
        self.assertEqual(self.tree.taxonomy_id_map[0][0].taxonomy_id_range, range(0, 6))
        self.assertEqual(self.tree.taxonomy_id_map[0][1].taxonomy_id_range, range(6, 7))


class TestTaxonomyDb(unittest.TestCase):
    def setUp(self):
        """
        Create a taxonomy database from a file-like object.
        """
        taxonomy_file = io.StringIO(TAXONOMY_SAMPLE)
        with fasta.FastaDbFactory("/tmp/test.fasta.db") as factory:
            factory.write_entries(fasta.entries(io.StringIO(FASTA_SAMPLE)))
        self.fasta_db = fasta.FastaDb(factory.path)
        with taxonomy.TaxonomyDbFactory("/tmp/test.tax.db", self.fasta_db, 7) as factory:
            self.taxonomy_entries = list(taxonomy.entries(taxonomy_file))
            factory.write_entries(self.taxonomy_entries)
        self.taxonomy_db = taxonomy.TaxonomyDb(factory.path)
        self.unique_labels = []
        for entry in self.taxonomy_entries:
            if entry.label not in self.unique_labels:
                self.unique_labels.append(entry.label)
        self.invalid_label = "d__Bacteria;p__Proteobacteria;c__XYZ;o__Acetobacterales;f__;g__;s__"
        self.test_label = "d__Bacteria;p__Proteobacteria;c__XYZ;o__;f__;g__;s__"

    def tearDown(self):
        """
        Remove the taxonomy database.
        """
        self.fasta_db.close()
        for f in self.fasta_db.path.iterdir():
            f.unlink()
        self.fasta_db.path.rmdir()
        self.taxonomy_db.close()
        for f in self.taxonomy_db.path.iterdir():
            f.unlink()
        self.taxonomy_db.path.rmdir()

    def test_num_sequences(self):
        """
        Check that the correct number of entries were inserted.
        """
        self.assertEqual(len(self.taxonomy_db), len(self.fasta_db))

    def test_num_labels(self):
        """
        Check that the number of labels in the taxonomoy DB is correct.
        """
        self.assertEqual(self.taxonomy_db.num_labels, len(self.unique_labels))

    def test_contains_fasta_id(self):
        """
        Check that the taxonomy database contains a FASTA sequence_id.
        """
        for entry in self.taxonomy_entries:
            self.assertIn(entry.sequence_id, self.taxonomy_db)

    def test_contains_label(self):
        """
        Check that the taxonomy database contains a label.
        """
        for entry in self.taxonomy_entries:
            self.assertTrue(self.taxonomy_db.has_taxonomy(entry.label))

    def test_iteration(self):
        """
        Check that the taxonomy database can be iterated over.
        """
        n = 0
        for entry in self.taxonomy_db:
            n += 1
            self.assertIn(entry.sequence_id, self.fasta_db)
        self.assertEqual(n, len(self.fasta_db))

    def test_sequence_index_to_id(self):
        """
        Check that the sequence index to ID mapping is correct.
        """
        for i, entry in enumerate(self.taxonomy_db):
            self.assertEqual(self.taxonomy_db.sequence_index_to_id(i), entry.sequence_id)

    def test_sequence_id_to_index(self):
        """
        Check that the sequence ID to index mapping is correct.
        """
        for i, entry in enumerate(self.taxonomy_db):
            self.assertEqual(self.taxonomy_db.sequence_id_to_index(entry.sequence_id), i)

    def test_sequence_id_to_label(self):
        """
        Check that the sequence ID to label mapping is correct.
        """
        for entry in self.taxonomy_entries:
            self.assertEqual(self.taxonomy_db[entry.sequence_id].taxonomy.taxonomy_label, entry.label)

    def test_sequence_index_to_label(self):
        """
        Check that the sequence index to label mapping is correct.
        """
        for i, entry in enumerate(self.taxonomy_entries):
            self.assertEqual(self.taxonomy_db[i].taxonomy.taxonomy_label, entry.label)

    def test_sequences_with_label(self):
        """
        Check that the sequences with label method returns the correct number of sequences.
        """
        n = 0
        for label in self.unique_labels:
           entries = list(self.taxonomy_db.sequences_with_taxonomy(label))
           for entry in entries:
                self.assertEqual(entry.taxonomy.taxonomy_label, label)
                n += 1
        self.assertEqual(n, len(self.fasta_db))

    def test_count(self):
        """
        Check that the count method returns the correct number of sequences.
        """
        for label in self.unique_labels:
            entries = list(self.taxonomy_db.sequences_with_taxonomy(label))
            self.assertEqual(len(entries), self.taxonomy_db.count(label))


if __name__ == "__main__":
    unittest.main()
