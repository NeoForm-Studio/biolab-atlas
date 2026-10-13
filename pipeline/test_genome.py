import json
import unittest
from pathlib import Path

from genome import composition, find_markers, parse_fasta, product

PANELS = json.loads((Path(__file__).parent / "markers.json").read_text(encoding="utf-8"))

FASTA = """>WP_1.1 MULTISPECIES: recombinase RecA [Deinococcus]
MDEK
>WP_2.1 DNA repair protein PprA [Deinococcus radiodurans]
IVYW
REL*
>WP_3.1 diaminobutyrate--2-oxoglutarate transaminase family protein [Halobacterium]
KK
"""


class GenomeTest(unittest.TestCase):
    def setUp(self):
        self.records = parse_fasta(FASTA)

    def test_parse_joins_wrapped_lines(self):
        self.assertEqual([s for _, s in self.records], ["MDEK", "IVYWREL*", "KK"])

    def test_product_strips_accession_multispecies_and_organism(self):
        self.assertEqual(product(self.records[0][0]), "recombinase RecA")

    def test_composition(self):
        c = composition(self.records)
        # 13 residues (stop removed): IVYWREL letters = I,V,Y,W,R,E,L + E in MDEK = 8
        self.assertEqual(c["proteins"], 3)
        self.assertAlmostEqual(c["ivywrel"], round(8 / 13, 4))
        self.assertEqual(c["predicted_ogt_c"], round(937 * 8 / 13 - 335))
        # D,E = 3 (D,E in MDEK; E in REL); K,R = 4 (K in MDEK, R, K, K)
        self.assertAlmostEqual(c["acidic_excess"], round(3 / 13 - 4 / 13, 4))

    def test_empty_proteome(self):
        self.assertIsNone(composition([]))

    def test_markers_exact_names_only(self):
        m = find_markers(self.records, PANELS)
        self.assertEqual(m["radiation"]["found"], 1)
        self.assertEqual(m["radiation"]["genes"][0]["accessions"], ["WP_2.1"])
        # a "... family protein" must not count as the real enzyme
        self.assertEqual(m["salt"]["found"], 0)


class NoGenomeTest(unittest.TestCase):
    def test_species_without_assembly_reports_relatives(self):
        from unittest import mock
        import atlas
        taxon = {"tax_id": 9, "organism_name": "Testus nullus", "rank": "SPECIES", "lineage": [1, 2, 8]}
        for name, value in [("resolve", taxon), ("inventory", {"assemblies": 0}), ("pick_assembly", None)]:
            p = mock.patch.object(atlas, name, return_value=value)
            p.start()
            self.addCleanup(p.stop)
        with mock.patch.object(atlas, "relatives", return_value=[{"name": "Testus plenus"}]) as rel:
            e = atlas.profile("Testus nullus")
        self.assertEqual(e["status"], "no_genome")
        self.assertEqual(e["taxon"]["domain"], "Bacteria")
        rel.assert_called_once_with(8)  # relatives come from the parent taxon


if __name__ == "__main__":
    unittest.main()
