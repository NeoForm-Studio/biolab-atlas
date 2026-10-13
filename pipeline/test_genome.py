import json
import shutil
import unittest
from pathlib import Path

from genome import best_hits, composition, find_markers, parse_fasta, passes

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

    def test_best_hits_keeps_highest_bitscore(self):
        tsv = "q1\ts1\t40.0\t90\t1e-20\t100\nq1\ts2\t90.0\t95\t1e-50\t300\n"
        self.assertEqual(best_hits(tsv)["q1"]["accession"], "s2")

    def test_passes_applies_coverage_identity_and_per_gene_override(self):
        hit = {"evalue": 1e-30, "coverage": 98, "identity": 39.0}
        self.assertTrue(passes(hit, {}))
        self.assertFalse(passes(hit, {"min_identity": 45}))  # look-alike aminotransferase vs EctB
        self.assertFalse(passes({**hit, "coverage": 47}, {}))  # partial-length family hit
        self.assertFalse(passes({**hit, "evalue": 1e-6}, {}))
        self.assertFalse(passes(None, {}))

    @unittest.skipUnless(shutil.which("blastp"), "BLAST+ not installed")
    def test_markers_find_themselves(self):
        faa = Path(__file__).parent / "markers.faa"
        m = find_markers(parse_fasta(faa.read_text()), PANELS, faa)
        self.assertEqual({t: v["found"] for t, v in m.items()}, {t: len(p["genes"]) for t, p in PANELS.items()})

    @unittest.skipUnless(shutil.which("blastp"), "BLAST+ not installed")
    def test_unrelated_proteome_finds_nothing(self):
        m = find_markers(self.records, PANELS, Path(__file__).parent / "markers.faa")
        self.assertEqual(sum(v["found"] for v in m.values()), 0)


class TemperatureRangeTest(unittest.TestCase):
    def test_parses_single_values_ranges_and_negatives(self):
        from atlas import temperature_range
        self.assertEqual(temperature_range("26"), (26.0, 26.0))
        self.assertEqual(temperature_range("30-37"), (30.0, 37.0))
        self.assertEqual(temperature_range("-5-10"), (-5.0, 10.0))
        self.assertIsNone(temperature_range("mesophilic"))


class ScoreTest(unittest.TestCase):
    def test_percentile_and_auc(self):
        from score import auc, percentile
        self.assertEqual(percentile([10, 20, 30, 40], 30), 62)
        self.assertEqual(auc([3, 4], [1, 2]), 1.0)
        self.assertEqual(auc([1], [1]), 0.5)
        self.assertIsNone(auc([], [1]))

    def test_score_library_validates_heat_against_measured_temperatures(self):
        from score import score_library
        rows = [{"predicted_ogt_c": t, "acidic_excess": 0.01, "markers": {"salt": "0/3", "radiation": "0/6"},
                 "measured_growth_c": [m, m]} for t, m in [(95, 90), (80, 75), (35, 37), (30, 30), (15, 10)]]
        v = score_library(rows)
        self.assertEqual(v["heat"]["auc"], 1.0)
        self.assertEqual((v["heat"]["positives"], v["cold"]["positives"]), (2, 1))
        self.assertEqual(rows[0]["scores"]["heat"], 90)


class NoGenomeTest(unittest.TestCase):
    def test_species_without_assembly_reports_relatives(self):
        from unittest import mock
        import atlas
        taxon = {"tax_id": 9, "organism_name": "Testus nullus", "rank": "SPECIES", "lineage": [1, 2, 8]}
        for name, value in [("resolve", taxon), ("inventory", {"assemblies": 0}), ("pick_assembly", None), ("phenotype", None)]:
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
