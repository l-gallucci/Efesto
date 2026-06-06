"""Tests for src/efesto/scoring.py"""

import math
import pytest

from efesto.scoring import (
    cluster_confidence,
    co_occurrence_score,
    distance_decay,
    hmm_weight,
    uniop_pair_score,
)


class TestDistanceDecay:
    def test_zero_gap(self):
        assert distance_decay(0) == pytest.approx(1.0)

    def test_half_life(self):
        assert distance_decay(500) == pytest.approx(0.5, rel=1e-6)

    def test_double_half_life(self):
        assert distance_decay(1000) == pytest.approx(0.25, rel=1e-6)

    def test_negative_gap_clamped(self):
        assert distance_decay(-100) == pytest.approx(1.0)

    def test_custom_half_life(self):
        assert distance_decay(1000, half_life=1000) == pytest.approx(0.5, rel=1e-6)

    def test_large_gap_near_zero(self):
        assert distance_decay(10000) < 0.001


class TestCoOccurrenceScore:
    def test_single_gene_always_one(self):
        coords = {"a": {"start": 100, "end": 300, "contig": "c1"}}
        assert co_occurrence_score(["a"], coords, {"c1": 50000}) == 1.0

    def test_single_gene_with_canonical_size(self):
        coords = {"a": {"start": 100, "end": 300, "contig": "c1"}}
        assert co_occurrence_score(["a"], coords, {"c1": 50000}, canonical_size=3) == 1.0

    def test_empty_coords_no_penalty(self):
        s = co_occurrence_score(["a", "b"], {}, {})
        assert s == pytest.approx(1.0)

    def test_close_genes_no_edge(self):
        coords = {
            "a": {"start": 5000, "end": 5300, "contig": "c1"},
            "b": {"start": 5350, "end": 5650, "contig": "c1"},
        }
        s = co_occurrence_score(["a", "b"], coords, {"c1": 50000})
        assert 0.9 < s <= 1.0

    def test_distant_genes_penalized(self):
        coords = {
            "a": {"start": 5000, "end": 5300, "contig": "c1"},
            "b": {"start": 15000, "end": 15300, "contig": "c1"},
        }
        s = co_occurrence_score(["a", "b"], coords, {"c1": 50000})
        assert s < 0.1

    def test_edge_penalty_applied(self):
        coords = {
            "a": {"start": 100, "end": 400, "contig": "c1"},
            "b": {"start": 500, "end": 800, "contig": "c1"},
        }
        s_edge    = co_occurrence_score(["a", "b"], coords, {"c1": 50000})
        coords_mid = {
            "a": {"start": 10000, "end": 10300, "contig": "c1"},
            "b": {"start": 10400, "end": 10700, "contig": "c1"},
        }
        s_mid = co_occurrence_score(["a", "b"], coords_mid, {"c1": 50000})
        assert s_edge < s_mid

    def test_completeness_penalty(self):
        coords = {
            "a": {"start": 5000, "end": 5300, "contig": "c1"},
            "b": {"start": 5400, "end": 5700, "contig": "c1"},
        }
        s_full    = co_occurrence_score(["a", "b"], coords, {"c1": 50000}, canonical_size=2)
        s_partial = co_occurrence_score(["a", "b"], coords, {"c1": 50000}, canonical_size=4)
        assert s_full > s_partial

    def test_completeness_capped_at_one(self):
        coords = {
            "a": {"start": 5000, "end": 5300, "contig": "c1"},
            "b": {"start": 5400, "end": 5700, "contig": "c1"},
        }
        s = co_occurrence_score(["a", "b"], coords, {"c1": 50000}, canonical_size=1)
        assert s <= 1.0

    def test_score_in_unit_interval(self):
        coords = {
            "a": {"start": 100, "end": 400, "contig": "c1"},
            "b": {"start": 500, "end": 800, "contig": "c1"},
            "c": {"start": 900, "end": 1200, "contig": "c1"},
        }
        s = co_occurrence_score(list(coords), coords, {"c1": 5000}, canonical_size=5)
        assert 0.0 <= s <= 1.0


class TestUniOPPairScore:
    def test_empty_probs_neutral(self):
        assert uniop_pair_score(["a", "b"], {}) == 1.0

    def test_single_gene_neutral(self):
        assert uniop_pair_score(["a"], {("a", "b"): 0.3}) == 1.0

    def test_weakest_link(self):
        pp = {("a", "b"): 0.9, ("a", "c"): 0.6, ("b", "c"): 0.8}
        assert uniop_pair_score(["a", "b", "c"], pp) == pytest.approx(0.6)

    def test_reverse_key_lookup(self):
        pp = {("b", "a"): 0.7}
        assert uniop_pair_score(["a", "b"], pp) == pytest.approx(0.7)

    def test_partial_coverage_uses_available(self):
        pp = {("a", "b"): 0.8}
        s = uniop_pair_score(["a", "b", "c"], pp)
        assert s == pytest.approx(0.8)

    def test_no_matching_pairs_neutral(self):
        pp = {("x", "y"): 0.1}
        assert uniop_pair_score(["a", "b"], pp) == 1.0


class TestHmmWeight:
    def test_empty_cluster(self):
        assert hmm_weight([]) == 1.0

    def test_all_calibrated(self):
        rows = [{"confidence": "calibrated"}] * 4
        assert hmm_weight(rows) == pytest.approx(1.0)

    def test_all_low_confidence(self):
        rows = [{"confidence": "low_confidence"}] * 4
        assert hmm_weight(rows) == pytest.approx(0.5)

    def test_mixed(self):
        rows = [{"confidence": "calibrated"}, {"confidence": "low_confidence"}]
        assert hmm_weight(rows) == pytest.approx(0.75)

    def test_unknown_confidence_treated_as_low(self):
        rows = [{"confidence": "unknown_value"}]
        assert hmm_weight(rows) == pytest.approx(0.5)


class TestClusterConfidence:
    def test_all_ones(self):
        assert cluster_confidence(1.0, 1.0, 1.0) == pytest.approx(1.0)

    def test_multiplicative(self):
        assert cluster_confidence(0.8, 0.9, 0.7) == pytest.approx(0.8 * 0.9 * 0.7)

    def test_bgc_boost_applied(self):
        assert cluster_confidence(1.0, 1.0, 1.0, bgc_boost=1.2) == pytest.approx(1.2)

    def test_capped_at_1_2(self):
        assert cluster_confidence(1.0, 1.0, 1.0, bgc_boost=2.0) == pytest.approx(1.2)

    def test_default_bgc_boost_is_one(self):
        assert cluster_confidence(0.5, 0.6, 0.8) == pytest.approx(0.5 * 0.6 * 0.8)
