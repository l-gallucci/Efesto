"""Tests for src/efesto/clustering.py"""

from efesto.clustering import build_clusters, cluster_by_index


class TestCatalogModeSingletons:
    def test_catalog_mode_ignores_index_adjacency(self):
        # Headers share a stem (as colliding catalog IDs can after dereplication)
        # and are numerically adjacent — without catalog_mode these would merge.
        orf_hits = {f"contigA_{i}": {"hmm_stem": "x"} for i in range(1, 6)}
        clusters = build_clusters(
            "genome", orf_hits, orf_coords={}, max_gap=5, max_bp_gap=5000,
            strand_aware=False, catalog_mode=True)
        assert sorted(clusters) == sorted([[o] for o in orf_hits])

    def test_catalog_mode_ignores_real_coordinates_too(self):
        orf_hits = {"a": {"hmm_stem": "x"}, "b": {"hmm_stem": "y"}}
        orf_coords = {
            "a": {"contig": "c1", "start": 100, "end": 200, "strand": "+"},
            "b": {"contig": "c1", "start": 210, "end": 300, "strand": "+"},
        }
        clusters = build_clusters(
            "genome", orf_hits, orf_coords, max_gap=5, max_bp_gap=5000,
            strand_aware=False, catalog_mode=True)
        assert sorted(clusters) == [["a"], ["b"]]

    def test_non_catalog_mode_unaffected(self):
        # Same colliding-stem input, catalog_mode=False: original index
        # clustering behavior (still assumes header adjacency == proximity).
        orf_hits = {f"contigA_{i}": {"hmm_stem": "x"} for i in range(1, 6)}
        clusters = build_clusters(
            "genome", orf_hits, orf_coords={}, max_gap=5, max_bp_gap=5000,
            strand_aware=False, catalog_mode=False)
        assert len(clusters) == 1
        assert len(clusters[0]) == 5


class TestClusterByIndex:
    def test_gap_splits_cluster(self):
        orfs = {"c_1", "c_2", "c_8"}
        clusters = cluster_by_index(orfs, max_gap=5)
        stems = sorted(sorted(c) for c in clusters)
        assert stems == [["c_1", "c_2"], ["c_8"]]
