# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Tests for uniop.py — _parse_uniop_pred pair_probs return value."""

import pytest

from efesto.uniop import _parse_uniop_pred, _parse_uniop_operon


def _write_pred(tmp_path, lines):
    p = tmp_path / "uniop.pred"
    p.write_text("Gene A\tGene B\tPrediction\n" + "\n".join(lines) + "\n")
    return str(p)


def _idx(pairs):
    """Build idx_to_orf from {idx: name} dict."""
    return pairs


class TestParseUniopPred:
    def test_returns_tuple(self, tmp_path):
        p = _write_pred(tmp_path, ["1\t2\t0.9"])
        idx = {1: "orf1", 2: "orf2"}
        result = _parse_uniop_pred(str(p), idx)
        assert isinstance(result, tuple) and len(result) == 2

    def test_pair_probs_stored(self, tmp_path):
        p = _write_pred(tmp_path, ["1\t2\t0.85"])
        idx = {1: "orf1", 2: "orf2"}
        orf_to_op, pair_probs = _parse_uniop_pred(str(p), idx)
        assert len(pair_probs) == 1
        key = ("orf1", "orf2") if ("orf1", "orf2") in pair_probs else ("orf2", "orf1")
        assert pair_probs[key] == pytest.approx(0.85)

    def test_pair_probs_below_threshold_still_stored(self, tmp_path):
        p = _write_pred(tmp_path, ["1\t2\t0.2"])
        idx = {1: "orf1", 2: "orf2"}
        orf_to_op, pair_probs = _parse_uniop_pred(str(p), idx)
        assert len(pair_probs) == 1
        assert "orf1" not in orf_to_op

    def test_canonical_key_order(self, tmp_path):
        p = _write_pred(tmp_path, ["2\t1\t0.7"])
        idx = {1: "orf1", 2: "orf2"}
        _, pair_probs = _parse_uniop_pred(str(p), idx)
        assert ("orf1", "orf2") in pair_probs or ("orf2", "orf1") in pair_probs
        stored_prob = pair_probs.get(("orf1", "orf2")) or pair_probs.get(("orf2", "orf1"))
        assert stored_prob == pytest.approx(0.7)

    def test_operon_membership_above_threshold(self, tmp_path):
        p = _write_pred(tmp_path, ["1\t2\t0.9", "2\t3\t0.8"])
        idx = {1: "orf1", 2: "orf2", 3: "orf3"}
        orf_to_op, _ = _parse_uniop_pred(str(p), idx)
        assert orf_to_op["orf1"] == orf_to_op["orf2"] == orf_to_op["orf3"]

    def test_operon_membership_below_threshold_excluded(self, tmp_path):
        p = _write_pred(tmp_path, ["1\t2\t0.3"])
        idx = {1: "orf1", 2: "orf2"}
        orf_to_op, _ = _parse_uniop_pred(str(p), idx)
        assert "orf1" not in orf_to_op
        assert "orf2" not in orf_to_op

    def test_empty_prediction_field_skipped(self, tmp_path):
        p = _write_pred(tmp_path, ["1\t2\t", "1\t3\t0.9"])
        idx = {1: "orf1", 2: "orf2", 3: "orf3"}
        _, pair_probs = _parse_uniop_pred(str(p), idx)
        assert len(pair_probs) == 1

    def test_missing_file_returns_empty(self, tmp_path):
        orf_to_op, pair_probs = _parse_uniop_pred(
            str(tmp_path / "nonexistent.pred"), {1: "orf1"})
        assert orf_to_op == {}
        assert pair_probs == {}

    def test_multiple_pairs(self, tmp_path):
        lines = ["1\t2\t0.9", "1\t3\t0.6", "2\t3\t0.8"]
        p = _write_pred(tmp_path, lines)
        idx = {1: "a", 2: "b", 3: "c"}
        _, pair_probs = _parse_uniop_pred(str(p), idx)
        assert len(pair_probs) == 3


@pytest.mark.parametrize('probability', ['NaN', 'inf', '-0.1', '1.1'])
def test_invalid_probabilities_not_used(tmp_path, probability):
    path=_write_pred(tmp_path,[f'1\t2\t{probability}'])
    assert _parse_uniop_pred(path,{1:'a',2:'b'})==({}, {})


@pytest.mark.parametrize('use_operon',[True,False])
def test_pair_probabilities_read_with_or_without_membership(tmp_path,use_operon):
    from efesto.uniop import run_uniop
    faa=tmp_path/'sample.faa';faa.write_text('>a\nMAA\n>b\nMKK\n')
    work=tmp_path/'_uniop/sample';work.mkdir(parents=True)
    # Existing empty membership also exercises fallback assembly from pred.
    (work/'uniop.operon').write_text('idx_genes,idx_op\n'+('"[1, 2]",7\n' if use_operon else ''))
    (work/'uniop.pred').write_text('Gene A\tGene B\tPrediction\n1\t2\t0.8\n')
    ops,probs=run_uniop([faa],None,tmp_path,'unused',prodigal_faa_dir=tmp_path)
    assert probs=={'sample.faa':{('a','b'):.8}}
    assert ops['sample.faa']['a']==ops['sample.faa']['b']
    if use_operon:assert ops['sample.faa']['a']=='OP0007'
