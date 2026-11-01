import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hic3d_demo as h  # noqa: E402


def test_dataset_is_reproducible_and_symmetric():
    a = h.build_dataset(3)
    b = h.build_dataset(3)
    assert np.array_equal(a["control"], b["control"])
    assert np.array_equal(a["rearranged"], b["rearranged"])
    for key in ("control", "rearranged"):
        assert np.array_equal(a[key], a[key].T)
        assert np.all(np.diag(a[key]) == 0)


def test_rearrangement_is_a_permutation_with_known_breakpoints():
    rng = np.random.default_rng(0)
    bnd = h.make_domains(600, rng)
    r = h.apply_rearrangements(600, bnd, rng, n_inv=3, n_swap=1)
    assert sorted(r.ref_idx) == list(range(600))          # no bin lost or duplicated
    assert len(r.breakpoints) >= 8
    assert np.all((r.breakpoints > 0) & (r.breakpoints < 600))


def test_snap_places_breakpoints_on_tad_boundaries():
    rng = np.random.default_rng(1)
    bnd = h.make_domains(600, rng)
    r = h.apply_rearrangements(600, bnd, rng, snap=True)
    allowed = set(bnd.tolist())
    assert all(int(p) in allowed for p in r.breakpoints)


def test_match_and_prf():
    assert h.match([10, 50], [11, 90], tol=2) == (1, 1, 1)
    res = h.prf([10, 50], [10, 50], tol=0)
    assert res["precision"] == res["recall"] == res["f1"] == 1.0


def test_tad_and_compartment_recovery_beats_chance():
    scores, comp = [], []
    for seed in range(3):
        d = h.build_dataset(seed)
        tads, _ = h.call_tad_boundaries(d["control"])
        scores.append(h.prf(tads, d["boundaries"], tol=2)["f1"])
        comp.append(h.compartment_accuracy(h.call_compartments(d["control"]), d["compartments"]))
    assert np.mean(scores) > 0.8
    assert np.mean(comp) > 0.8


def test_3d_reconstruction_correlates_with_truth_and_beats_shuffle():
    d = h.build_dataset(2)
    rec = h.reconstruct_3d(d["control"])
    q = h.reconstruction_quality(d["xyz"], rec)
    assert q["spearman"] > 0.7
    assert q["procrustes_disparity"] < 0.6 * q["null_disparity_mean"]


def test_breakpoint_caller_finds_most_true_breakpoints_somewhere():
    d = h.build_dataset(4)
    bps, _ = h.call_breakpoints(d["rearranged"])
    assert h.prf(bps, d["breakpoints"], tol=3)["recall"] >= 0.3


def test_enrichment_test_detects_and_rejects():
    n = 600
    boundaries = np.arange(20, 580, 30)
    on_target = boundaries[:8] + 1
    res = h.enrichment_test(on_target, boundaries, n, tol=3, n_perm=2000)
    assert res["observed"] == 8 and res["p_value"] < 0.01
    rng = np.random.default_rng(0)
    pvals = [h.enrichment_test(rng.integers(1, n, 8), boundaries, n, tol=3, n_perm=500, seed=s)["p_value"]
             for s in range(20)]
    assert np.mean(np.array(pvals) < 0.05) < 0.25          # no systematic false positives


def test_cli_enrich(tmp_path, capsys):
    bp, bd = tmp_path / "bp.txt", tmp_path / "bd.txt"
    bp.write_text("21\n51\n81\n")
    bd.write_text("20\n50\n80\n110\n")
    h.main(["enrich", "--breakpoints", str(bp), "--boundaries", str(bd), "--n-bins", "600"])
    assert '"observed": 3' in capsys.readouterr().out
