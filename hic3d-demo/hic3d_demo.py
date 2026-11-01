#!/usr/bin/env python3
"""hic3d_demo: simulate a 3D genome, derive Hi-C contacts, and recover its structure.

The pipeline has a known ground truth at every step, so each analysis can be scored:

1. Simulate a reference genome made of TADs (compact domains) and A/B compartments.
2. Simulate a 3D conformation, turn pairwise distances into Hi-C contact counts.
3. Introduce rearrangements (inversions, segment swaps) in the sampled genome and map the
   reads back to the reference, which produces the off-diagonal Hi-C signatures of
   structural variants.
4. Recover: TAD boundaries (insulation score), A/B compartments (PCA of the O/E matrix),
   rearrangement breakpoints (long-range contact-profile change points) and the 3D
   structure (shortest-path distances + classical MDS, in the spirit of ShRec3D).
5. Test whether independently known breakpoints (e.g. synteny breaks) fall on called TAD
   boundaries more often than chance (permutation test).

Everything is synthetic. The point is to show the methods work and to measure how well.
Usage:
    python hic3d_demo.py run --seed 1 --outdir results
    python hic3d_demo.py benchmark --seeds 20
    python hic3d_demo.py enrich --breakpoints bp.txt --boundaries tads.txt --n-bins 600
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import signal, stats
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist, squareform


# --------------------------------------------------------------------------------------
# 1. Reference genome: domains and compartments
# --------------------------------------------------------------------------------------
def make_domains(n_bins: int, rng: np.random.Generator, size_range=(20, 40)) -> np.ndarray:
    """Return the start bin of every TAD after the first one (the true TAD boundaries)."""
    starts, pos = [], 0
    while True:
        pos += int(rng.integers(size_range[0], size_range[1] + 1))
        if pos >= n_bins - size_range[0]:
            break
        starts.append(pos)
    return np.array(starts, dtype=int)


def domain_ids(n_bins: int, boundaries: np.ndarray) -> np.ndarray:
    ids = np.zeros(n_bins, dtype=int)
    for b in boundaries:
        ids[b:] += 1
    return ids


def make_compartments(n_bins: int, rng: np.random.Generator, size_range=(40, 90)) -> np.ndarray:
    """Alternating A (1) / B (0) blocks along the reference genome."""
    comp = np.zeros(n_bins, dtype=int)
    pos, state = 0, int(rng.integers(0, 2))
    while pos < n_bins:
        size = int(rng.integers(size_range[0], size_range[1] + 1))
        comp[pos:pos + size] = state
        pos += size
        state = 1 - state
    return comp


# --------------------------------------------------------------------------------------
# 2. Rearrangements
# --------------------------------------------------------------------------------------
@dataclass
class Rearranged:
    ref_idx: np.ndarray               # ref_idx[s] = reference bin of sample bin s
    breakpoints: np.ndarray           # true breakpoints, in reference coordinates
    events: list = field(default_factory=list)


def _free(intervals, a, b, gap=5):
    return all(b + gap <= x or a >= y + gap for x, y in intervals)


def apply_rearrangements(n_bins, boundaries, rng, n_inv=3, n_swap=1, snap=False) -> Rearranged:
    """Apply non-overlapping inversions and segment swaps to the reference order.

    snap=True places every breakpoint on a true TAD boundary; snap=False places them at
    random bins. Breakpoint p means: the junction lies between reference bins p-1 and p.
    """
    cand = np.concatenate([[0], boundaries, [n_bins]]) if snap else None
    used, events = [], []

    def pick(lo, hi):
        for _ in range(5000):
            if snap:
                a, b = sorted(rng.choice(cand, 2, replace=False))
            else:
                a = int(rng.integers(1, n_bins - hi))
                b = a + int(rng.integers(lo, hi + 1))
            if lo <= b - a <= hi and _free(used, a, b):
                return int(a), int(b)
        raise RuntimeError("could not place rearrangement")

    ref_idx = np.arange(n_bins)
    bps = []
    for _ in range(n_inv):
        a, b = pick(50, 100)
        used.append((a, b))
        ref_idx[a:b] = ref_idx[a:b][::-1].copy()
        bps += [a, b]
        events.append(("inversion", a, b))
    for _ in range(n_swap):
        a1, b1 = pick(30, 70)
        used.append((a1, b1))
        a2, b2 = pick(30, 70)          # pick() guarantees it does not overlap (a1, b1)
        used.append((a2, b2))
        if a2 < a1:
            a1, b1, a2, b2 = a2, b2, a1, b1
        A, M, B = (ref_idx[a1:b1].copy(), ref_idx[b1:a2].copy(), ref_idx[a2:b2].copy())
        ref_idx[a1:b2] = np.concatenate([B, M, A])
        bps += [a1, b1, a2, b2]
        events.append(("swap", a1, b1, a2, b2))
    bps = np.array(sorted(set(p for p in bps if 0 < p < n_bins)), dtype=int)
    return Rearranged(ref_idx, bps, events)


# --------------------------------------------------------------------------------------
# 3. 3D conformation and Hi-C contacts
# --------------------------------------------------------------------------------------
def simulate_conformation(ref_idx, dom_ref, rng, step=0.5, pull=0.15, link=4.0, confine=0.03):
    """3D coordinates of the sample genome: each run of bins from one TAD is a compact blob."""
    n = len(ref_idx)
    seg = np.zeros(n, dtype=int)
    for s in range(1, n):
        contiguous = abs(ref_idx[s] - ref_idx[s - 1]) == 1 and dom_ref[ref_idx[s]] == dom_ref[ref_idx[s - 1]]
        seg[s] = seg[s - 1] + (0 if contiguous else 1)
    n_seg = seg[-1] + 1
    centres = np.zeros((n_seg, 3))
    for k in range(1, n_seg):
        d = rng.normal(size=3)
        centres[k] = centres[k - 1] * (1 - confine) + link * d / np.linalg.norm(d)
    xyz = np.zeros((n, 3))
    for s in range(n):
        if s == 0 or seg[s] != seg[s - 1]:
            xyz[s] = centres[seg[s]]
        else:
            xyz[s] = xyz[s - 1] + step * rng.normal(size=3) - pull * (xyz[s - 1] - centres[seg[s]])
    return xyz


def simulate_contacts(xyz, comp_sample, rng, depth=60.0, gamma=0.3):
    """Poisson contact counts with expectation ~ 1/distance, modulated by compartment type."""
    d = squareform(pdist(xyz))
    np.fill_diagonal(d, np.inf)
    d0 = np.median(np.linalg.norm(np.diff(xyz, axis=0), axis=1))
    lam = depth * d0 / d
    same = comp_sample[:, None] == comp_sample[None, :]
    lam = lam * np.where(same, 1 + gamma, 1 - gamma)
    counts = np.triu(rng.poisson(lam), 1)
    return counts + counts.T


def build_dataset(seed: int, n_bins: int = 600, n_inv: int = 3, n_swap: int = 1, snap: bool = False):
    """Return a control map (no rearrangements) and a rearranged map, with all ground truth."""
    rng = np.random.default_rng(seed)
    bnd = make_domains(n_bins, rng)
    dom = domain_ids(n_bins, bnd)
    comp = make_compartments(n_bins, rng)

    ident = np.arange(n_bins)
    xyz_c = simulate_conformation(ident, dom, rng)
    control = simulate_contacts(xyz_c, comp, rng)

    rearr = apply_rearrangements(n_bins, bnd, rng, n_inv, n_swap, snap)
    xyz_r = simulate_conformation(rearr.ref_idx, dom, rng)
    m_sample = simulate_contacts(xyz_r, comp[rearr.ref_idx], rng)
    rearranged = np.zeros_like(m_sample)
    rearranged[np.ix_(rearr.ref_idx, rearr.ref_idx)] = m_sample   # map reads to the reference
    return dict(control=control, rearranged=rearranged, boundaries=bnd, compartments=comp,
                xyz=xyz_c, breakpoints=rearr.breakpoints, events=rearr.events, n_bins=n_bins)


# --------------------------------------------------------------------------------------
# 4. Analyses
# --------------------------------------------------------------------------------------
def match(pred, true, tol):
    """One-to-one matching of predicted to true positions within tol; returns (tp, fp, fn)."""
    pred, true = list(sorted(pred)), list(sorted(true))
    used, tp = set(), 0
    for p in pred:
        best, best_d = None, tol + 1
        for k, t in enumerate(true):
            if k not in used and abs(p - t) <= tol and abs(p - t) < best_d:
                best, best_d = k, abs(p - t)
        if best is not None:
            used.add(best); tp += 1
    return tp, len(pred) - tp, len(true) - tp


def prf(pred, true, tol):
    tp, fp, fn = match(pred, true, tol)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return dict(precision=prec, recall=rec, f1=f1, n_called=len(pred), n_true=len(true))


def insulation_score(M, w=10):
    """log2 mean contact across each bin, in a w x w square (low = insulating boundary)."""
    n = M.shape[0]
    s = np.full(n, np.nan)
    for i in range(w, n - w):
        s[i] = M[i - w + 1:i + 1, i + 1:i + w + 1].mean()
    return np.log2(s / np.nanmedian(s))


def call_tad_boundaries(M, w=10, prominence=0.8, min_dist=15):
    """Boundaries are local minima of the insulation score. Returns first bin of each new TAD."""
    s = insulation_score(M, w)
    filled = np.where(np.isnan(s), np.nanmax(s), s)
    peaks, _ = signal.find_peaks(-filled, prominence=prominence, distance=min_dist)
    return peaks + 1, s


def observed_over_expected(M):
    n = M.shape[0]
    oe = np.zeros_like(M, dtype=float)
    for k in range(1, n):
        diag = np.diagonal(M, k)
        e = diag.mean() if diag.mean() > 0 else 1.0
        idx = np.arange(n - k)
        oe[idx, idx + k] = diag / e
        oe[idx + k, idx] = diag / e
    np.fill_diagonal(oe, 1.0)
    return oe


def call_compartments(M):
    """PC1 of the correlation matrix of the O/E map. Sign is arbitrary."""
    corr = np.corrcoef(observed_over_expected(M))
    vals, vecs = np.linalg.eigh(corr)
    return vecs[:, -1]


def compartment_accuracy(pc1, truth):
    pred = (pc1 > 0).astype(int)
    acc = (pred == truth).mean()
    return max(acc, 1 - acc)           # eigenvector sign carries no information


def breakpoint_scores(M, k=4, far=40):
    """Change-point score between bin i and i+1.

    Adjacent windows of k bins should interact with the same distant bins, so their
    contact profiles over columns farther than `far` bins should agree. A rearrangement
    breakpoint joins different partners on the two sides, which breaks that agreement.
    Returns a robust z-score (median / MAD) for each junction i|i+1.
    """
    n = M.shape[0]
    cols = np.arange(n)
    raw = np.full(n - 1, np.nan)
    for i in range(k - 1, n - k):
        left = M[i - k + 1:i + 1].sum(0)
        right = M[i + 1:i + k + 1].sum(0)
        mask = (cols < i - k + 1 - far) | (cols > i + k + far)
        if mask.sum() < 30:
            continue
        l, r = left[mask], right[mask]
        raw[i] = np.mean((l - r) ** 2 / (l + r + 1))
    med = np.nanmedian(raw)
    mad = np.nanmedian(np.abs(raw - med)) * 1.4826
    return (raw - med) / mad


def call_breakpoints(M, z_thresh=8.0, min_dist=8, **kw):
    z = breakpoint_scores(M, **kw)
    filled = np.nan_to_num(z, nan=0.0)
    peaks, _ = signal.find_peaks(filled, height=z_thresh, distance=min_dist)
    return peaks + 1, z


def reconstruct_3d(M, min_count=1):
    """ShRec3D-style: contacts -> distances -> graph shortest paths -> classical MDS (3D)."""
    n = M.shape[0]
    w = np.where(M >= min_count, 1.0 / np.maximum(M, 1), 0.0)
    for i in range(n - 1):                                 # keep the chain connected
        w[i, i + 1] = w[i + 1, i] = 1.0 / max(M[i, i + 1], 1)
    d = shortest_path(csr_matrix(w), directed=False)
    j = np.eye(n) - 1.0 / n
    b = -0.5 * j @ (d ** 2) @ j
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1][:3]
    return vecs[:, order] * np.sqrt(np.maximum(vals[order], 0))


def reconstruction_quality(xyz_true, xyz_rec, rng=None, n_null=200):
    """Spearman correlation of pairwise distances, and Procrustes disparity vs a shuffled null."""
    rho = stats.spearmanr(pdist(xyz_true), pdist(xyz_rec)).statistic
    disparity = procrustes(xyz_true, xyz_rec)[2]
    rng = rng or np.random.default_rng(0)
    null = [procrustes(xyz_true, xyz_rec[rng.permutation(len(xyz_rec))])[2] for _ in range(n_null)]
    return dict(spearman=float(rho), procrustes_disparity=float(disparity),
                null_disparity_mean=float(np.mean(null)))


# --------------------------------------------------------------------------------------
# 5. Enrichment of breakpoints at TAD boundaries
# --------------------------------------------------------------------------------------
def enrichment_test(breakpoints, boundaries, n_bins, tol=3, n_perm=10000, seed=0):
    """Permutation test: are breakpoints closer to TAD boundaries than uniformly random ones?"""
    rng = np.random.default_rng(seed)
    bnd = np.asarray(boundaries)
    bps = np.asarray(breakpoints)

    def hits(pos):
        return int(np.sum(np.min(np.abs(pos[:, None] - bnd[None, :]), axis=1) <= tol)) if len(bnd) else 0

    obs = hits(bps)
    null = np.array([hits(rng.integers(1, n_bins, size=len(bps))) for _ in range(n_perm)])
    p = (1 + np.sum(null >= obs)) / (1 + n_perm)
    return dict(observed=obs, n_breakpoints=len(bps), null_mean=float(null.mean()), p_value=float(p))


# --------------------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------------------
def analyse(seed: int, snap: bool = False, n_bins: int = 600, n_perm: int = 5000):
    data = build_dataset(seed, n_bins=n_bins, snap=snap)
    tads, _ = call_tad_boundaries(data["control"])
    pc1 = call_compartments(data["control"])
    rec = reconstruct_3d(data["control"])
    bps, _ = call_breakpoints(data["rearranged"])
    result = dict(
        seed=seed,
        snap=snap,
        tad=prf(tads, data["boundaries"], tol=2),
        compartment_accuracy=float(compartment_accuracy(pc1, data["compartments"])),
        reconstruction=reconstruction_quality(data["xyz"], rec, np.random.default_rng(seed)),
        breakpoints=prf(bps, data["breakpoints"], tol=3),
        # The enrichment test uses the true breakpoints as a stand-in for independent calls
        # (e.g. synteny breaks from genome comparison), so it does not inherit the Hi-C
        # breakpoint caller's own bias towards domain boundaries (see README).
        enrichment=enrichment_test(data["breakpoints"], tads, n_bins, n_perm=n_perm, seed=seed),
    )
    return result, data, dict(tads=tads, pc1=pc1, rec=rec, bps=bps)


def make_figure(data, calls, result, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = data["n_bins"]
    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.28)

    def show(ax, M, title):
        ax.imshow(np.log1p(M), cmap="Reds", origin="upper")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("bin"); ax.set_ylabel("bin")

    ax = fig.add_subplot(gs[0, 0])
    show(ax, data["control"], "Control contact map (log counts)")
    for b in data["boundaries"]:
        ax.plot([b, b], [0, 14], color="tab:blue", lw=0.8)
    for b in calls["tads"]:
        ax.plot([0, 14], [b, b], color="k", lw=0.8)
    ax.text(0.02, 0.02, "blue: true TAD boundaries (top)\nblack: called (left)", transform=ax.transAxes,
            fontsize=8, va="bottom", bbox=dict(fc="w", alpha=0.8, lw=0))

    ax = fig.add_subplot(gs[0, 1])
    pc1 = calls["pc1"] * (1 if np.corrcoef(calls["pc1"], data["compartments"])[0, 1] > 0 else -1)
    ax.plot(pc1 / np.abs(pc1).max(), color="k", lw=1, label="recovered PC1")
    ax.fill_between(range(n), -1, 1, where=data["compartments"] == 1, color="tab:red", alpha=0.15,
                    label="true A")
    ax.set_title(f"A/B compartments (accuracy {result['compartment_accuracy']:.2f})", fontsize=11)
    ax.set_xlabel("bin"); ax.legend(loc="lower right", fontsize=8)

    ax = fig.add_subplot(gs[0, 2], projection="3d")
    true_al, rec_al, _ = procrustes(data["xyz"], calls["rec"])
    ax.plot(*true_al.T, color="tab:blue", lw=0.8, label="truth")
    ax.plot(*rec_al.T, color="tab:orange", lw=0.8, label="reconstruction (aligned)")
    ax.set_title(f"3D reconstruction (Spearman {result['reconstruction']['spearman']:.2f})", fontsize=11)
    ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
    ax.legend(fontsize=8, loc="upper left")

    ax = fig.add_subplot(gs[1, 0])
    show(ax, data["rearranged"], "Rearranged sample mapped to reference")
    for b in data["breakpoints"]:
        ax.plot([b, b], [0, 14], color="tab:blue", lw=1.2)
    for b in calls["bps"]:
        ax.plot([0, 14], [b, b], color="k", lw=1.2)
    ax.text(0.02, 0.02, "blue: true breakpoints (top)\nblack: called (left)", transform=ax.transAxes,
            fontsize=8, va="bottom", bbox=dict(fc="w", alpha=0.8, lw=0))

    ax = fig.add_subplot(gs[1, 1])
    z = breakpoint_scores(data["rearranged"])
    ax.plot(np.arange(1, n), z, color="k", lw=0.8)
    for b in data["breakpoints"]:
        ax.axvline(b, color="tab:blue", alpha=0.6)
    ax.axhline(8, color="tab:red", ls="--", lw=0.8)
    ax.set_title("Breakpoint score (robust z); blue = truth", fontsize=11)
    ax.set_xlabel("junction (bin)"); ax.set_ylabel("z")

    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    r = result
    lines = [
        "Recovery against ground truth", "",
        f"TAD boundaries   P {r['tad']['precision']:.2f}  R {r['tad']['recall']:.2f}  F1 {r['tad']['f1']:.2f}",
        f"Compartments     accuracy {r['compartment_accuracy']:.2f}",
        f"3D structure     Spearman {r['reconstruction']['spearman']:.2f}",
        f"                 Procrustes disparity {r['reconstruction']['procrustes_disparity']:.2f}",
        f"                 (shuffled null {r['reconstruction']['null_disparity_mean']:.2f})",
        f"Breakpoints(HiC) P {r['breakpoints']['precision']:.2f}  R {r['breakpoints']['recall']:.2f}  F1 {r['breakpoints']['f1']:.2f}",
        f"Breakpoints at TAD boundaries: {r['enrichment']['observed']}/{r['enrichment']['n_breakpoints']} hits, "
        f"p = {r['enrichment']['p_value']:.3f}",
        "", f"seed {r['seed']}, snap-to-TAD = {r['snap']}", "synthetic data",
    ]
    ax.text(0, 1, "\n".join(lines), va="top", family="monospace", fontsize=10)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def cmd_run(args):
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    result, data, calls = analyse(args.seed, snap=args.snap)
    (out / f"metrics_seed{args.seed}.json").write_text(json.dumps(result, indent=2))
    make_figure(data, calls, result, out / f"overview_seed{args.seed}.png")
    print(json.dumps(result, indent=2))


def cmd_benchmark(args):
    rows = {False: [], True: []}
    for snap in (False, True):
        for seed in range(args.seeds):
            rows[snap].append(analyse(seed, snap=snap, n_perm=2000)[0])
    out = {}
    for snap, res in rows.items():
        out["snap_to_tad" if snap else "random_breakpoints"] = dict(
            tad_f1=float(np.mean([r["tad"]["f1"] for r in res])),
            compartment_accuracy=float(np.mean([r["compartment_accuracy"] for r in res])),
            reconstruction_spearman=float(np.mean([r["reconstruction"]["spearman"] for r in res])),
            breakpoint_precision=float(np.mean([r["breakpoints"]["precision"] for r in res])),
            breakpoint_recall=float(np.mean([r["breakpoints"]["recall"] for r in res])),
            enrichment_median_p=float(np.median([r["enrichment"]["p_value"] for r in res])),
            enrichment_frac_p_lt_0_05=float(np.mean([r["enrichment"]["p_value"] < 0.05 for r in res])),
        )
    out["n_seeds"] = args.seeds
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    (Path(args.outdir) / "benchmark.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


def cmd_enrich(args):
    bps = np.loadtxt(args.breakpoints, dtype=int, ndmin=1)
    bnd = np.loadtxt(args.boundaries, dtype=int, ndmin=1)
    print(json.dumps(enrichment_test(bps, bnd, args.n_bins, tol=args.tol), indent=2))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="one simulated dataset: metrics + figure")
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--snap", action="store_true", help="place rearrangement breakpoints on TAD boundaries")
    r.add_argument("--outdir", default="results")
    r.set_defaults(func=cmd_run)
    b = sub.add_parser("benchmark", help="many seeds, both breakpoint placements")
    b.add_argument("--seeds", type=int, default=20)
    b.add_argument("--outdir", default="results")
    b.set_defaults(func=cmd_benchmark)
    e = sub.add_parser("enrich", help="permutation test on your own bin positions")
    e.add_argument("--breakpoints", required=True, help="text file, one bin index per line")
    e.add_argument("--boundaries", required=True, help="text file, one bin index per line")
    e.add_argument("--n-bins", type=int, required=True)
    e.add_argument("--tol", type=int, default=3)
    e.set_defaults(func=cmd_enrich)
    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
