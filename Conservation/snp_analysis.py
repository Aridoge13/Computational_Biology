"""
snp_population_analysis.py
==========================
Population-genomic analysis of a multi-sample SNP VCF, designed for datasets
built from public sequencing data (NCBI SRA / BioProject / BioSample).

What it does
------------
1. QC      Biallelic SNP selection, optional QUAL filter, per-sample and
           per-site missingness, minor-allele-frequency filter, removal of
           sex/mitochondrial/unplaced contigs. A step-by-step QC table is
           written so every filtering decision can be audited.
2. PCA     LD-pruned principal component analysis (Patterson scaling).
3. Fst     Pairwise Hudson Fst between populations, block-jackknife SE.
4. Diversity  Per-population observed/expected heterozygosity and the
           inbreeding coefficient F; per-sample heterozygosity.
5. ROH     Runs of homozygosity called with a two-state HMM (Viterbi), ROH
           length classes and genomic inbreeding coefficient F_ROH.
6. Kinship KING-robust pairwise kinship, with related pairs flagged.
7. ADMIXTURE (optional) Writes a PLINK .bed/.bim/.fam from the LD-pruned
           SNPs, runs the external `admixture` binary over a range of K with
           cross-validation, and plots the ancestry proportions.

Input
-----
--vcf        Multi-sample VCF (plain or bgzip) of SNPs called against ONE
             reference assembly, e.g. produced by this route:
                 SRA run accessions  ->  fasterq-dump / prefetch
                 -> bwa-mem2 / minimap2 alignment -> GATK or bcftools calling
                 -> joint-called, filtered VCF
             Sample names in the VCF must match the metadata file. If the VCF
             carries SRR accessions use them; otherwise rename with
             `bcftools reheader -s`.
--metadata   Table with one row per sample (TSV or CSV, auto-detected).
             An NCBI SRA Run Selector export (SraRunTable.csv) works directly:
                 --sample-col Run --pop-col breed
             (use whichever attribute column holds the population/breed).

Example
-------
    python snp_population_analysis.py \\
        --vcf cohort.snps.vcf.gz --metadata SraRunTable.csv \\
        --sample-col Run --pop-col breed --out results/ \\
        --chroms NC_056054.1,NC_056055.1 \\
        --admixture-k 2 6 --threads 8

Notes and caveats
-----------------
* NCBI reference sequences are accessions (NC_...), not "1", "2", ... Pass
  the autosomes explicitly with --chroms, or exclude X/Y/MT with
  --exclude-chroms. Without either, contigs with fewer than
  --min-contig-variants SNPs (unplaced scaffolds) are dropped.
* ROH are called on variant sites only. Because invariant sites are absent
  from a SNP VCF, heterozygosity per site is much higher than genome-wide
  heterozygosity; the HMM uses per-site rates estimated from the data and
  ROH are positioned in base pairs, so lengths remain meaningful, but
  absolute numbers depend on SNP density and the MAF filter. Report the
  thresholds used (they are saved in run_parameters.json).
* KING-robust kinship assumes limited population structure. With several
  breeds in one VCF, run the script per breed (or subset the metadata)
  for interpretable relatedness.
* Memory: genotypes are held as int8 (variants x samples). 10 million SNPs
  x 100 samples is about 1 GB. Pre-filter very large VCFs with
  `bcftools view -m2 -M2 -v snps` and subset samples if needed.

Requirements: numpy, pandas, matplotlib, scikit-allel, numba (recommended
for ROH speed). ADMIXTURE (https://dalexander.github.io/admixture/) is only
needed for --admixture-k.
"""

import argparse
import json
import logging
import math
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import allel  # noqa: E402

try:
    from numba import njit

    _jit = njit(cache=False)
    HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    def _jit(f):
        return f

    HAVE_NUMBA = False

log = logging.getLogger("snp_population_analysis")

DEFAULT_EXCLUDE = "X,Y,MT,M,chrX,chrY,chrM,chrMT,Z,W,chrZ,chrW"
ROH_CLASS_EDGES_MB = [0, 1, 2, 4, 8, 16, np.inf]
ROH_CLASS_LABELS = ["<1 Mb", "1-2 Mb", "2-4 Mb", "4-8 Mb", "8-16 Mb", ">16 Mb"]


# --------------------------------------------------------------------------
# Input
# --------------------------------------------------------------------------
def read_metadata(path, sample_col, pop_col):
    meta = pd.read_csv(path, sep=None, engine="python", dtype=str)
    for col in (sample_col, pop_col):
        if col not in meta.columns:
            sys.exit(f"Metadata column '{col}' not found. Columns present: {list(meta.columns)}")
    meta = meta[[sample_col, pop_col]].rename(columns={sample_col: "sample", pop_col: "population"})
    meta = meta.dropna().drop_duplicates("sample")
    meta["population"] = meta["population"].str.strip().str.replace(r"\s+", "_", regex=True)
    return meta.set_index("sample")


def vcf_contig_lengths(headers):
    """Contig lengths from ##contig header lines (empty dict if absent)."""
    lengths = {}
    for line in getattr(headers, "headers", []):
        m = re.match(r"##contig=<ID=([^,>]+).*?length=(\d+)", line)
        if m:
            lengths[m.group(1)] = int(m.group(2))
    return lengths


def impute_mean(gn):
    """Return float32 genotype matrix with missing calls replaced by the site mean."""
    x = gn.astype(np.float32)
    miss = gn < 0
    x[miss] = np.nan
    with np.errstate(all="ignore"):
        mean = np.nanmean(x, axis=1)
    mean = np.where(np.isnan(mean), 0.0, mean)
    rows, _ = np.where(miss)
    x[miss] = mean[rows]
    return x


def load_genotypes(args, samples_keep, qc_rows):
    """
    Stream the VCF in chunks, applying site-level filters per chunk so peak
    memory stays close to the size of the *filtered* matrix.
    Returns chrom (object array), pos (int64), gn (int8, -1 = missing),
    samples (list) and contig lengths.
    """
    headers = allel.read_vcf_headers(args.vcf)
    vcf_samples = list(headers.samples)
    contig_len = vcf_contig_lengths(headers)
    use = [s for s in vcf_samples if s in samples_keep]
    if len(use) < 2:
        sys.exit("Fewer than 2 VCF samples match the metadata. Check --sample-col and VCF sample names.")
    missing_in_vcf = sorted(set(samples_keep) - set(vcf_samples))
    if missing_in_vcf:
        log.warning("%d metadata samples not in VCF (ignored): %s", len(missing_in_vcf), missing_in_vcf[:10])
    log.info("Reading %d samples from %s", len(use), args.vcf)

    fields = ["variants/CHROM", "variants/POS", "variants/numalt", "variants/is_snp", "calldata/GT"]
    if args.min_qual is not None:
        fields.append("variants/QUAL")
    _, _, _, chunks = allel.iter_vcf_chunks(
        args.vcf, fields=fields, alt_number=1, samples=use, chunk_length=100_000
    )

    include = set(args.chroms.split(",")) if args.chroms else None
    exclude = set(args.exclude_chroms.split(",")) if args.exclude_chroms else set()

    keep_chrom, keep_pos, keep_gn = [], [], []
    n_total = n_snp = n_chrom = n_qual = n_miss = 0
    for chunk, _, _, _ in chunks:
        chrom = chunk["variants/CHROM"].astype(str)
        n_total += len(chrom)
        m = (chunk["variants/numalt"] == 1) & chunk["variants/is_snp"]
        n_snp += int(m.sum())
        if include is not None:
            c = np.isin(chrom, list(include))
        else:
            c = ~np.isin(chrom, list(exclude))
        m &= c
        n_chrom += int(m.sum())
        if args.min_qual is not None:
            q = chunk["variants/QUAL"]
            m &= np.nan_to_num(q, nan=-1) >= args.min_qual
        n_qual += int(m.sum())
        if not m.any():
            continue
        gt = chunk["calldata/GT"][m]  # (v, s, 2)
        gn = np.where((gt < 0).any(axis=2), -1, gt.sum(axis=2)).astype(np.int8)
        miss_frac = (gn < 0).mean(axis=1)
        ok = miss_frac <= args.max_site_missing_pre
        n_miss += int(ok.sum())
        keep_chrom.append(chrom[m][ok])
        keep_pos.append(chunk["variants/POS"][m][ok].astype(np.int64))
        keep_gn.append(gn[ok])

    if not keep_gn:
        sys.exit("No variants survived the initial filters.")
    qc_rows += [
        ("records in VCF", n_total),
        ("biallelic SNPs", n_snp),
        ("on selected chromosomes", n_chrom),
        ("QUAL filter" if args.min_qual is not None else "QUAL filter (not applied)", n_qual),
        (f"site missingness <= {args.max_site_missing_pre} (pre-filter)", n_miss),
    ]
    return (np.concatenate(keep_chrom), np.concatenate(keep_pos), np.vstack(keep_gn), use, contig_len)


def qc_filter(args, chrom, pos, gn, samples, qc_rows):
    """Sample missingness -> site missingness -> MAF -> contig size filters."""
    s_miss = (gn < 0).mean(axis=0)
    drop = s_miss > args.max_sample_missing
    if drop.any():
        log.warning("Dropping %d samples with missingness > %.2f: %s",
                    int(drop.sum()), args.max_sample_missing, [s for s, d in zip(samples, drop) if d])
    gn = gn[:, ~drop]
    samples = [s for s, d in zip(samples, drop) if not d]
    qc_rows.append((f"samples kept (missingness <= {args.max_sample_missing})", len(samples)))
    if len(samples) < 2:
        sys.exit("Fewer than 2 samples left after QC.")

    keep = (gn < 0).mean(axis=1) <= args.max_site_missing
    chrom, pos, gn = chrom[keep], pos[keep], gn[keep]
    qc_rows.append((f"site missingness <= {args.max_site_missing}", len(pos)))

    called = (gn >= 0).sum(axis=1)
    alt = np.where(gn >= 0, gn, 0).sum(axis=1)
    freq = alt / (2.0 * np.maximum(called, 1))
    maf = np.minimum(freq, 1 - freq)
    keep = maf >= args.maf
    chrom, pos, gn = chrom[keep], pos[keep], gn[keep]
    qc_rows.append((f"MAF >= {args.maf}", len(pos)))

    if args.chroms is None:
        uniq, counts = np.unique(chrom, return_counts=True)
        small = set(uniq[counts < args.min_contig_variants])
        if small:
            keep = ~np.isin(chrom, list(small))
            chrom, pos, gn = chrom[keep], pos[keep], gn[keep]
            log.info("Dropped %d contigs with < %d SNPs", len(small), args.min_contig_variants)
        qc_rows.append((f"contigs with >= {args.min_contig_variants} SNPs", len(pos)))
    if len(pos) == 0:
        sys.exit("No variants left after QC; relax the thresholds.")
    return chrom, pos, gn, samples


def chrom_slices(chrom):
    """Index ranges of contiguous chromosome blocks (VCF must be sorted)."""
    change = np.flatnonzero(chrom[1:] != chrom[:-1]) + 1
    starts = np.concatenate([[0], change])
    ends = np.concatenate([change, [len(chrom)]])
    return [(chrom[s], s, e) for s, e in zip(starts, ends)]


def ld_prune(chrom, gn, window, step, r2):
    """LD pruning per chromosome with scikit-allel; returns boolean mask."""
    mask = np.zeros(len(chrom), dtype=bool)
    for name, s, e in chrom_slices(chrom):
        block = np.rint(impute_mean(gn[s:e])).astype(np.int8)
        mask[s:e] = allel.locate_unlinked(block, size=window, step=step, threshold=r2)
    return mask


# --------------------------------------------------------------------------
# Analyses
# --------------------------------------------------------------------------
def run_pca(gn_pruned, samples, pops, n_comp, out):
    n_comp = min(n_comp, len(samples) - 1)
    x = impute_mean(gn_pruned)
    # drop invariant columns introduced by sample filtering
    x = x[x.std(axis=1) > 0]
    coords, model = allel.pca(x, n_components=n_comp, scaler="patterson")
    cols = [f"PC{i + 1}" for i in range(n_comp)]
    df = pd.DataFrame(coords, columns=cols, index=samples)
    df.insert(0, "population", pops)
    df.index.name = "sample"
    df.to_csv(out / "pca_coordinates.tsv", sep="\t")
    var = pd.Series(model.explained_variance_ratio_, index=cols, name="explained_variance_ratio")
    var.to_csv(out / "pca_explained_variance.tsv", sep="\t")

    colors = pop_colors(pops)
    fig, axes = plt.subplots(1, 2 if n_comp >= 4 else 1, figsize=(12 if n_comp >= 4 else 6, 5), squeeze=False)
    for ax, (i, j) in zip(axes[0], [(0, 1), (2, 3)]):
        for p in sorted(set(pops)):
            sel = np.array(pops) == p
            ax.scatter(df.iloc[sel, i + 1], df.iloc[sel, j + 1], s=35, label=p, color=colors[p], edgecolor="k", lw=0.3)
        ax.set_xlabel(f"PC{i + 1} ({var.iloc[i] * 100:.1f}%)")
        ax.set_ylabel(f"PC{j + 1} ({var.iloc[j] * 100:.1f}%)")
    axes[0][0].legend(fontsize=7, frameon=False, bbox_to_anchor=(1.0, 1.0))
    fig.suptitle(f"PCA on {x.shape[0]:,} LD-pruned SNPs")
    fig.tight_layout()
    fig.savefig(out / "pca.png", dpi=200)
    plt.close(fig)


def allele_counts(gn_sub):
    called = gn_sub >= 0
    an = 2 * called.sum(axis=1)
    alt = np.where(called, gn_sub, 0).sum(axis=1)
    return np.stack([an - alt, alt], axis=1).astype(np.int32), an


def run_fst(gn, samples, pops, min_pop_size, blen, out):
    pops_arr = np.array(pops)
    names = [p for p in sorted(set(pops)) if (pops_arr == p).sum() >= min_pop_size]
    if len(names) < 2:
        log.warning("Fst skipped: need >= 2 populations with >= %d samples", min_pop_size)
        return
    counts = {}
    for p in names:
        counts[p] = allele_counts(gn[:, pops_arr == p])
    mat = pd.DataFrame(np.nan, index=names, columns=names)
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ac1, an1 = counts[a]
            ac2, an2 = counts[b]
            ok = (an1 >= 2) & (an2 >= 2)
            fst, se, _, _ = allel.average_hudson_fst(ac1[ok], ac2[ok], blen=blen)
            mat.loc[a, b] = mat.loc[b, a] = fst
            rows.append((a, b, fst, se, int(ok.sum())))
    for p in names:
        mat.loc[p, p] = 0.0
    mat.to_csv(out / "fst_matrix.tsv", sep="\t")
    pd.DataFrame(rows, columns=["pop1", "pop2", "fst_hudson", "jackknife_se", "n_sites"]).to_csv(
        out / "fst_pairs.tsv", sep="\t", index=False)
    heatmap(mat.values, names, names, "Pairwise Hudson Fst", out / "fst_heatmap.png", annotate=len(names) <= 12)


def run_diversity(gn, samples, pops, min_pop_size, out):
    pops_arr = np.array(pops)
    rows = []
    for p in sorted(set(pops)):
        idx = pops_arr == p
        if idx.sum() < max(2, min_pop_size):
            continue
        sub = gn[:, idx]
        ac, an = allele_counts(sub)
        called = (sub >= 0).sum(axis=1)
        het = (sub == 1).sum(axis=1)
        with np.errstate(all="ignore"):
            ho = np.where(called > 0, het / called, np.nan)
            he = allel.mean_pairwise_difference(ac)
            poly = (ac[:, 0] > 0) & (ac[:, 1] > 0)
        ho_m, he_m = np.nanmean(ho), np.nanmean(he)
        rows.append({
            "population": p, "n_samples": int(idx.sum()), "n_sites": len(ho),
            "polymorphic_fraction": poly.mean(), "Ho": ho_m, "He": he_m, "F": 1 - ho_m / he_m,
        })
    pd.DataFrame(rows).to_csv(out / "diversity_by_population.tsv", sep="\t", index=False)


def sample_summary(gn, samples, pops):
    called = (gn >= 0).sum(axis=0)
    het = (gn == 1).sum(axis=0)
    return pd.DataFrame({
        "sample": samples, "population": pops, "n_called": called,
        "missingness": 1 - called / gn.shape[0], "het_rate": het / np.maximum(called, 1),
    }).set_index("sample")


# ---- ROH ------------------------------------------------------------------
@_jit
def _viterbi(het, pos, lp_het0, lp_hom0, lp_het1, lp_hom1, rate_in, rate_out):
    """
    Two-state HMM. State 0 = non-ROH, state 1 = ROH.
    Emission: heterozygous / homozygous call at each variant site.
    Transition probabilities depend on the physical distance d (bp) between
    consecutive sites: P(switch) = 1 - exp(-rate * d).
    Returns the most likely state path (int8 array).
    """
    n = het.shape[0]
    back = np.zeros((n, 2), dtype=np.int8)
    pi1 = rate_in / (rate_in + rate_out)
    v0 = math.log(1.0 - pi1) + (lp_het0 if het[0] == 1 else lp_hom0)
    v1 = math.log(pi1) + (lp_het1 if het[0] == 1 else lp_hom1)
    for i in range(1, n):
        d = pos[i] - pos[i - 1]
        if d < 1:
            d = 1
        stay0 = -rate_in * d
        to1 = math.log(-math.expm1(-rate_in * d))
        stay1 = -rate_out * d
        to0 = math.log(-math.expm1(-rate_out * d))
        a0 = v0 + stay0
        b0 = v1 + to0
        a1 = v1 + stay1
        b1 = v0 + to1
        if a0 >= b0:
            n0 = a0
            back[i, 0] = 0
        else:
            n0 = b0
            back[i, 0] = 1
        if a1 >= b1:
            n1 = a1
            back[i, 1] = 1
        else:
            n1 = b1
            back[i, 1] = 0
        v0 = n0 + (lp_het0 if het[i] == 1 else lp_hom0)
        v1 = n1 + (lp_het1 if het[i] == 1 else lp_hom1)
    path = np.zeros(n, dtype=np.int8)
    state = 0 if v0 >= v1 else 1
    path[n - 1] = state
    for i in range(n - 1, 0, -1):
        state = back[i, state]
        path[i - 1] = state
    return path


def call_roh(args, chrom, pos, gn, samples, pops, het_rates, out):
    p_roh = args.roh_p_het_roh
    segments = []
    blocks = chrom_slices(chrom)
    for si, s in enumerate(samples):
        p_non = args.roh_p_het_non if args.roh_p_het_non else float(np.clip(het_rates[si], 0.05, 0.5))
        lp = (math.log(p_non), math.log(1 - p_non), math.log(p_roh), math.log(1 - p_roh))
        for name, a, b in blocks:
            g = gn[a:b, si]
            ok = g >= 0
            if ok.sum() < args.roh_min_snps:
                continue
            het = (g[ok] == 1).astype(np.int8)
            p = pos[a:b][ok]
            path = _viterbi(het, p, lp[0], lp[1], lp[2], lp[3], args.roh_rate_in, args.roh_rate_out)
            d = np.diff(np.concatenate([[0], path, [0]]))
            starts, ends = np.flatnonzero(d == 1), np.flatnonzero(d == -1) - 1
            for i0, i1 in zip(starts, ends):
                length = int(p[i1] - p[i0] + 1)
                n_snps = int(i1 - i0 + 1)
                if length >= args.roh_min_length and n_snps >= args.roh_min_snps:
                    segments.append((s, pops[si], name, int(p[i0]), int(p[i1]), length, n_snps,
                                     int(het[i0:i1 + 1].sum())))
    roh = pd.DataFrame(segments, columns=["sample", "population", "chrom", "start", "end",
                                         "length_bp", "n_snps", "n_het_in_roh"])
    roh.to_csv(out / "roh_segments.tsv", sep="\t", index=False)
    return roh


def summarise_roh(roh, samples, pops, genome_len, out):
    cls = pd.cut(roh["length_bp"] / 1e6, bins=ROH_CLASS_EDGES_MB, labels=ROH_CLASS_LABELS, right=False) \
        if len(roh) else pd.Series(dtype="category")
    summ = pd.DataFrame({"sample": samples, "population": pops}).set_index("sample")
    g = roh.groupby("sample")
    summ["n_roh"] = g.size().reindex(summ.index).fillna(0).astype(int)
    summ["total_roh_bp"] = g["length_bp"].sum().reindex(summ.index).fillna(0).astype(np.int64)
    summ["F_ROH"] = summ["total_roh_bp"] / genome_len
    for lab in ROH_CLASS_LABELS:
        sub = roh[cls == lab] if len(roh) else roh
        summ[f"roh_bp_{lab}"] = sub.groupby("sample")["length_bp"].sum().reindex(summ.index).fillna(0).astype(np.int64)
    summ.to_csv(out / "roh_summary.tsv", sep="\t")

    colors = pop_colors(pops)
    order = summ.sort_values(["population", "F_ROH"]).index
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [2.5, 1]})
    bottom = np.zeros(len(order))
    cmap = plt.get_cmap("viridis", len(ROH_CLASS_LABELS))
    for k, lab in enumerate(ROH_CLASS_LABELS):
        v = summ.loc[order, f"roh_bp_{lab}"].values / 1e6
        axes[0].bar(range(len(order)), v, bottom=bottom, color=cmap(k), label=lab, width=0.85)
        bottom += v
    axes[0].set_xticks(range(len(order)))
    axes[0].set_xticklabels([f"{i}\n{summ.loc[i, 'population']}" for i in order], rotation=90, fontsize=6)
    axes[0].set_ylabel("Total ROH length (Mb)")
    axes[0].legend(title="ROH length", fontsize=7, frameon=False)
    popnames = sorted(set(pops))
    data = [summ.loc[summ.population == p, "F_ROH"].values for p in popnames]
    bp = axes[1].boxplot(data, patch_artist=True, showfliers=False)
    for patch, p in zip(bp["boxes"], popnames):
        patch.set_facecolor(colors[p])
    axes[1].set_xticks(range(1, len(popnames) + 1))
    axes[1].set_xticklabels(popnames, rotation=60, ha="right", fontsize=7)
    axes[1].set_ylabel("F_ROH")
    fig.tight_layout()
    fig.savefig(out / "roh.png", dpi=200)
    plt.close(fig)


# ---- Kinship ----------------------------------------------------------------
def king_robust(gn_pruned):
    """
    KING-robust kinship (Manichaikul et al. 2010, between-family estimator):
        phi_ij = (N_AaAa - 2 * N_AAaa) / (N_Aa(i) + N_Aa(j))
    computed for all pairs by matrix products on genotype indicators.
    """
    g = gn_pruned
    het = (g == 1).astype(np.float32)
    h0 = (g == 0).astype(np.float32)
    h2 = (g == 2).astype(np.float32)
    n_aa_aa = het.T @ het
    n_opposite = h0.T @ h2 + h2.T @ h0
    n_het = het.sum(axis=0)
    denom = n_het[:, None] + n_het[None, :]
    with np.errstate(all="ignore"):
        phi = (n_aa_aa - 2 * n_opposite) / denom
    np.fill_diagonal(phi, 0.5)
    return phi


def run_kinship(gn_pruned, samples, pops, threshold, out):
    phi = king_robust(gn_pruned)
    df = pd.DataFrame(phi, index=samples, columns=samples)
    df.to_csv(out / "kinship_king_matrix.tsv", sep="\t")
    rows = []
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            if phi[i, j] >= threshold:
                degree = ("duplicate/MZ" if phi[i, j] > 0.354 else "1st degree" if phi[i, j] > 0.177
                          else "2nd degree" if phi[i, j] > 0.0884 else "3rd degree")
                rows.append((samples[i], pops[i], samples[j], pops[j], phi[i, j], degree))
    pd.DataFrame(rows, columns=["sample1", "pop1", "sample2", "pop2", "kinship", "inferred_degree"]).to_csv(
        out / "kinship_related_pairs.tsv", sep="\t", index=False)
    pop_order = np.argsort(pops, kind="stable")
    heatmap(phi[np.ix_(pop_order, pop_order)], [samples[i] for i in pop_order], [samples[i] for i in pop_order],
            "KING-robust kinship", out / "kinship_heatmap.png", annotate=False, vmin=-0.1, vmax=0.5, cmap="magma")


# ---- ADMIXTURE ---------------------------------------------------------------
def write_plink_bed(prefix, chrom, pos, gn, samples):
    """
    Write a PLINK 1 binary fileset (SNP-major). A1 = ALT allele, A2 = REF.
    Chromosomes are written as integers because ADMIXTURE needs numeric codes.
    Genotype codes: 00 = hom A1, 01 = missing, 10 = het, 11 = hom A2.
    """
    n_var, n_s = gn.shape
    code = np.full(gn.shape, 1, dtype=np.uint8)
    code[gn == 2] = 0
    code[gn == 1] = 2
    code[gn == 0] = 3
    pad = (-n_s) % 4
    if pad:
        code = np.hstack([code, np.zeros((n_var, pad), dtype=np.uint8)])
    c = code.reshape(n_var, -1, 4)
    packed = (c[:, :, 0] | (c[:, :, 1] << 2) | (c[:, :, 2] << 4) | (c[:, :, 3] << 6)).astype(np.uint8)
    with open(f"{prefix}.bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        fh.write(packed.tobytes())
    chrom_codes = {c_: i + 1 for i, c_ in enumerate(dict.fromkeys(chrom))}
    with open(f"{prefix}.bim", "w") as fh:
        for c_, p_ in zip(chrom, pos):
            fh.write(f"{chrom_codes[c_]}\t{c_}:{p_}\t0\t{p_}\tA\tG\n")
    with open(f"{prefix}.fam", "w") as fh:
        for s in samples:
            fh.write(f"{s}\t{s}\t0\t0\t0\t-9\n")


def run_admixture(args, chrom_p, pos_p, gn_p, samples, pops, out):
    exe = args.admixture_exe
    adm = out / "admixture"
    adm.mkdir(exist_ok=True)
    prefix = adm / "pruned"
    write_plink_bed(str(prefix), chrom_p, pos_p, gn_p, samples)
    kmin, kmax = args.admixture_k
    cv = {}
    for k in range(kmin, kmax + 1):
        cmd = [exe, f"--cv={args.cv_folds}", f"-j{args.threads}", f"-s{args.seed}", "pruned.bed", str(k)]
        log.info("Running: %s", " ".join(cmd))
        try:
            res = subprocess.run(cmd, cwd=adm, capture_output=True, text=True, check=True)
        except FileNotFoundError:
            log.error("ADMIXTURE executable '%s' not found; skipping ancestry analysis. "
                      "PLINK files were still written to %s", exe, adm)
            return
        except subprocess.CalledProcessError as e:
            log.error("ADMIXTURE failed for K=%d:\n%s", k, e.stdout[-2000:] + e.stderr[-2000:])
            return
        (adm / f"log{k}.out").write_text(res.stdout)
        m = re.search(r"CV error \(K=\d+\): ([0-9.eE+-]+)", res.stdout)
        if m:
            cv[k] = float(m.group(1))
    if not cv:
        return
    cv_df = pd.Series(cv, name="cv_error").rename_axis("K")
    cv_df.to_csv(out / "admixture_cv_error.tsv", sep="\t")
    best = int(cv_df.idxmin())
    log.info("ADMIXTURE lowest CV error at K=%d", best)

    order = np.lexsort((np.arange(len(samples)), np.array(pops)))
    ks = sorted(cv)
    fig, axes = plt.subplots(len(ks) + 1, 1, figsize=(max(8, len(samples) * 0.18), 2 + 1.6 * len(ks)),
                             gridspec_kw={"height_ratios": [1.2] + [1] * len(ks)})
    axes[0].plot(ks, [cv[k] for k in ks], "o-")
    axes[0].set_ylabel("CV error")
    axes[0].set_xticks(ks)
    for ax, k in zip(axes[1:], ks):
        q = np.loadtxt(adm / f"pruned.{k}.Q")[order]
        cmap = plt.get_cmap("tab10")
        bottom = np.zeros(len(order))
        for c in range(k):
            ax.bar(range(len(order)), q[:, c], bottom=bottom, color=cmap(c % 10), width=1.0, lw=0)
            bottom += q[:, c]
        ax.set_ylabel(f"K={k}")
        ax.set_xlim(-0.5, len(order) - 0.5)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        pd.DataFrame(np.loadtxt(adm / f"pruned.{k}.Q"), index=samples,
                     columns=[f"anc{c + 1}" for c in range(k)]).assign(population=pops).to_csv(
            out / f"admixture_K{k}_ancestry.tsv", sep="\t")
    pop_sorted = np.array(pops)[order]
    axes[-1].set_xticks(range(len(order)))
    axes[-1].set_xticklabels([f"{samples[i]} | {pops[i]}" for i in order], rotation=90, fontsize=5)
    fig.tight_layout()
    fig.savefig(out / "admixture.png", dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------
# Plot helpers
# --------------------------------------------------------------------------
def pop_colors(pops):
    names = sorted(set(pops))
    cmap = plt.get_cmap("tab20" if len(names) > 10 else "tab10")
    return {p: cmap(i % cmap.N) for i, p in enumerate(names)}


def heatmap(mat, xlabels, ylabels, title, path, annotate=True, vmin=None, vmax=None, cmap="viridis"):
    fig, ax = plt.subplots(figsize=(max(5, 0.5 * len(xlabels) + 3), max(4, 0.5 * len(ylabels) + 2)))
    im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(xlabels)))
    ax.set_yticks(range(len(ylabels)))
    fs = 8 if len(xlabels) <= 30 else 4
    ax.set_xticklabels(xlabels, rotation=90, fontsize=fs)
    ax.set_yticklabels(ylabels, fontsize=fs)
    if annotate:
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if not np.isnan(mat[i, j]):
                    ax.text(j, i, f"{mat[i, j]:.3f}", ha="center", va="center", fontsize=6, color="w")
    fig.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True)
    p.add_argument("--metadata", required=True, help="TSV/CSV with sample and population columns")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--sample-col", default="sample")
    p.add_argument("--pop-col", default="population")
    p.add_argument("--analyses", default="pca,diversity,fst,roh,kinship",
                   help="Comma-separated subset of: pca,diversity,fst,roh,kinship")
    # selection and QC
    p.add_argument("--chroms", default=None, help="Comma-separated contigs to analyse (e.g. the autosomes)")
    p.add_argument("--exclude-chroms", default=DEFAULT_EXCLUDE, help="Ignored when --chroms is given")
    p.add_argument("--min-contig-variants", type=int, default=1000,
                   help="Without --chroms, drop contigs with fewer SNPs (unplaced scaffolds)")
    p.add_argument("--min-qual", type=float, default=None)
    p.add_argument("--max-site-missing-pre", type=float, default=0.5,
                   help="Loose per-site missingness applied while streaming the VCF")
    p.add_argument("--max-site-missing", type=float, default=0.1)
    p.add_argument("--max-sample-missing", type=float, default=0.2)
    p.add_argument("--maf", type=float, default=0.01)
    # pruning / PCA / Fst
    p.add_argument("--ld-window", type=int, default=50, help="LD window in SNPs")
    p.add_argument("--ld-step", type=int, default=10, help="LD step in SNPs")
    p.add_argument("--ld-r2", type=float, default=0.2)
    p.add_argument("--pca-components", type=int, default=10)
    p.add_argument("--min-pop-size", type=int, default=3)
    p.add_argument("--fst-block", type=int, default=1000, help="Jackknife block length in SNPs")
    # ROH
    p.add_argument("--roh-min-length", type=int, default=500_000, help="Minimum ROH length (bp)")
    p.add_argument("--roh-min-snps", type=int, default=50)
    p.add_argument("--roh-p-het-roh", type=float, default=0.005,
                   help="P(het call | ROH): residual genotyping error")
    p.add_argument("--roh-p-het-non", type=float, default=None,
                   help="P(het call | non-ROH); default is each sample's own het rate (clipped 0.05-0.5)")
    p.add_argument("--roh-rate-in", type=float, default=2e-8, help="non-ROH -> ROH switch rate per bp")
    p.add_argument("--roh-rate-out", type=float, default=5e-7, help="ROH -> non-ROH switch rate per bp")
    # kinship
    p.add_argument("--kinship-threshold", type=float, default=0.0884)
    # admixture
    p.add_argument("--admixture-k", type=int, nargs=2, metavar=("KMIN", "KMAX"), default=None)
    p.add_argument("--admixture-exe", default="admixture")
    p.add_argument("--cv-folds", type=int, default=5)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    analyses = {a.strip() for a in args.analyses.split(",") if a.strip()}
    unknown = analyses - {"pca", "diversity", "fst", "roh", "kinship"}
    if unknown:
        sys.exit(f"Unknown analyses: {unknown}")
    if "roh" in analyses and not HAVE_NUMBA:
        log.warning("numba not installed: ROH calling will be very slow on large VCFs (pip install numba)")
    t0 = time.time()

    meta = read_metadata(args.metadata, args.sample_col, args.pop_col)
    qc_rows = []
    chrom, pos, gn, samples, contig_len = load_genotypes(args, set(meta.index), qc_rows)
    chrom, pos, gn, samples = qc_filter(args, chrom, pos, gn, samples, qc_rows)
    pops = [meta.loc[s, "population"] for s in samples]
    pd.DataFrame(qc_rows, columns=["step", "remaining"]).to_csv(out / "qc_summary.tsv", sep="\t", index=False)
    log.info("After QC: %d SNPs, %d samples, %d populations (%s)", len(pos), len(samples),
             len(set(pops)), dict(pd.Series(pops).value_counts()))

    ss = sample_summary(gn, samples, pops)
    ss.to_csv(out / "sample_summary.tsv", sep="\t")

    pruned = None
    if analyses & {"pca", "kinship"} or args.admixture_k:
        pruned = ld_prune(chrom, gn, args.ld_window, args.ld_step, args.ld_r2)
        qc_rows.append((f"LD-pruned (r2<{args.ld_r2}, {args.ld_window}/{args.ld_step})", int(pruned.sum())))
        pd.DataFrame(qc_rows, columns=["step", "remaining"]).to_csv(out / "qc_summary.tsv", sep="\t", index=False)
        log.info("LD pruning kept %d / %d SNPs", pruned.sum(), len(pruned))

    if "pca" in analyses:
        log.info("PCA")
        run_pca(gn[pruned], samples, pops, args.pca_components, out)
    if "diversity" in analyses:
        log.info("Diversity")
        run_diversity(gn, samples, pops, args.min_pop_size, out)
    if "fst" in analyses:
        log.info("Fst")
        run_fst(gn, samples, pops, args.min_pop_size, args.fst_block, out)
    if "roh" in analyses:
        log.info("ROH")
        roh = call_roh(args, chrom, pos, gn, samples, pops, ss["het_rate"].values, out)
        used = dict.fromkeys(chrom)
        if contig_len and all(c in contig_len for c in used):
            genome_len = sum(contig_len[c] for c in used)
        else:
            genome_len = sum(int(pos[e - 1] - pos[s] + 1) for _, s, e in chrom_slices(chrom))
            log.warning("Contig lengths not in VCF header; F_ROH uses the SNP-covered span (%d bp)", genome_len)
        summarise_roh(roh, samples, pops, genome_len, out)
        log.info("%d ROH segments in %d samples", len(roh), roh["sample"].nunique() if len(roh) else 0)
    if "kinship" in analyses:
        log.info("Kinship")
        run_kinship(gn[pruned], samples, pops, args.kinship_threshold, out)
    if args.admixture_k:
        log.info("ADMIXTURE")
        run_admixture(args, chrom[pruned], pos[pruned], gn[pruned], samples, pops, out)

    params = vars(args).copy()
    params.update({"n_snps_after_qc": int(len(pos)), "n_samples": len(samples),
                   "versions": {"numpy": np.__version__, "pandas": pd.__version__, "scikit-allel": allel.__version__},
                   "runtime_seconds": round(time.time() - t0, 1)})
    (out / "run_parameters.json").write_text(json.dumps(params, indent=2, default=str))
    log.info("Done in %.1fs. Results in %s", time.time() - t0, out)


if __name__ == "__main__":
    main()