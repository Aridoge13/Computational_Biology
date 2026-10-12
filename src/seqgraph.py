#!/usr/bin/env python3
"""
seqgraph.py - sequence variation graphs for viral populations, in one file.

Takes a multiple sequence alignment of a viral population (for example the
haplotypes of a quasispecies, optionally grouped by sample), builds a
sequence variation graph, and runs graph algorithms on it:

  * graph construction from an alignment, with linear-chain compaction
  * topological ordering, heaviest path (consensus) and k-best path enumeration
  * counting of source-to-sink paths (the haplotype space implied by the graph)
  * bubble (variant site) detection on the DAG, with allele tallies per sample
  * per-site entropy, hotspot windows, codon impact for SNPs in annotated genes
  * per-sample diversity and rank-correlation with numeric metadata
    (sample-level diversity and allele-level frequencies, BH-adjusted)
  * GFA 1.0 export and import, with a round-trip check against the input
  * a self-contained interactive HTML viewer (layered graph, entropy track,
    gene track, per-sample highlighting)
  * a simulator that plants known hotspots and a known metadata association,
    so the whole pipeline can be benchmarked against ground truth

Only the Python standard library is used (Python 3.9 or newer).

Usage
-----
  python seqgraph.py demo --out demo/demo           # simulate, build, benchmark
  python seqgraph.py build aln.fasta --out run/x    # build and analyse
        [--ref NAME] [--genes genes.tsv] [--meta meta.tsv] [--k 5] [--window 30]
  python seqgraph.py check run/x.gfa aln.fasta      # verify GFA spells the input
  python seqgraph.py simulate --out sim/x           # write a simulated dataset
  python seqgraph.py test                           # run the unit tests

Input conventions
-----------------
  FASTA alignment, equal-length records. Optional header tokens:
      >hap1 sample=S01 count=37
  `count` (or `size`) is the abundance of the haplotype, `sample` groups
  haplotypes into samples. The first record is the reference unless --ref is set.
  Gaps are '-' or '.'; any other non-ACGT symbol is treated as missing.
  genes.tsv: gene<TAB>start<TAB>end (1-based, reference coordinates, + strand).
  meta.tsv: header row with 'sample' first; numeric columns are analysed.

Variants in the first or last alignment column cannot form a bubble (no flanking
sequence on both sides) and are not reported. Limits: pure Python, intended for alignments up to a few thousand columns and a
few thousand haplotypes. Codon impact assumes one change per codon.
"""
from __future__ import annotations

import argparse
import collections
import heapq
import html
import io
import json
import math
import os
import random
import sys
import unittest

GAPS = "-."
BASES = "ACGT"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {
    a + b + c: _AA[16 * i + 4 * j + k]
    for i, a in enumerate("TCAG")
    for j, b in enumerate("TCAG")
    for k, c in enumerate("TCAG")
}


# --------------------------------------------------------------------------
# Input
# --------------------------------------------------------------------------
class Rec:
    __slots__ = ("name", "seq", "count", "sample")

    def __init__(self, name, seq, count=1.0, sample=None):
        self.name, self.seq, self.count = name, seq, count
        self.sample = sample if sample is not None else name


def parse_fasta(fh):
    """Yield (header, sequence) pairs from an open FASTA handle."""
    head, buf = None, []
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if head is not None:
                yield head, "".join(buf)
            head, buf = line[1:].strip(), []
        else:
            buf.append(line.upper().replace("U", "T"))
    if head is not None:
        yield head, "".join(buf)


def read_msa(fh):
    """Read an alignment. Returns (records, number_of_masked_symbols)."""
    recs, masked, seen = [], 0, collections.Counter()
    for head, seq in parse_fasta(fh):
        toks = head.split()
        kv = dict(t.split("=", 1) for t in toks[1:] if "=" in t)
        name = toks[0] if toks else f"seq{len(recs) + 1}"
        seen[name] += 1
        if seen[name] > 1:
            name = f"{name}.{seen[name]}"
        clean = []
        for ch in seq:
            if ch in GAPS or ch in BASES:
                clean.append(ch)
            else:
                clean.append("-")
                masked += 1
        count = float(kv.get("count", kv.get("size", 1)))
        recs.append(Rec(name, "".join(clean), count, kv.get("sample", name)))
    if not recs:
        raise ValueError("no sequences found")
    if len({len(r.seq) for r in recs}) != 1:
        raise ValueError("records are not aligned: lengths differ")
    return recs, masked


def load_genes(path):
    genes = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or line.startswith("#") or not p[1].strip().isdigit():
                continue
            genes.append((p[0], int(p[1]), int(p[2])))
    return genes


def load_meta(path):
    with open(path) as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip() and not l.startswith("#")]
    head = rows[0]
    data = {r[0]: dict(zip(head[1:], r[1:])) for r in rows[1:]}
    numeric = []
    for c in head[1:]:
        try:
            [float(d[c]) for d in data.values()]
            numeric.append(c)
        except (ValueError, KeyError):
            pass
    return data, numeric, [c for c in head[1:] if c not in numeric]


# --------------------------------------------------------------------------
# Graph
# --------------------------------------------------------------------------
class Graph:
    """Directed graph of sequence segments with weighted edges and named paths."""

    def __init__(self):
        self.seq = {}                                   # node -> sequence
        self.cols = {}                                  # node -> (first, last) alignment column
        self.out = collections.defaultdict(dict)        # u -> {v: weight}
        self.inn = collections.defaultdict(dict)        # v -> {u: weight}
        self.paths = {}                                 # path name -> [nodes]
        self.pmeta = {}                                 # path name -> (count, sample)

    def add_edge(self, u, v, w):
        self.out[u][v] = self.out[u].get(v, 0.0) + w
        self.inn[v][u] = self.inn[v].get(u, 0.0) + w

    def n_edges(self):
        return sum(len(d) for d in self.out.values())

    def coverage(self):
        cov = collections.defaultdict(float)
        for name, p in self.paths.items():
            c = self.pmeta[name][0]
            for n in p:
                cov[n] += c
        return cov

    def topo_order(self):
        """Kahn's algorithm. Ties are broken by alignment column, then name."""
        indeg = {n: len(self.inn[n]) for n in self.seq}
        key = lambda n: (self.cols.get(n, (0, 0))[0], n)
        heap = [key(n) + (n,) for n, d in indeg.items() if d == 0]
        heapq.heapify(heap)
        order = []
        while heap:
            n = heapq.heappop(heap)[-1]
            order.append(n)
            for v in self.out[n]:
                indeg[v] -= 1
                if indeg[v] == 0:
                    heapq.heappush(heap, key(v) + (v,))
        if len(order) != len(self.seq):
            raise ValueError("graph contains a cycle; a DAG is required")
        return order

    def spell(self, name):
        return "".join(self.seq[n] for n in self.paths[name])


def build_graph(recs):
    """Alignment -> graph with one node per (column, base); gaps become edge bypasses."""
    g, ids = Graph(), {}
    names = []
    for r in recs:
        path = []
        for c, ch in enumerate(r.seq):
            if ch in GAPS:
                continue
            key = (c, ch)
            if key not in ids:
                nid = f"n{len(ids)}"
                ids[key] = nid
                g.seq[nid] = ch
                g.cols[nid] = (c, c)
            path.append(ids[key])
        if not path:
            continue
        names.append(r.name)
        g.paths[r.name] = path
        g.pmeta[r.name] = (r.count, r.sample)
        for a, b in zip(path, path[1:]):
            g.add_edge(a, b, r.count)
    return g, names


def compact(g):
    """Merge unbranched chains (u -> v where u has one successor and v one predecessor)."""
    starts = collections.Counter(p[0] for p in g.paths.values())
    ends = collections.Counter(p[-1] for p in g.paths.values())
    nodes = sorted(g.seq, key=lambda n: (g.cols[n], g.seq[n]))
    nxt = {}
    for u in nodes:
        if len(g.out[u]) == 1 and u not in ends:
            (v,) = g.out[u]
            if len(g.inn[v]) == 1 and v not in starts:
                nxt[u] = v
    has_pred = set(nxt.values())
    rep, ng, k = {}, Graph(), 0
    for u in nodes:
        if u in has_pred:
            continue
        chain = [u]
        while chain[-1] in nxt:
            chain.append(nxt[chain[-1]])
        k += 1
        nid = f"s{k}"
        for x in chain:
            rep[x] = nid
        ng.seq[nid] = "".join(g.seq[x] for x in chain)
        ng.cols[nid] = (g.cols[chain[0]][0], g.cols[chain[-1]][1])
    for u, d in g.out.items():
        for v, w in d.items():
            if rep[u] != rep[v]:
                ng.add_edge(rep[u], rep[v], w)
    for name, p in g.paths.items():
        q = []
        for x in p:
            if not q or q[-1] != rep[x]:
                q.append(rep[x])
        ng.paths[name] = q
        ng.pmeta[name] = g.pmeta[name]
    return ng


# --------------------------------------------------------------------------
# Graph algorithms
# --------------------------------------------------------------------------
def heaviest_path(g, order, cov):
    """Path maximising the summed node coverage (a graph-aware consensus)."""
    best, back = {}, {}
    for n in order:
        cand = [(best[p], p) for p in g.inn[n]]
        b, p = max(cand, key=lambda x: x[0]) if cand else (0.0, None)
        best[n], back[n] = b + cov[n], p
    end = max(order, key=lambda n: best[n])
    path = []
    while end is not None:
        path.append(end)
        end = back[end]
    return path[::-1]


def count_paths(g, order):
    """Number of source-to-sink paths: the haplotype space the graph can express."""
    cnt = {}
    for n in order:
        cnt[n] = sum(cnt[p] for p in g.inn[n]) if g.inn[n] else 1
    return sum(cnt[n] for n in order if not g.out[n])


def kbest_paths(g, order, cov, k):
    """The k highest-scoring source-to-sink paths (dynamic programming on the DAG)."""
    best = {}
    for n in order:
        cands = []
        if not g.inn[n]:
            cands.append((cov[n], None, 0))
        for p in g.inn[n]:
            for r, (sc, _, _) in enumerate(best[p]):
                cands.append((sc + cov[n], p, r))
        cands.sort(key=lambda x: -x[0])
        best[n] = cands[:k]
    ends = [(sc, n, r) for n in order if not g.out[n] for r, (sc, _, _) in enumerate(best[n])]
    ends.sort(key=lambda x: -x[0])
    result = []
    for sc, n, r in ends[:k]:
        path = []
        while n is not None:
            path.append(n)
            _, pn, pr = best[n][r]
            n, r = pn, pr
        result.append((sc, path[::-1]))
    return result


def find_bubbles(g, order):
    """Find (source, sink) pairs of bubbles in a DAG.

    From each branching node, nodes are visited in topological order. The
    bubble closes at the node where all pending branches have converged,
    provided no visited node has a predecessor outside the bubble and no
    branch ends early (a tip). Worst case O(n * bubble size).
    """
    idx = {n: i for i, n in enumerate(order)}
    found = []
    for s in order:
        if len(g.out[s]) < 2:
            continue
        inside = {s}
        heap = [(idx[v], v) for v in g.out[s]]
        heapq.heapify(heap)
        queued = {v for _, v in heap}
        sink = None
        while heap:
            _, v = heapq.heappop(heap)
            if any(p not in inside for p in g.inn[v]):
                break
            inside.add(v)
            if not heap:
                sink = v
                break
            if not g.out[v]:
                break
            for w in g.out[v]:
                if w not in queued:
                    queued.add(w)
                    heapq.heappush(heap, (idx[w], w))
        if sink is not None:
            found.append((s, sink))
    return found


# --------------------------------------------------------------------------
# Statistics helpers (no numpy/scipy)
# --------------------------------------------------------------------------
def entropy(counts):
    tot = sum(counts)
    if tot <= 0:
        return 0.0
    return -sum((c / tot) * math.log2(c / tot) for c in counts if c > 0)


def ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for t in range(i, j + 1):
            r[order[t]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def _betacf(a, b, x):
    tiny, qab, qap, qam = 1e-300, a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        d = d if abs(d) > tiny else tiny
        c = 1 + aa / c
        c = c if abs(c) > tiny else tiny
        d = 1 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        d = d if abs(d) > tiny else tiny
        c = 1 + aa / c
        c = c if abs(c) > tiny else tiny
        d = 1 / d
        de = d * c
        h *= de
        if abs(de - 1) < 3e-14:
            break
    return h


def betai(a, b, x):
    """Regularised incomplete beta function I_x(a, b)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                  + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * _betacf(a, b, x) / a
    return 1 - bt * _betacf(b, a, 1 - x) / b


def corr_pvalue(r, n):
    """Two-sided p-value for a correlation r from n pairs (t approximation)."""
    if n < 4 or r != r:
        return float("nan")
    if abs(r) >= 1:
        return 0.0
    df = n - 2
    t2 = r * r * df / (1 - r * r)
    return betai(df / 2, 0.5, df / (df + t2))


def bh_adjust(pvals):
    """Benjamini-Hochberg q-values; NaN p-values stay NaN."""
    idx = [i for i, p in enumerate(pvals) if p == p]
    idx.sort(key=lambda i: pvals[i])
    q, m, prev = [float("nan")] * len(pvals), len(idx), 1.0
    for rank in range(m, 0, -1):
        i = idx[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        q[i] = prev
    return q


# --------------------------------------------------------------------------
# Bubbles: variant sites with allele tallies
# --------------------------------------------------------------------------
class Bubble:
    def __init__(self, source, sink):
        self.source, self.sink = source, sink
        self.alleles = collections.Counter()                         # allele string -> weight
        self.samples = collections.defaultdict(collections.Counter)  # sample -> allele counts
        self.nodes = {}                                              # allele string -> node tuple
        self.ref_allele = None
        self.pos = 0
        self.end = 0
        self.insertion = False
        self.kind = ""
        self.entropy = 0.0
        self.maf = 0.0
        self.gene = ""
        self.impact = []

    @property
    def coverage(self):
        return sum(self.alleles.values())


def col_to_ref(ref_aln):
    """For each alignment column: (reference position, is_insertion_relative_to_reference)."""
    m, pos = [], 0
    for ch in ref_aln:
        if ch in GAPS:
            m.append((pos, True))
        else:
            pos += 1
            m.append((pos, False))
    return m


def collect_bubbles(g, order, ref_name, c2r):
    pos_in = {name: {n: i for i, n in enumerate(p)} for name, p in g.paths.items()}
    bubbles = []
    for s, t in find_bubbles(g, order):
        b = Bubble(s, t)
        for name, p in g.paths.items():
            ps = pos_in[name]
            if s in ps and t in ps and ps[s] < ps[t]:
                inner = tuple(p[ps[s] + 1:ps[t]])
                a = "".join(g.seq[n] for n in inner)
                cnt, sample = g.pmeta[name]
                b.alleles[a] += cnt
                b.samples[sample][a] += cnt
                b.nodes.setdefault(a, inner)
                if name == ref_name:
                    b.ref_allele = a
        live = [a for a, c in b.alleles.items() if c > 0]
        if len(live) < 2:
            continue
        lens = {len(a) for a in live}
        if lens == {1}:
            b.kind = "SNP"
        elif len(lens) == 1:
            b.kind = "MNP"
        else:
            b.kind = "indel"
        b.entropy = entropy(list(b.alleles.values()))
        b.maf = 1 - max(b.alleles.values()) / b.coverage
        anchor = min(g.cols[n][0] for inner in b.nodes.values() for n in inner)
        b.pos, b.insertion = c2r[anchor]
        b.end = c2r[max(g.cols[n][1] for inner in b.nodes.values() for n in inner)][0]
        bubbles.append(b)
    bubbles.sort(key=lambda b: (b.pos, b.source))
    return bubbles


def annotate(bubbles, refseq, genes):
    """Attach gene names and codon-level impact for SNPs and indels in genes."""
    for b in bubbles:
        for name, gs, ge in genes:
            if gs <= b.pos <= ge:
                b.gene = name
                break
        if not b.gene or b.ref_allele is None:
            continue
        gs, ge = next((s, e) for n, s, e in genes if n == b.gene)
        ref = b.ref_allele
        for alt in b.alleles:
            if alt == ref or b.alleles[alt] <= 0:
                continue
            if b.kind == "SNP" and len(ref) == 1 and not b.insertion:
                off = (b.pos - gs) % 3
                cs = b.pos - off
                if cs + 2 > ge:
                    continue
                codon = refseq[cs - 1:cs + 2]
                new = codon[:off] + alt + codon[off + 1:]
                a1, a2 = CODON.get(codon, "X"), CODON.get(new, "X")
                if a1 == a2:
                    label = "synonymous"
                elif a2 == "*":
                    label = "stop_gained"
                else:
                    label = "missense"
                b.impact.append(f"{ref}>{alt}:{label}({a1}>{a2})")
            elif len(alt) != len(ref):
                d = len(alt) - len(ref)
                b.impact.append(f"{ref or '-'}>{alt or '-'}:" + ("frameshift" if d % 3 else "inframe_indel"))


def hotspot_windows(bubbles, window, top):
    """Non-overlapping windows with the highest summed site entropy."""
    pts = sorted((b.pos, b.entropy) for b in bubbles)
    scored = [(sum(e for q, e in pts if p <= q < p + window), p) for p, _ in pts]
    scored.sort(reverse=True)
    chosen = []
    for sc, p in scored:
        if all(abs(p - q) >= window for _, q in chosen):
            chosen.append((sc, p))
        if len(chosen) == top:
            break
    return [(p, p + window - 1, sc) for sc, p in chosen]


# --------------------------------------------------------------------------
# Sample-level statistics and metadata association
# --------------------------------------------------------------------------
def sample_stats(g, bubbles):
    totals = collections.defaultdict(float)
    haps = collections.defaultdict(collections.Counter)
    for name, p in g.paths.items():
        c, s = g.pmeta[name]
        if c > 0:
            totals[s] += c
            haps[s][tuple(p)] += c
    rows = {}
    for s, tot in totals.items():
        pi, poly = 0.0, 0
        for b in bubbles:
            cnt = b.samples.get(s)
            if not cnt:
                continue
            n = sum(cnt.values())
            f2 = sum((c / n) ** 2 for c in cnt.values())
            pi += 1 - f2
            if len([c for c in cnt.values() if c > 0]) > 1:
                poly += 1
        rows[s] = {
            "sample": s,
            "weight": tot,
            "haplotypes": len(haps[s]),
            "hap_entropy": entropy(list(haps[s].values())),
            "pi": pi,
            "poly_sites": poly,
        }
    return rows


def associations(bubbles, stats, meta, numeric, min_n=5, freq_range=(0.02, 0.98)):
    """Spearman correlation of sample metrics and allele frequencies with numeric metadata."""
    out_sample, out_allele = [], []
    samples = [s for s in stats if s in meta]
    for col in numeric:
        for metric in ("hap_entropy", "pi", "haplotypes", "poly_sites"):
            xs = [float(meta[s][col]) for s in samples]
            ys = [float(stats[s][metric]) for s in samples]
            if len(samples) >= min_n:
                r = spearman(xs, ys)
                out_sample.append({"metadata": col, "metric": metric, "n": len(samples),
                                   "rho": r, "p": corr_pvalue(r, len(samples))})
        tests = []
        for bi, b in enumerate(bubbles):
            tot = sum(b.alleles.values())
            for a, c in b.alleles.items():
                if not (freq_range[0] <= c / tot <= freq_range[1]):
                    continue
                xs, ys = [], []
                for s in samples:
                    cnt = b.samples.get(s)
                    if cnt and sum(cnt.values()) > 0:
                        xs.append(float(meta[s][col]))
                        ys.append(cnt.get(a, 0.0) / sum(cnt.values()))
                if len(xs) < min_n:
                    continue
                r = spearman(xs, ys)
                tests.append({"metadata": col, "bubble": bi, "pos": b.pos, "allele": a or "-",
                              "ref": b.ref_allele == a, "gene": b.gene, "n": len(xs),
                              "rho": r, "p": corr_pvalue(r, len(xs))})
        qs = bh_adjust([t["p"] for t in tests])
        for t, q in zip(tests, qs):
            t["q"] = q
        out_allele.extend(tests)
    out_allele.sort(key=lambda t: (t["q"] if t["q"] == t["q"] else 2, t["p"] if t["p"] == t["p"] else 2))
    return out_sample, out_allele


# --------------------------------------------------------------------------
# GFA 1.0 export and import
# --------------------------------------------------------------------------
def _fmt(x):
    return str(int(x)) if float(x).is_integer() else f"{x:.6g}"


def write_gfa(g, fh):
    cov = g.coverage()
    fh.write("H\tVN:Z:1.0\n")
    for n in g.topo_order():
        fh.write(f"S\t{n}\t{g.seq[n]}\tDP:f:{_fmt(cov[n])}\tCL:Z:{g.cols[n][0]}-{g.cols[n][1]}\n")
    for u in g.topo_order():
        for v, w in sorted(g.out[u].items()):
            fh.write(f"L\t{u}\t+\t{v}\t+\t0M\tRC:f:{_fmt(w)}\n")
    for name, p in g.paths.items():
        c, s = g.pmeta[name]
        fh.write(f"P\t{name}\t{','.join(n + '+' for n in p)}\t*\tWT:f:{_fmt(c)}\tSM:Z:{s}\n")


def read_gfa(fh):
    g = Graph()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        tags = {t[:2]: t[5:] for t in f if len(t) > 5 and t[2] == ":" and t[4] == ":"}
        if f[0] == "S":
            g.seq[f[1]] = f[2]
            cl = tags.get("CL", "0-0").split("-")
            g.cols[f[1]] = (int(cl[0]), int(cl[1]))
        elif f[0] == "L":
            g.add_edge(f[1], f[3], float(tags.get("RC", 1)))
        elif f[0] == "P":
            g.paths[f[1]] = [x[:-1] for x in f[2].split(",")]
            g.pmeta[f[1]] = (float(tags.get("WT", 1)), tags.get("SM", f[1]))
    for n in g.seq:
        g.inn[n], g.out[n]
    return g


def check_roundtrip(gfa_path, aln_path):
    """Verify that every path in the GFA spells exactly the ungapped input sequence."""
    with open(aln_path) as fh:
        recs, _ = read_msa(fh)
    with open(gfa_path) as fh:
        g = read_gfa(fh)
    bad = []
    for r in recs:
        want = "".join(ch for ch in r.seq if ch in BASES)
        if r.name not in g.paths:
            bad.append((r.name, "missing path"))
        elif g.spell(r.name) != want:
            bad.append((r.name, "sequence differs"))
    g.topo_order()
    return len(recs), bad


# --------------------------------------------------------------------------
# Simulator with planted ground truth
# --------------------------------------------------------------------------
def _poisson(rng, lam):
    L, k, p = math.exp(-lam), 0, 1.0
    while True:
        p *= rng.random()
        if p <= L:
            return k
        k += 1


def simulate(length=600, n_samples=16, haps=12, seed=1, n_hot=6):
    """Simulate quasispecies haplotypes in alignment form.

    Truth that the pipeline should recover:
      * `n_hot` hotspot sites that mutate about 25 times faster than the rest
      * one SNP whose alternative-allele frequency tracks the metadata value `temp`
      * higher `temp` gives a larger quasispecies (more mutations per haplotype)
    """
    rng = random.Random(seed)
    root = [rng.choice(BASES) for _ in range(length)]
    ins_after = sorted(rng.sample(range(30, length - 30, 7), 3))
    ins_seq = {i: [rng.choice(BASES) for _ in range(3)] for i in ins_after}
    layout, col_of = [], {}
    for i in range(length):
        col_of[i] = len(layout)
        layout.append(("r", i))
        if i in ins_after:
            layout += [("i", i, k) for k in range(3)]
    template = [root[x[1]] if x[0] == "r" else "-" for x in layout]
    sites = list(range(12, length - 12))
    hot = sorted(rng.sample(sites, n_hot))
    weights = {i: (25.0 if i in hot else 1.0) for i in sites}
    assoc = rng.choice([i for i in range(length // 2 + 20, length - 30) if i not in hot])
    assoc_alt = rng.choice([b for b in BASES if b != root[assoc]])
    temps = [0.05 + 0.9 * i / (n_samples - 1) for i in range(n_samples)]
    rng.shuffle(temps)
    pick = list(weights)
    wlist = [weights[i] for i in pick]

    def mutate(row, n_mut, prob_assoc):
        for _ in range(n_mut):
            u = rng.random()
            i = rng.choices(pick, wlist)[0]
            c = col_of[i]
            if u < 0.88:
                if row[c] != "-":
                    row[c] = rng.choice([b for b in BASES if b != row[c]])
            elif u < 0.95:
                for k in range(rng.randint(1, 3)):
                    if i + k < length - 12 and row[col_of[i + k]] != "-":
                        row[col_of[i + k]] = "-"
            else:
                site = rng.choice(ins_after)
                for k in range(3):
                    row[col_of[site] + 1 + k] = ins_seq[site][k]
        if rng.random() < prob_assoc and row[col_of[assoc]] != "-":
            row[col_of[assoc]] = assoc_alt

    lines = [f">ref sample=reference count=0\n{''.join(template)}\n"]
    meta = ["sample\ttemp\tph"]
    for si, temp in enumerate(temps):
        sname = f"S{si + 1:02d}"
        meta.append(f"{sname}\t{10 + 20 * temp:.2f}\t{rng.uniform(5.5, 8.0):.2f}")
        for h in range(haps):
            row = list(template)
            mutate(row, 1 + _poisson(rng, 0.3 + 3.5 * temp), temp)
            cnt = 1 + int(rng.expovariate(1 / 20))
            lines.append(f">{sname}_h{h + 1:02d} sample={sname} count={cnt}\n{''.join(row)}\n")
    third = (length // 2 // 3) * 3
    genes = [("geneA", 1, third), ("geneB", third + 1, (length // 3) * 3)]
    truth = {"hotspots": [i + 1 for i in hot], "assoc_pos": assoc + 1, "assoc_alt": assoc_alt,
             "assoc_metadata": "temp", "seed": seed}
    return "".join(lines), "\n".join(meta) + "\n", genes, truth


def write_simulation(prefix, **kw):
    fasta, meta, genes, truth = simulate(**kw)
    os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)
    with open(prefix + ".fasta", "w") as f:
        f.write(fasta)
    with open(prefix + ".meta.tsv", "w") as f:
        f.write(meta)
    with open(prefix + ".genes.tsv", "w") as f:
        f.write("gene\tstart\tend\n" + "".join(f"{n}\t{s}\t{e}\n" for n, s, e in genes))
    with open(prefix + ".truth.json", "w") as f:
        json.dump(truth, f, indent=2)
    return truth


# --------------------------------------------------------------------------
# HTML viewer
# --------------------------------------------------------------------------
CSS = """
:root{--bg:#fff;--fg:#1c2430;--mut:#667085;--line:#d5dae1;--node:#2f5fa7;--bb:#c2571a;--edge:#8793a5;
--syn:#2a8f84;--mis:#d98324;--ind:#7b4fb0;--oth:#8793a5;--card:#f6f8fa}
@media (prefers-color-scheme:dark){:root{--bg:#14181f;--fg:#e4e8ee;--mut:#9aa5b5;--line:#2c3441;
--node:#6c9be0;--bb:#f08a4b;--edge:#6b778a;--card:#1b212b}}
body{font:14px/1.45 system-ui,sans-serif;background:var(--bg);color:var(--fg);margin:0;padding:16px 20px}
h1{font-size:20px;margin:0 0 4px}h2{font-size:15px;margin:22px 0 6px}
.mut{color:var(--mut)}.scroll{overflow-x:auto;border:1px solid var(--line);border-radius:6px;background:var(--card)}
table{border-collapse:collapse;font-size:13px}td,th{padding:3px 10px;border-bottom:1px solid var(--line);text-align:left}
th{color:var(--mut);font-weight:600}select{font:inherit;padding:2px 6px}
svg text{fill:var(--fg);font-family:ui-monospace,monospace}
.node rect{fill:var(--node)}.node.bb rect{fill:var(--bb)}.node text{fill:#fff;font-size:10px;text-anchor:middle}
.edge{fill:none;stroke:var(--edge)}.cards span{display:inline-block;margin:0 18px 6px 0}
"""

JS = """
const D=JSON.parse(document.getElementById('d').textContent);
function setSample(s){const nf=D.node[s],ef=D.edge[s];
for(const id in nf){const e=document.getElementById('n_'+id);if(e)e.style.opacity=0.18+0.82*nf[id];}
for(const id in ef){const e=document.getElementById('e_'+id);if(e)e.style.strokeOpacity=0.1+0.9*ef[id];}}
document.getElementById('sel').addEventListener('change',e=>setSample(e.target.value));
setSample('all');
"""


def _layout(g, order, backbone, cov):
    layer = {}
    for n in order:
        layer[n] = max((layer[p] + 1 for p in g.inn[n]), default=0)
    by = collections.defaultdict(list)
    for n in order:
        by[layer[n]].append(n)
    y = {}
    for l in sorted(by):
        def key(n):
            ps = [y[p] for p in g.inn[n] if p in y]
            return (sum(ps) / len(ps) if ps else 0.0, -cov[n])
        nodes = sorted(by[l], key=key)
        for i, n in enumerate(nodes):
            y[n] = i - (len(nodes) - 1) / 2
    return layer, y, max(len(v) for v in by.values())


def render_html(title, g, order, cov, backbone, bubbles, genes, ref_len, windows, gene_rows,
                srows, assoc_sample, assoc_allele, summary):
    esc = html.escape
    LS, NW, NH, RS = 64, 50, 20, 30
    layer, y, widest = _layout(g, order, backbone, cov)
    width = (max(layer.values()) + 1) * LS + 40
    height = widest * RS + 50
    cy = height / 2
    maxw = max((w for d in g.out.values() for w in d.values()), default=1) or 1
    bset = set(backbone)
    svg = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']
    for u in order:
        for v, w in g.out[u].items():
            x1, y1 = layer[u] * LS + 20 + NW, cy + y[u] * RS
            x2, y2 = layer[v] * LS + 20, cy + y[v] * RS
            xm = (x1 + x2) / 2
            sw = 0.8 + 3.2 * math.sqrt(w / maxw)
            svg.append(f'<path id="e_{u}_{v}" class="edge" stroke-width="{sw:.2f}" '
                       f'd="M{x1:.1f},{y1:.1f} C{xm:.1f},{y1:.1f} {xm:.1f},{y2:.1f} {x2:.1f},{y2:.1f}">'
                       f'<title>{u} to {v}: weight {_fmt(w)}</title></path>')
    for n in order:
        s = g.seq[n]
        label = s if len(s) <= 5 else f"{len(s)} nt"
        x, yy = layer[n] * LS + 20, cy + y[n] * RS - NH / 2
        cls = "node bb" if n in bset else "node"
        svg.append(f'<g id="n_{n}" class="{cls}"><title>{n}: {esc(s[:60])}{"..." if len(s) > 60 else ""} '
                   f'(coverage {_fmt(cov[n])}, columns {g.cols[n][0] + 1}-{g.cols[n][1] + 1})</title>'
                   f'<rect x="{x}" y="{yy:.1f}" width="{NW}" height="{NH}" rx="4"/>'
                   f'<text x="{x + NW / 2}" y="{yy + 14:.1f}">{esc(label)}</text></g>')
    svg.append("</svg>")

    TW, TH = 1100, 90
    sx = lambda p: 10 + (TW - 20) * p / max(ref_len, 1)
    track = [f'<svg width="{TW}" height="{TH + 50}" xmlns="http://www.w3.org/2000/svg">']
    for a, b, sc in windows:
        track.append(f'<rect x="{sx(a):.1f}" y="2" width="{max(sx(b) - sx(a), 2):.1f}" height="{TH}" '
                     f'fill="var(--bb)" opacity="0.14"><title>hotspot window {a}-{b}, summed entropy {sc:.2f}</title></rect>')
    top = max((b.entropy for b in bubbles), default=1) or 1
    for b in bubbles:
        imp = " ".join(b.impact)
        colour = ("var(--mis)" if "missense" in imp or "stop" in imp or "frameshift" in imp else
                  "var(--syn)" if "synonymous" in imp else "var(--ind)" if b.kind == "indel" else "var(--oth)")
        h = (TH - 8) * b.entropy / top
        track.append(f'<rect x="{sx(b.pos) - 1.2:.1f}" y="{TH - h:.1f}" width="2.4" height="{h:.1f}" fill="{colour}">'
                     f'<title>pos {b.pos} {b.kind} {esc(b.gene)} entropy {b.entropy:.2f} bits; '
                     f'alleles {esc(json.dumps({k or "-": _fmt(v) for k, v in b.alleles.items()}))} {esc(imp)}</title></rect>')
    track.append(f'<line x1="10" x2="{TW - 10}" y1="{TH}" y2="{TH}" stroke="var(--line)"/>')
    for name, gs, ge in genes:
        track.append(f'<rect x="{sx(gs):.1f}" y="{TH + 8}" width="{sx(ge) - sx(gs):.1f}" height="16" rx="3" '
                     f'fill="var(--line)"/><text x="{sx(gs) + 6:.1f}" y="{TH + 20}" font-size="11">{esc(name)}</text>')
    for p in range(0, ref_len + 1, 100):
        track.append(f'<text x="{sx(p):.1f}" y="{TH + 44}" font-size="10" text-anchor="middle">{p}</text>')
    track.append("</svg>")

    total = collections.defaultdict(float)
    nw = collections.defaultdict(lambda: collections.defaultdict(float))
    ew = collections.defaultdict(lambda: collections.defaultdict(float))
    for name, p in g.paths.items():
        c, s = g.pmeta[name]
        if c <= 0:
            continue
        for tag in ("all", s):
            total[tag] += c
            for n in p:
                nw[tag][n] += c
            for a, b in zip(p, p[1:]):
                ew[tag][f"{a}_{b}"] += c
    data = {"node": {t: {n: nw[t][n] / total[t] for n in g.seq} for t in total},
            "edge": {t: {k: v / total[t] for k, v in ew[t].items()} for t in total}}
    opts = "".join(f'<option value="{esc(s)}">{esc(s)}</option>' for s in ["all"] + sorted(k for k in total if k != "all"))

    def table(head, rows):
        return ('<div class="scroll"><table><tr>' + "".join(f"<th>{esc(h)}</th>" for h in head) + "</tr>"
                + "".join("<tr>" + "".join(f"<td>{esc(str(c))}</td>" for c in r) + "</tr>" for r in rows)
                + "</table></div>")

    f4 = lambda x: "nan" if x != x else f"{x:.3g}"
    cards = "".join(f"<span><b>{esc(str(v))}</b> <span class='mut'>{esc(k)}</span></span>" for k, v in summary.items())
    parts = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
             f"<title>{esc(title)}</title><style>{CSS}</style></head><body>",
             f"<h1>{esc(title)}</h1><div class='cards'>{cards}</div>",
             "<h2>Variation graph</h2><p class='mut'>Orange nodes are the heaviest path. Edge width is the number of "
             "haplotype copies using the edge. Choose a sample to fade the graph to what that sample carries. "
             "Sample: <select id='sel'>" + opts + "</select></p>",
             "<div class='scroll'>" + "".join(svg) + "</div>",
             "<h2>Site entropy along the reference</h2><p class='mut'>Bars are variant sites (bubbles). Colour: "
             "orange missense, stop or frameshift; teal synonymous; purple indel; grey other. Shaded bands are hotspot windows.</p>",
             "<div class='scroll'>" + "".join(track) + "</div>",
             "<h2>Hotspot windows</h2>",
             table(["start", "end", "summed entropy (bits)"], [(a, b, f"{s:.2f}") for a, b, s in windows]),
             "<h2>Genes</h2>",
             table(["gene", "variant sites", "missense/stop/frameshift", "synonymous", "indel", "entropy per bp"],
                   [(r["gene"], r["sites"], r["nonsyn"], r["syn"], r["indel"], f"{r['density']:.4f}") for r in gene_rows]),
             "<h2>Samples</h2>",
             table(["sample", "weight", "haplotypes", "haplotype entropy", "diversity (pi)", "polymorphic sites"],
                   [(r["sample"], _fmt(r["weight"]), r["haplotypes"], f"{r['hap_entropy']:.2f}", f"{r['pi']:.2f}", r["poly_sites"])
                    for r in sorted(srows, key=lambda r: r["sample"])])]
    if assoc_sample:
        parts += ["<h2>Metadata and sample-level diversity</h2>",
                  table(["metadata", "metric", "n", "Spearman rho", "p"],
                        [(a["metadata"], a["metric"], a["n"], f4(a["rho"]), f4(a["p"])) for a in assoc_sample])]
    if assoc_allele:
        parts += ["<h2>Metadata and allele frequency (top 15 by q-value)</h2>",
                  table(["metadata", "position", "gene", "allele", "n", "Spearman rho", "p", "q (BH)"],
                        [(a["metadata"], a["pos"], a["gene"], a["allele"], a["n"], f4(a["rho"]), f4(a["p"]), f4(a["q"]))
                         for a in assoc_allele[:15]])]
    data_json = json.dumps(data).replace("</", "<\\/")
    parts.append(f"<script id='d' type='application/json'>{data_json}</script><script>{JS}</script></body></html>")
    return "".join(parts)


# --------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------
def human(n):
    return f"{n:,}" if n < 10 ** 9 else f"{n:.3e}"


def run_build(aln, out, ref=None, genes_path=None, meta_path=None, k=5, window=30, quiet=False):
    say = (lambda *a: None) if quiet else print
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(aln) as fh:
        recs, masked = read_msa(fh)
    ref_rec = next((r for r in recs if r.name == ref), None) if ref else recs[0]
    if ref_rec is None:
        raise SystemExit(f"reference '{ref}' not found in {aln}")
    ncol = len(recs[0].seq)
    say(f"alignment: {len(recs)} sequences x {ncol} columns" + (f" ({masked} non-ACGT symbols treated as gaps)" if masked else ""))

    g0, _ = build_graph(recs)
    g = compact(g0)
    order = g.topo_order()
    cov = g.coverage()
    say(f"graph: {len(g0.seq):,} column nodes -> {len(g.seq):,} nodes, {g.n_edges():,} edges after chain compaction")

    backbone = heaviest_path(g, order, cov)
    n_paths = count_paths(g, order)
    distinct = len({tuple(p) for n, p in g.paths.items() if g.pmeta[n][0] > 0})
    say(f"haplotype space: {human(n_paths)} source-to-sink paths vs {distinct} distinct observed haplotypes")
    observed = {tuple(p) for p in g.paths.values()}
    kb = kbest_paths(g, order, cov, k)
    unobserved = sum(1 for _, p in kb if tuple(p) not in observed)

    c2r = col_to_ref(ref_rec.seq)
    bubbles = collect_bubbles(g, order, ref_rec.name, c2r)
    refseq = "".join(ch for ch in ref_rec.seq if ch in BASES)
    genes = load_genes(genes_path) if genes_path else []
    annotate(bubbles, refseq, genes)
    kinds = collections.Counter(b.kind for b in bubbles)
    say(f"variant sites (bubbles): {len(bubbles)}  " + "  ".join(f"{a}={n}" for a, n in sorted(kinds.items())))
    say(f"top-{k} implied haplotypes: {unobserved} not observed in the data")

    windows = hotspot_windows(bubbles, window, 5)
    gene_rows = []
    for name, gs, ge in genes:
        bs = [b for b in bubbles if b.gene == name]
        imp = " ".join(" ".join(b.impact) for b in bs)
        gene_rows.append({"gene": name, "sites": len(bs),
                          "nonsyn": sum(imp.count(x) for x in ("missense", "stop_gained", "frameshift")),
                          "syn": imp.count("synonymous"), "indel": sum(1 for b in bs if b.kind == "indel"),
                          "density": sum(b.entropy for b in bs) / (ge - gs + 1)})
    stats = sample_stats(g, bubbles)

    a_sample, a_allele, numeric = [], [], []
    if meta_path:
        meta, numeric, skipped = load_meta(meta_path)
        a_sample, a_allele = associations(bubbles, stats, meta, numeric)
        if skipped:
            say(f"metadata: non-numeric columns ignored: {', '.join(skipped)}")

    with open(out + ".gfa", "w") as fh:
        write_gfa(g, fh)
    with open(out + ".bubbles.tsv", "w") as fh:
        fh.write("pos\tend\tkind\tgene\tentropy_bits\tmaf\tcoverage\talleles\timpact\n")
        for b in bubbles:
            al = ";".join(f"{a or '-'}:{_fmt(c)}" for a, c in b.alleles.most_common())
            fh.write(f"{b.pos}\t{b.end}\t{b.kind}\t{b.gene}\t{b.entropy:.4f}\t{b.maf:.4f}\t{_fmt(b.coverage)}\t{al}\t{' '.join(b.impact)}\n")
    with open(out + ".samples.tsv", "w") as fh:
        fh.write("sample\tweight\thaplotypes\thap_entropy\tpi\tpoly_sites\n")
        for r in sorted(stats.values(), key=lambda r: r["sample"]):
            fh.write(f"{r['sample']}\t{_fmt(r['weight'])}\t{r['haplotypes']}\t{r['hap_entropy']:.4f}\t{r['pi']:.4f}\t{r['poly_sites']}\n")
    if meta_path:
        with open(out + ".assoc.tsv", "w") as fh:
            fh.write("level\tmetadata\tpos\tgene\tallele\tn\trho\tp\tq\n")
            for a in a_sample:
                fh.write(f"sample\t{a['metadata']}\t{a['metric']}\t\t\t{a['n']}\t{a['rho']:.4f}\t{a['p']:.3g}\t\n")
            for a in a_allele:
                fh.write(f"allele\t{a['metadata']}\t{a['pos']}\t{a['gene']}\t{a['allele']}\t{a['n']}\t{a['rho']:.4f}\t{a['p']:.3g}\t{a['q']:.3g}\n")
    summary = {"sequences": len(recs), "columns": ncol, "nodes": len(g.seq), "edges": g.n_edges(),
               "variant sites": len(bubbles), "source-to-sink paths": human(n_paths),
               "observed haplotypes": distinct, "samples": len(stats)}
    with open(out + ".html", "w") as fh:
        fh.write(render_html(os.path.basename(out) + ": variation graph", g, order, cov, backbone, bubbles, genes,
                             len(refseq), windows, gene_rows, list(stats.values()), a_sample, a_allele, summary))
    say(f"wrote {out}.gfa .bubbles.tsv .samples.tsv" + (" .assoc.tsv" if meta_path else "") + " .html")
    return {"bubbles": bubbles, "stats": stats, "assoc_sample": a_sample, "assoc_allele": a_allele,
            "windows": windows, "n_paths": n_paths, "distinct": distinct, "graph": g}


def run_demo(out, seed=1):
    truth = write_simulation(out, seed=seed)
    print(f"simulated: hotspots at {truth['hotspots']}, planted association at "
          f"position {truth['assoc_pos']} ({truth['assoc_alt']} allele vs temp)")
    res = run_build(out + ".fasta", out, genes_path=out + ".genes.tsv", meta_path=out + ".meta.tsv")
    ranked = sorted(res["bubbles"], key=lambda b: -b.entropy)
    topn = ranked[:len(truth["hotspots"]) * 2]
    hit = [h for h in truth["hotspots"] if any(b.pos - 1 <= h <= b.end + 1 for b in topn)]
    print(f"\nbenchmark against simulated truth")
    print(f"  hotspots inside one of the top {len(truth['hotspots']) * 2} entropy variant sites: {len(hit)} of {len(truth['hotspots'])}")
    top = [a for a in res["assoc_allele"] if a["metadata"] == "temp"][:3]
    rank = next((i + 1 for i, a in enumerate(a for a in res["assoc_allele"] if a["metadata"] == "temp")
                 if a["pos"] == truth["assoc_pos"] and a["allele"] == truth["assoc_alt"]), None)
    if rank:
        a = [a for a in res["assoc_allele"] if a["metadata"] == "temp"][rank - 1]
        print(f"  planted association ranked #{rank} of {sum(1 for a in res['assoc_allele'] if a['metadata'] == 'temp')} "
              f"allele tests (rho={a['rho']:.2f}, q={a['q']:.2g})")
    else:
        print("  planted association not found among tested alleles")
    pi = next((a for a in res["assoc_sample"] if a["metadata"] == "temp" and a["metric"] == "pi"), None)
    if pi:
        print(f"  sample diversity (pi) vs temp: rho={pi['rho']:.2f}, p={pi['p']:.2g}")
    n, bad = check_roundtrip(out + ".gfa", out + ".fasta")
    print(f"  GFA round trip: {n - len(bad)} of {n} sequences reproduced exactly")
    print(f"\nopen {out}.html in a browser to explore the graph")
    return {"hit": len(hit), "n_hot": len(truth["hotspots"]), "rank": rank, "bad": len(bad)}


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------
class Tests(unittest.TestCase):
    TOY = ">a count=5\nACGTA\n>b count=3\nACCTA\n>c count=2\nAC-TA\n"

    def graph(self, text=None):
        recs, _ = read_msa(io.StringIO(text or self.TOY))
        g0, _ = build_graph(recs)
        return recs, g0, compact(g0)

    def test_compaction_and_spelling(self):
        recs, g0, g = self.graph()
        self.assertEqual(len(g0.seq), 6)
        self.assertEqual(sorted(g.seq.values()), ["AC", "C", "G", "TA"])
        for r in recs:
            self.assertEqual(g.spell(r.name), r.seq.replace("-", ""))

    def test_bubble_alleles_and_entropy(self):
        recs, g0, g = self.graph()
        order = g.topo_order()
        c2r = col_to_ref(recs[0].seq)
        (b,) = collect_bubbles(g, order, "a", c2r)
        self.assertEqual(dict(b.alleles), {"G": 5, "C": 3, "": 2})
        self.assertEqual(b.kind, "indel")
        self.assertAlmostEqual(b.entropy, entropy([5, 3, 2]))
        self.assertEqual(b.pos, 3)

    def test_two_separate_snps_give_two_bubbles(self):
        _, _, g = self.graph(">a\nAACAAGAA\n>b\nAAGAAGAA\n>c\nAACAATAA\n")
        self.assertEqual(len(find_bubbles(g, g.topo_order())), 2)

    def test_tip_is_not_a_bubble(self):
        _, _, g = self.graph(">a\nACGTAC\n>b\nACGTAC\n>c\nACG---\n>d\nACCTAC\n")
        order = g.topo_order()
        for s, t in find_bubbles(g, order):
            self.assertTrue(g.out[t] or g.inn[t])

    def test_path_count_and_heaviest(self):
        _, _, g = self.graph(">a count=9\nACG\n>b count=1\nTCC\n")
        order = g.topo_order()
        self.assertEqual(count_paths(g, order), 4)
        best = heaviest_path(g, order, g.coverage())
        self.assertEqual("".join(g.seq[n] for n in best), "ACG")
        kb = kbest_paths(g, order, g.coverage(), 4)
        self.assertEqual(len(kb), 4)
        self.assertEqual(kb[0][1], best)

    def test_gfa_round_trip(self):
        recs, _, g = self.graph()
        buf = io.StringIO()
        write_gfa(g, buf)
        g2 = read_gfa(io.StringIO(buf.getvalue()))
        for r in recs:
            self.assertEqual(g2.spell(r.name), r.seq.replace("-", ""))
        self.assertEqual(g2.n_edges(), g.n_edges())
        self.assertEqual(g2.pmeta["a"], (5.0, "a"))

    def test_statistics(self):
        self.assertAlmostEqual(spearman([1, 2, 3, 4, 5], [2, 4, 6, 8, 10]), 1.0)
        self.assertAlmostEqual(spearman([1, 2, 3], [3, 2, 1]), -1.0)
        self.assertEqual(ranks([10, 20, 20, 30]), [1.0, 2.5, 2.5, 4.0])
        self.assertAlmostEqual(corr_pvalue(0.6, 10), 0.0667, delta=0.002)
        self.assertAlmostEqual(corr_pvalue(0.0, 20), 1.0, places=6)
        q = bh_adjust([0.01, 0.04, 0.03, 0.2])
        self.assertAlmostEqual(q[0], 0.04)
        self.assertAlmostEqual(q[2], 0.0533333, places=5)

    def test_codon_impact(self):
        recs, _ = read_msa(io.StringIO(">r\nATGAAACC\n>m\nATGAAGCC\n>n\nATGGAACC\n"))
        g0, _ = build_graph(recs)
        g = compact(g0)
        bs = collect_bubbles(g, g.topo_order(), "r", col_to_ref(recs[0].seq))
        annotate(bs, "ATGAAACC", [("g", 1, 6)])
        labels = " ".join(" ".join(b.impact) for b in bs)
        self.assertIn("A>G:synonymous(K>K)", labels)
        self.assertIn("A>G:missense(K>E)", labels)

    def test_simulation_is_deterministic(self):
        self.assertEqual(simulate(seed=3)[0], simulate(seed=3)[0])
        self.assertNotEqual(simulate(seed=3)[0], simulate(seed=4)[0])

    def test_end_to_end_recovers_truth(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            r = run_demo(os.path.join(d, "t"), seed=1)
        self.assertEqual(r["bad"], 0)
        self.assertGreaterEqual(r["hit"], r["n_hot"] - 2)
        self.assertIsNotNone(r["rank"])
        self.assertLessEqual(r["rank"], 3)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="Sequence variation graphs for viral populations.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="build and analyse a graph from an alignment")
    b.add_argument("alignment")
    b.add_argument("--out", required=True, help="output prefix")
    b.add_argument("--ref", help="reference record name (default: first record)")
    b.add_argument("--genes", help="genes TSV: gene, start, end (reference coordinates)")
    b.add_argument("--meta", help="metadata TSV with 'sample' as first column")
    b.add_argument("--k", type=int, default=5, help="number of best paths to report")
    b.add_argument("--window", type=int, default=30, help="hotspot window size (reference bp)")
    c = sub.add_parser("check", help="verify that a GFA spells the sequences of an alignment")
    c.add_argument("gfa")
    c.add_argument("alignment")
    s = sub.add_parser("simulate", help="write a simulated dataset with planted truth")
    s.add_argument("--out", required=True)
    s.add_argument("--length", type=int, default=600)
    s.add_argument("--samples", type=int, default=16)
    s.add_argument("--haps", type=int, default=12)
    s.add_argument("--seed", type=int, default=1)
    d = sub.add_parser("demo", help="simulate, build, and benchmark end to end")
    d.add_argument("--out", default="seqgraph_demo/demo")
    d.add_argument("--seed", type=int, default=1)
    sub.add_parser("test", help="run the unit tests")
    a = ap.parse_args(argv)
    if a.cmd == "build":
        run_build(a.alignment, a.out, a.ref, a.genes, a.meta, a.k, a.window)
    elif a.cmd == "check":
        n, bad = check_roundtrip(a.gfa, a.alignment)
        print(f"{n - len(bad)} of {n} sequences reproduced exactly")
        for name, why in bad:
            print(f"  {name}: {why}")
        return 1 if bad else 0
    elif a.cmd == "simulate":
        t = write_simulation(a.out, length=a.length, n_samples=a.samples, haps=a.haps, seed=a.seed)
        print(json.dumps(t))
    elif a.cmd == "demo":
        run_demo(a.out, a.seed)
    elif a.cmd == "test":
        res = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
        return 0 if res.wasSuccessful() else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())