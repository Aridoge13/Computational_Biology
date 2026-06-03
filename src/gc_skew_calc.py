from Bio import SeqIO, Entrez
import numpy as np
import matplotlib.pyplot as plt
import os

# Configuration
Entrez.email = "example@email.com"

GENOMES = {
    # Asgard archaea
    "Lokiarchaeota": {"accession": "GCA_001940045.1", "group": "Asgard"},
    "Thorarchaeota": {"accession": "GCA_002775275.1", "group": "Asgard"},
    "Heimdallarchaeota": {"accession": "GCA_002728275.1", "group": "Asgard"},
    # TACK archaea
    "Sulfolobus_solfataricus": {"accession": "GCA_000007005.1", "group": "TACK"},
    "Nitrosopumilus_maritimus": {"accession": "GCA_000018465.1", "group": "TACK"},
    "Pyrobaculum_aerophilum": {"accession": "GCA_000007225.1", "group": "TACK"},
    # Euryarchaeota
    "Methanobacterium_thermo": {"accession": "GCA_000008545.1", "group": "Euryarchaeota"},
    "Haloarcula_marismortui": {"accession": "GCA_000011085.1", "group": "Euryarchaeota"},
    "Archaeoglobus_fulgidus": {"accession": "GCA_000007685.1", "group": "Euryarchaeota"},
    # Bacteria (reference comparison)
    "Escherichia_coli_K12": {"accession": "GCA_000005845.2", "group": "Bacteria"},
    "Bacillus_subtilis": {"accession": "GCA_000009045.1", "group": "Bacteria"},
}

GROUP_COLORS = {
    "Asgard": "#E63946",
    "TACK": "#2A9D8F",
    "Euryarchaeota": "#E9C46A",
    "Bacteria": "#AAAAAA",
}

GENOME_DIR = "genomes"
FIGURES_DIR = "figures"
os.makedirs(GENOME_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


# Core Functions

def gc_content(seq):
    """Calculate GC content of a DNA sequence."""
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0


def gc_skew(seq, window=10000, step=5000):
    """Compute GC skew in sliding windows across a sequence."""
    skews = []
    positions = []
    seq = str(seq).upper()
    for i in range(0, len(seq) - window, step):
        window_seq = seq[i:i + window]
        g = window_seq.count('G')
        c = window_seq.count('C')
        denom = g + c
        skew = (g - c) / denom if denom > 0 else 0
        skews.append(skew)
        positions.append(i + window // 2)
    return np.array(positions), np.array(skews)


def cumulative_gc_skew(skews):
    """Compute cumulative GC skew — reveals origin/terminus positions."""
    return np.cumsum(skews)


# Data Acquisition

def download_genome(name, accession, genome_dir=GENOME_DIR):
    """
    Download genome FASTA from NCBI if not already present.
    Returns path to the FASTA file.
    """
    fasta_path = os.path.join(genome_dir, f"{name}.fasta")
    if os.path.exists(fasta_path):
        print(f"  [cached] {name}")
        return fasta_path

    print(f"  [downloading] {name} ({accession})...")
    try:
        search_handle = Entrez.esearch(db="assembly", term=accession)
        search_record = Entrez.read(search_handle)
        search_handle.close()

        # Fetch via nuccore linked to assembly
        link_handle = Entrez.elink(
            dbfrom="assembly",
            db="nuccore",
            id=search_record["IdList"][0],
            linkname="assembly_nuccore_refseq"
        )
        link_record = Entrez.read(link_handle)
        link_handle.close()

        nuc_ids = [
            link["Id"]
            for link in link_record[0]["LinkSetDb"][0]["Link"]
        ]

        # Fetch only the largest sequence (main chromosome)
        fetch_handle = Entrez.efetch(
            db="nuccore",
            id=nuc_ids[0],
            rettype="fasta",
            retmode="text"
        )
        with open(fasta_path, "w") as f:
            f.write(fetch_handle.read())
        fetch_handle.close()
        print(f"  [saved] {fasta_path}")
        return fasta_path

    except Exception as e:
        print(f"  [error] Could not download {name}: {e}")
        return None


def load_primary_sequence(fasta_path):
    """
    Load the longest sequence from a FASTA file.
    For multi-chromosome genomes, this returns the primary chromosome.
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return None
    # Return longest sequence — most likely the main chromosome
    return max(records, key=lambda r: len(r.seq))


# Analysis

def analyse_genome(name, fasta_path):
    """
    Run GC content, GC skew, and cumulative skew analysis for one genome.
    Returns a results dictionary.
    """
    record = load_primary_sequence(fasta_path)
    if record is None:
        print(f"  [warning] No sequence found for {name}")
        return None

    seq = str(record.seq).upper()
    genome_len = len(seq)

    positions, skews = gc_skew(seq)
    cum_skews = cumulative_gc_skew(skews)

    # Normalise positions to genome fraction for cross-species comparison
    norm_positions = positions / genome_len

    # Summary statistics
    gc = gc_content(seq)
    skew_amplitude = cum_skews.max() - cum_skews.min()
    predicted_origin_fraction = norm_positions[np.argmin(cum_skews)]

    return {
        "name": name,
        "genome_length": genome_len,
        "gc_content": gc,
        "norm_positions": norm_positions,
        "skews": skews,
        "cum_skews": cum_skews,
        "skew_amplitude": skew_amplitude,
        "predicted_origin_fraction": predicted_origin_fraction,
    }


# Visualisation

def plot_cumulative_skew(results_list, genome_metadata):
    """
    Figure 1: Cumulative GC skew across all genomes.
    One subplot per genome, coloured by lineage.
    """
    n = len(results_list)
    ncols = 3
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(15, nrows * 3.5))
    axes = axes.flatten()

    for idx, result in enumerate(results_list):
        ax = axes[idx]
        name = result["name"]
        group = genome_metadata[name]["group"]
        color = GROUP_COLORS[group]

        ax.plot(
            result["norm_positions"],
            result["cum_skews"],
            color=color,
            linewidth=1.5
        )
        ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_title(
            f"{name.replace('_', ' ')}\n({group})",
            fontsize=8,
            fontweight="bold"
        )
        ax.set_xlabel("Genome position (fraction)", fontsize=7)
        ax.set_ylabel("Cumulative GC skew", fontsize=7)
        ax.tick_params(labelsize=7)

    # Hide unused subplots
    for idx in range(len(results_list), len(axes)):
        axes[idx].set_visible(False)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=color, label=group)
        for group, color in GROUP_COLORS.items()
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower right",
        fontsize=9,
        title="Lineage"
    )

    fig.suptitle(
        "Cumulative GC Skew Across Archaeal and Bacterial Genomes",
        fontsize=13,
        fontweight="bold",
        y=1.01
    )
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, "figure1_cumulative_gc_skew.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    print(f"[saved] {path}")
    plt.close()


def plot_summary_scatter(results_list, genome_metadata):
    """
    Figure 2: Skew amplitude vs predicted origin position.
    Each point is one genome, coloured by lineage, labelled by name.
    """
    fig, ax = plt.subplots(figsize=(10, 7))

    for result in results_list:
        name = result["name"]
        group = genome_metadata[name]["group"]
        color = GROUP_COLORS[group]

        ax.scatter(
            result["predicted_origin_fraction"],
            result["skew_amplitude"],
            color=color,
            s=100,
            zorder=3,
            edgecolors="black",
            linewidths=0.5
        )
        ax.annotate(
            name.replace("_", " "),
            (result["predicted_origin_fraction"], result["skew_amplitude"]),
            fontsize=7,
            xytext=(5, 3),
            textcoords="offset points"
        )

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=color, label=group)
        for group, color in GROUP_COLORS.items()
    ]
    ax.legend(handles=legend_elements, title="Lineage", fontsize=9)

    ax.set_xlabel(
        "Predicted origin position (genome fraction)",
        fontsize=11
    )
    ax.set_ylabel(
        "Cumulative skew amplitude\n(proxy for replication symmetry)",
        fontsize=11
    )
    ax.set_title(
        "Replication-Associated GC Skew: Summary Across Lineages",
        fontsize=12,
        fontweight="bold"
    )
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, "figure2_skew_summary.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"[saved] {path}")
    plt.close()


def print_summary_table(results_list, genome_metadata):
    """Print a summary table of key metrics for all genomes."""
    print(f"\n{'Name':<35} {'Group':<15} {'Length (Mb)':<12} "
          f"{'GC%':<8} {'Amplitude':<12} {'Pred. Origin':<12}")
    print("-" * 95)
    for r in results_list:
        group = genome_metadata[r['name']]['group']
        print(
            f"{r['name'].replace('_', ' '):<35} "
            f"{group:<15} "
            f"{r['genome_length']/1e6:<12.2f} "
            f"{r['gc_content']:<8.1f} "
            f"{r['skew_amplitude']:<12.3f} "
            f"{r['predicted_origin_fraction']:<12.3f}"
        )


# Main

if __name__ == "__main__":
    print("=== Archaeal GC Skew Analysis ===\n")

    # Download all genomes
    print("[1] Acquiring genomes...")
    fasta_paths = {}
    for name, meta in GENOMES.items():
        path = download_genome(name, meta["accession"])
        if path:
            fasta_paths[name] = path

    # Run analysis
    print("\n[2] Running GC skew analysis...")
    results = []
    for name, path in fasta_paths.items():
        print(f"  Analysing {name}...")
        result = analyse_genome(name, path)
        if result:
            results.append(result)

    # Print summary
    print("\n[3] Summary statistics:")
    print_summary_table(results, GENOMES)

    # Generate figures
    print("\n[4] Generating figures...")
    plot_cumulative_skew(results, GENOMES)
    plot_summary_scatter(results, GENOMES)

    print("\n=== Analysis complete ===")
    print(f"Figures saved to ./{FIGURES_DIR}/")