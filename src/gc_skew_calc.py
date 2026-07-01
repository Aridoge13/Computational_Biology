from Bio import SeqIO, Entrez
import numpy as np
import matplotlib.pyplot as plt
import os

# Configuration
Entrez.email = "aritra.mukherjee98@gmail.com"  # Change to your email

GENOMES = {
    # ASGARD
    "Lokiarchaeota": {"accession": "GCA_001940045.1", "group": "Asgard"},
    "Thorarchaeota": {"accession": "GCA_002775275.1", "group": "Asgard"},
    "Heimdallarchaeota": {"accession": "GCA_002728275.1", "group": "Asgard"},
    "Odinarchaeota_LCB_4": {"accession": "GCA_002499645.1", "group": "Asgard"},
    "Prometheoarchaeum_syntrophicum": {"accession": "GCA_014524485.1", "group": "Asgard"},
    "Candidatus_Hodarchaeia": {"accession": "GCA_017607445.1", "group": "Asgard"},
    "Candidatus_Jordarchaeia": {"accession": "GCA_017607455.1", "group": "Asgard"},
    "Hermodarchaeota": {"accession": "GCA_012841095.1", "group": "Asgard"},
    "Baldrarchaeota": {"accession": "GCA_017607435.1", "group": "Asgard"},
    "Tyrarchaeota": {"accession": "GCA_017607425.1", "group": "Asgard"},

    # DPANN
    "Nanoarchaeum_equitans": {"accession": "GCA_000091665.1", "group": "DPANN"},
    "Candidatus_Parvarchaeum_acidophilus_ARMAN5": {"accession": "GCA_002412085.1", "group": "DPANN"},
    "archaeon_GW2011_AR5": {"accession": "GCA_000806115.1", "group": "DPANN"},
    "Candidatus_Micrarchaeum_harzensis": {"accession": "GCA_000819665.1", "group": "DPANN"},
    "Candidatus_Mancarchaeum_acidiphilum": {"accession": "GCA_001563885.1", "group": "DPANN"},
    "Candidatus_Woesearchaeota_archaeon_AR20": {"accession": "GCA_002495055.1", "group": "DPANN"},
    "Candidatus_Pacearchaeota_archaeon": {"accession": "GCA_002495065.1", "group": "DPANN"},
    "Candidatus_Diapherotrites_archaeon": {"accession": "GCA_002495075.1", "group": "DPANN"},
    "Candidatus_Aenigmarchaeota_archaeon": {"accession": "GCA_002495085.1", "group": "DPANN"},
    "Candidatus_Huberarchaeum_crystalense": {"accession": "GCA_002503445.1", "group": "DPANN"},

    # TACK
    "Sulfolobus_solfataricus": {"accession": "GCA_000007005.1", "group": "TACK"},
    "Sulfolobus_acidocaldarius": {"accession": "GCA_000009965.1", "group": "TACK"},
    "Pyrobaculum_aerophilum": {"accession": "GCA_000007225.1", "group": "TACK"},
    "Nitrosopumilus_maritimus": {"accession": "GCA_000018465.1", "group": "TACK"},
    "Ignicoccus_hospitalis": {"accession": "GCA_000021665.1", "group": "TACK"},
    "Thermoproteus_tenax": {"accession": "GCA_000021505.1", "group": "TACK"},
    "Metallosphaera_sedula": {"accession": "GCA_000011205.1", "group": "TACK"},
    "Aeropyrum_pernix": {"accession": "GCA_000011125.1", "group": "TACK"},
    "Thermofilum_pendens": {"accession": "GCA_000024705.1", "group": "TACK"},
    "Caldivirga_maquilingensis": {"accession": "GCA_000204585.1", "group": "TACK"},

    # EURYARCHAEOTA
    "Methanobacterium_thermo": {"accession": "GCA_000008545.1", "group": "Euryarchaeota"},
    "Haloarcula_marismortui": {"accession": "GCA_000011085.1", "group": "Euryarchaeota"},
    "Archaeoglobus_fulgidus": {"accession": "GCA_000007685.1", "group": "Euryarchaeota"},
    "Methanocaldococcus_jannaschii": {"accession": "GCA_000091665.2", "group": "Euryarchaeota"},
    "Methanosarcina_mazei": {"accession": "GCA_000007345.1", "group": "Euryarchaeota"},
    "Methanococcus_maripaludis": {"accession": "GCA_000011585.1", "group": "Euryarchaeota"},
    "Methanopyrus_kandleri": {"accession": "GCA_000011725.1", "group": "Euryarchaeota"},
    "Halobacterium_salinarum": {"accession": "GCA_000006805.1", "group": "Euryarchaeota"},
    "Methanobrevibacter_smithii": {"accession": "GCA_000016525.1", "group": "Euryarchaeota"},
    "Thermococcus_kodakarensis": {"accession": "GCA_000009965.2", "group": "Euryarchaeota"},

    # REDUCED BACTERIAL SYMBIONTS
    "Buchnera_aphidicola_APS": {"accession": "GCA_000009605.1", "group": "Symbiont"},
    "Wigglesworthia_glossinidia": {"accession": "GCA_000009785.1", "group": "Symbiont"},
    "Blochmannia_floridanus": {"accession": "GCA_000009625.1", "group": "Symbiont"},
    "Carsonella_ruddii": {"accession": "GCA_000012245.1", "group": "Symbiont"},
    "Tremblaya_princeps": {"accession": "GCA_000014525.1", "group": "Symbiont"},
    "Hodgkinia_cicadicola": {"accession": "GCA_000155385.1", "group": "Symbiont"},
    "Baumannia_cicadellinicola": {"accession": "GCA_000009845.1", "group": "Symbiont"},
    "Portiera_aleyrodidarum": {"accession": "GCA_000147015.1", "group": "Symbiont"},
    "Nasuia_deltocephalinicola": {"accession": "GCA_000334075.1", "group": "Symbiont"},
    "Sulcia_muelleri": {"accession": "GCA_000010105.1", "group": "Symbiont"},
}


GROUP_COLORS = {
    "DPANN": "#9B5DE5",
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
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0

def gc_skew(seq, window=10000, step=5000):
    seq = str(seq).upper()
    if len(seq) < window:
        window = max(1, len(seq) // 2)
        step = window // 2 if window > 1 else 1
        print(f"    [notice] Sequence length {len(seq)} < default window; using window={window}")

    skews = []
    positions = []
    for i in range(0, len(seq) - window, step):
        window_seq = seq[i:i + window]
        g = window_seq.count('G')
        c = window_seq.count('C')
        denom = g + c
        skew = (g - c) / denom if denom > 0 else 0
        skews.append(skew)
        positions.append(i + window // 2)

    if not positions:
        positions = [0]
        skews = [0]
    return np.array(positions), np.array(skews)

def cumulative_gc_skew(skews):
    return np.cumsum(skews)

# Data Acquisition

def download_genome(name, accession, genome_dir=GENOME_DIR):
    fasta_path = os.path.join(genome_dir, f"{name}.fasta")
    if os.path.exists(fasta_path):
        print(f"  [cached] {name}")
        return fasta_path

    print(f"  [downloading] {name} ({accession})...")
    try:
        # Directly fetch the full assembly from the Assembly database
        fetch_handle = Entrez.efetch(
            db="assembly",
            id=accession,
            rettype="fasta",
            retmode="text"
        )
        content = fetch_handle.read()
        fetch_handle.close()

        if not content.strip():
            raise ValueError("Empty FASTA content")

        with open(fasta_path, "w") as f:
            f.write(content)
        print(f"  [saved] {fasta_path}")
        return fasta_path

    except Exception as e:
        print(f"  [error] Could not download {name}: {e}")
        return None

def load_primary_sequence(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return None
    return max(records, key=lambda r: len(r.seq))

# Analysis

def analyse_genome(name, fasta_path):
    record = load_primary_sequence(fasta_path)
    if record is None:
        print(f"  [warning] No sequence found for {name}")
        return None

    seq = str(record.seq).upper()
    genome_len = len(seq)
    if genome_len < 100:
        print(f"  [warning] Genome too short ({genome_len} bp) for meaningful skew; skipping")
        return None

    positions, skews = gc_skew(seq)
    cum_skews = cumulative_gc_skew(skews)

    if cum_skews.size == 0:
        print(f"  [warning] No skew values computed for {name}; skipping")
        return None

    norm_positions = positions / genome_len
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
    results_list = [r for r in results_list if r is not None]
    if not results_list:
        print("No valid results to plot cumulative skew.")
        return

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

        ax.plot(result["norm_positions"], result["cum_skews"],
                color=color, linewidth=1.5)
        ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_title(f"{name.replace('_', ' ')}\n({group})", fontsize=8, fontweight="bold")
        ax.set_xlabel("Genome position (fraction)", fontsize=7)
        ax.set_ylabel("Cumulative GC skew", fontsize=7)
        ax.tick_params(labelsize=7)

    for idx in range(len(results_list), len(axes)):
        axes[idx].set_visible(False)

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, label=group) for group, color in GROUP_COLORS.items()]
    fig.legend(handles=legend_elements, loc="lower right", fontsize=9, title="Lineage")

    fig.suptitle("Cumulative GC Skew Across Archaeal and Bacterial Genomes", fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, "figure1_cumulative_gc_skew.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    print(f"[saved] {path}")
    plt.close()

def plot_summary_scatter(results_list, genome_metadata):
    results_list = [r for r in results_list if r is not None]
    if not results_list:
        print("No valid results to plot summary scatter.")
        return

    fig, ax = plt.subplots(figsize=(10, 7))
    for result in results_list:
        name = result["name"]
        group = genome_metadata[name]["group"]
        color = GROUP_COLORS[group]

        ax.scatter(result["predicted_origin_fraction"], result["skew_amplitude"],
                   color=color, s=100, zorder=3, edgecolors="black", linewidths=0.5)
        ax.annotate(name.replace("_", " "),
                    (result["predicted_origin_fraction"], result["skew_amplitude"]),
                    fontsize=7, xytext=(5, 3), textcoords="offset points")

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, label=group) for group, color in GROUP_COLORS.items()]
    ax.legend(handles=legend_elements, title="Lineage", fontsize=9)

    ax.set_xlabel("Predicted origin position (genome fraction)", fontsize=11)
    ax.set_ylabel("Cumulative skew amplitude\n(proxy for replication symmetry)", fontsize=11)
    ax.set_title("Replication-Associated GC Skew: Summary Across Lineages", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, "figure2_skew_summary.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"[saved] {path}")
    plt.close()

def print_summary_table(results_list, genome_metadata):
    results_list = [r for r in results_list if r is not None]
    if not results_list:
        print("No valid results to print.")
        return

    print(f"\n{'Name':<35} {'Group':<15} {'Length (Mb)':<12} "
          f"{'GC%':<8} {'Amplitude':<12} {'Pred. Origin':<12}")
    print("-" * 95)
    for r in results_list:
        group = genome_metadata[r['name']]['group']
        print(f"{r['name'].replace('_', ' '):<35} "
              f"{group:<15} "
              f"{r['genome_length']/1e6:<12.2f} "
              f"{r['gc_content']:<8.1f} "
              f"{r['skew_amplitude']:<12.3f} "
              f"{r['predicted_origin_fraction']:<12.3f}")

# Main

if __name__ == "__main__":
    print("=== Archaeal GC Skew Analysis ===\n")
    print("[1] Acquiring genomes...")
    fasta_paths = {}
    for name, meta in GENOMES.items():
        path = download_genome(name, meta["accession"])
        if path:
            fasta_paths[name] = path

    print("\n[2] Running GC skew analysis...")
    results = []
    for name, path in fasta_paths.items():
        print(f"  Analysing {name}...")
        result = analyse_genome(name, path)
        if result:
            results.append(result)

    print("\n[3] Summary statistics:")
    print_summary_table(results, GENOMES)

    print("\n[4] Generating figures...")
    plot_cumulative_skew(results, GENOMES)
    plot_summary_scatter(results, GENOMES)

    print("\n=== Analysis complete ===")
    print(f"Figures saved to ./{FIGURES_DIR}/")