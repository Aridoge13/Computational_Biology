import os
import urllib.request
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from adjustText import adjust_text

# ==================== CONFIGURATION ====================
Entrez.email = "aritra.mukherjee98@gmail.com"

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

    # EURYARCHAEOTA (duplicates removed)
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
    "Symbiont": "#AAAAAA",
}

GENOME_DIR = "genomes"
FIG_DIR = "figures"

os.makedirs(GENOME_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


# ==================== 1. ROBUST DOWNLOAD ====================
def download_genome(name, accession):
    """
    Download full assembly FASTA from NCBI using the GenBank FTP path.
    """
    fasta_path = os.path.join(GENOME_DIR, f"{name}.fna")

    if os.path.exists(fasta_path):
        print(f"[cached] {name}")
        return fasta_path

    try:
        # Search for assembly ID
        search = Entrez.esearch(db="assembly", term=accession)
        search_record = Entrez.read(search)
        if not search_record["IdList"]:
            raise ValueError(f"No assembly found for: {accession}")

        assembly_id = search_record["IdList"][0]
        summary = Entrez.esummary(db="assembly", id=assembly_id)
        summary_record = Entrez.read(summary)

        # Prefer GenBank FTP path, fallback to RefSeq
        ftp_path = summary_record["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_GenBank"]
        if not ftp_path:
            ftp_path = summary_record["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]

        if not ftp_path:
            raise ValueError(f"No remote FTP path available for assembly ID: {assembly_id}")

        # Build download URL
        assembly_name = os.path.basename(ftp_path)
        download_url = f"{ftp_path}/{assembly_name}_genomic.fna.gz"

        gz_path = fasta_path + ".gz"
        print(f"[fetching FTP archive] {name} -> {download_url}")

        # Download and decompress
        req = urllib.request.Request(download_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response, open(gz_path, 'wb') as out_file:
            out_file.write(response.read())

        with gzip.open(gz_path, "rt") as f_in, open(fasta_path, "w") as f_out:
            f_out.write(f_in.read())

        os.remove(gz_path)
        print(f"[downloaded full assembly] {name}")
        return fasta_path

    except Exception as e:
        print(f"[ERROR downloading] {name}: {e}")
        return None


# ==================== 2. ANALYSIS CORE ====================
def gc_content(seq):
    """Calculate GC percentage for a nucleotide sequence."""
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0


def analyse_genome(name, fasta_path):
    """
    Parse a FASTA file (full assembly) and extract genome‑level metrics.
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) == 0:
        return None

    # Total assembly length
    total_bp = sum(len(record.seq) for record in records)

    # Global GC content over the whole assembly
    full_sequence = "".join(str(record.seq) for record in records)
    avg_gc = gc_content(full_sequence)

    # Count features from CDS/FASTA headers
    hypothetical_count = 0
    atpase_count = 0
    n_genes = 0
    cds_lengths = []

    for record in records:
        desc = record.description.lower()
        if "hypothetical" in desc:
            hypothetical_count += 1
        if "atp" in desc or "atpase" in desc or "synthase" in desc:
            atpase_count += 1
        if len(record.seq) > 0:
            n_genes += 1
            cds_lengths.append(len(record.seq))

    if n_genes == 0:
        n_genes = 1  # safety

    avg_protein_length = np.mean(cds_lengths) if cds_lengths else 0
    hypothetical_fraction = hypothetical_count / n_genes
    coding_density_proxy = (n_genes * avg_protein_length) / total_bp if total_bp > 0 else 0

    # Genome Dependency Index (exploratory heuristic)
    gdi = (
        (1 / np.log10(total_bp + 10)) * 100
        + (hypothetical_fraction * 25)
        - (1 if atpase_count == 0 else 0) * 10
    )

    return {
        "name": name,
        "group": GENOMES[name]["group"],
        "genome_size_mb": total_bp / 1e6,
        "genes": n_genes,
        "avg_gc": avg_gc,
        "avg_gene_length": avg_protein_length,
        "hypothetical_fraction": hypothetical_fraction,
        "atpase_genes": atpase_count,
        "coding_density": coding_density_proxy,
        "GDI": gdi
    }


# ==================== 3. PCA WITH ADJUSTED TEXT ====================
def plot_pca(df):
    """
    Perform PCA on selected genomic features and generate a labelled scatter plot.
    Labels are repelled using adjustText to avoid overlap.
    """
    # Features used for PCA (GDI excluded to avoid data leakage)
    features = [
        "genome_size_mb",
        "genes",
        "avg_gc",
        "avg_gene_length",
        "hypothetical_fraction",
        "atpase_genes"
    ]

    X = df[features]
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(X_scaled)

    df["PC1"] = coords[:, 0]
    df["PC2"] = coords[:, 1]

    fig, ax = plt.subplots(figsize=(11, 9))

    # Colour mapping
    colors = df["group"].map(GROUP_COLORS).fillna("#AAAAAA")
    ax.scatter(
        df["PC1"],
        df["PC2"],
        s=140,
        c=colors,
        edgecolors="black",
        zorder=3
    )

    # Add labels and collect them for adjustText
    texts = []
    for _, row in df.iterrows():
        clean_name = row["name"].replace("_", " ")
        texts.append(ax.text(row["PC1"], row["PC2"], clean_name,
                             fontsize=8, fontweight='medium'))

    # Repel labels to reduce overlap
    adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

    # Axis labels with variance explained
    var_exp = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}% Variance Explained)", fontsize=11)
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}% Variance Explained)", fontsize=11)
    ax.set_title("Comparative Genome Reduction Landscape\n(Pure Feature PCA – Corrected Dataset Mapping)",
                 fontsize=13, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.3, zorder=1)

    # Print feature loadings to console for reference
    print("\n=== Feature Loadings (Eigenvectors) ===")
    loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=features)
    print(loadings.to_string())

    plt.tight_layout()
    path = os.path.join(FIG_DIR, "genome_reduction_pca_clean.png")
    plt.savefig(path, dpi=300)
    print(f"\n[saved cluster plot] {path}")


# ==================== 4. SUMMARY TABLE ====================
def print_summary(df):
    """Print a clean summary table of all analysed genomes."""
    print("\n=== Genome Reduction Summary ===\n")
    print(df.to_string(index=False))


# ==================== MAIN ====================
if __name__ == "__main__":
    print("\n=== Comparative Genome Reduction Pipeline ===")

    # Download all genomes
    fasta_paths = {}
    for name, meta in GENOMES.items():
        path = download_genome(name, meta["accession"])
        if path:
            fasta_paths[name] = path

    # Analyse each
    results = []
    for name, path in fasta_paths.items():
        print(f"\n[analysing] {name}")
        r = analyse_genome(name, path)
        if r:
            results.append(r)

    if results:
        df = pd.DataFrame(results)
        print_summary(df)
        plot_pca(df)

        # Save results to CSV
        csv_path = "genome_reduction_results.csv"
        df.to_csv(csv_path, index=False)
        print(f"[saved summary spreadsheet] {csv_path}")

    print("\n=== Pipeline Run Finalized ===")