from Bio import Entrez, SeqIO
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# CONFIGURATION


Entrez.email = "your_email@example.com"

GENOMES = {
    # DPANN-like / streamlined
    "Nanoarchaeum_equitans": {
        "accession": "GCA_000091665.1",
        "group": "DPANN"
    },

    # Asgard
    "Lokiarchaeota": {
        "accession": "GCA_001940045.1",
        "group": "Asgard"
    },

    # Free-living archaea
    "Sulfolobus_solfataricus": {
        "accession": "GCA_000007005.1",
        "group": "TACK"
    },

    # Bacteria
    "Escherichia_coli_K12": {
        "accession": "GCA_000005845.2",
        "group": "Bacteria"
    }
}

GROUP_COLORS = {
    "DPANN": "#D62728",
    "Asgard": "#1F77B4",
    "TACK": "#2CA02C",
    "Bacteria": "#7F7F7F"
}

GENOME_DIR = "genomes"
FIG_DIR = "figures"

os.makedirs(GENOME_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


# DOWNLOAD GENOMES

def download_genome(name, accession):

    fasta_path = os.path.join(GENOME_DIR, f"{name}.fna")

    if os.path.exists(fasta_path):
        print(f"[cached] {name}")
        return fasta_path

    try:
        search = Entrez.esearch(db="assembly", term=accession)
        search_record = Entrez.read(search)
        assembly_id = search_record["IdList"][0]

        links = Entrez.elink(
            dbfrom="assembly",
            db="nuccore",
            id=assembly_id,
            linkname="assembly_nuccore_refseq"
        )

        link_record = Entrez.read(links)

        nuc_id = link_record[0]["LinkSetDb"][0]["Link"][0]["Id"]

        fetch = Entrez.efetch(
            db="nuccore",
            id=nuc_id,
            rettype="fasta_cds_na",
            retmode="text"
        )

        with open(fasta_path, "w") as f:
            f.write(fetch.read())

        print(f"[downloaded] {name}")

        return fasta_path

    except Exception as e:
        print(f"[ERROR] {name}: {e}")
        return None


# CORE ANALYSIS

def gc_content(seq):

    seq = seq.upper()

    gc = seq.count("G") + seq.count("C")

    return (gc / len(seq)) * 100 if len(seq) > 0 else 0


def analyse_genome(name, fasta_path):

    records = list(SeqIO.parse(fasta_path, "fasta"))

    if len(records) == 0:
        return None

    cds_lengths = []
    gc_values = []

    hypothetical_count = 0
    atpase_count = 0

    total_bp = 0

    for record in records:

        seq = str(record.seq)

        total_bp += len(seq)

        cds_lengths.append(len(seq))

        gc_values.append(gc_content(seq))

        desc = record.description.lower()

        if "hypothetical protein" in desc:
            hypothetical_count += 1

        if "atp synthase" in desc:
            atpase_count += 1

    n_genes = len(records)

    coding_density = total_bp / total_bp

    avg_protein_length = np.mean(cds_lengths)

    hypothetical_fraction = hypothetical_count / n_genes

    avg_gc = np.mean(gc_values)

    
    # GENOME DEPENDENCY INDEX (heuristic)
    # Higher = more streamlined / host-dependent

    gdi = (
        (1 / np.log10(total_bp + 10))
        * 100
        + hypothetical_fraction * 25
        - atpase_count * 2
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
        "GDI": gdi
    }


# PCA VISUALISATION

def plot_pca(df):

    features = [
        "genome_size_mb",
        "genes",
        "avg_gc",
        "avg_gene_length",
        "hypothetical_fraction",
        "atpase_genes",
        "GDI"
    ]

    X = df[features]

    scaler = StandardScaler()

    X_scaled = scaler.fit_transform(X)

    pca = PCA(n_components=2)

    coords = pca.fit_transform(X_scaled)

    df["PC1"] = coords[:, 0]
    df["PC2"] = coords[:, 1]

    fig, ax = plt.subplots(figsize=(9, 7))

    for _, row in df.iterrows():

        ax.scatter(
            row["PC1"],
            row["PC2"],
            s=120,
            color=GROUP_COLORS[row["group"]],
            edgecolors="black"
        )

        ax.text(
            row["PC1"] + 0.05,
            row["PC2"] + 0.05,
            row["name"].replace("_", " "),
            fontsize=8
        )

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

    ax.set_title(
        "Comparative Genome Reduction Landscape",
        fontsize=13,
        fontweight="bold"
    )

    plt.tight_layout()

    path = os.path.join(FIG_DIR, "genome_reduction_pca.png")

    plt.savefig(path, dpi=300)

    print(f"[saved] {path}")


# SUMMARY TABLE

def print_summary(df):

    print("\n=== Genome Reduction Summary ===\n")

    print(df.to_string(index=False))


# MAIN

if __name__ == "__main__":

    print("\n=== Comparative Genome Reduction Pipeline ===\n")

    fasta_paths = {}

    for name, meta in GENOMES.items():

        path = download_genome(name, meta["accession"])

        if path:
            fasta_paths[name] = path

    results = []

    for name, path in fasta_paths.items():

        print(f"\n[analysing] {name}")

        r = analyse_genome(name, path)

        if r:
            results.append(r)

    df = pd.DataFrame(results)

    print_summary(df)

    plot_pca(df)

    csv_path = "genome_reduction_results.csv"

    df.to_csv(csv_path, index=False)

    print(f"\n[saved] {csv_path}")

    print("\n=== Pipeline complete ===")

