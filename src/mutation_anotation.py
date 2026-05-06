import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os

# Output directory
os.makedirs("figures", exist_ok=True)

# Plotting style
sns.set_theme(style="whitegrid", palette="muted", font_scale=1.1)


def plot_mutation_types(mutations_df: pd.DataFrame, save_path: str = "figures/mutation_types.png") -> None:
    """
    Bar plot of high-impact mutation consequence types.
    Shows how many mutations fall into each VEP consequence category.
    """
    if mutations_df.empty:
        print("No mutations to plot.")
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    order = mutations_df["impact"].value_counts().index  # most frequent first

    sns.countplot(
        data=mutations_df,
        x="impact",
        order=order,
        palette="Blues_d",
        edgecolor="white",
        ax=ax
    )

    ax.set_title("High-Impact Mutation Consequence Types", fontsize=13, fontweight="bold")
    ax.set_xlabel("VEP Consequence", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.tick_params(axis="x", rotation=40)

    # Annotate bar heights
    for container in ax.containers:
        ax.bar_label(container, padding=3, fontsize=9)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved: {save_path}")


def plot_chromosomal_distribution(mutations_df: pd.DataFrame, save_path: str = "figures/chromosomal_distribution.png") -> None:
    """
    Bar plot of mutation counts per chromosome.
    Useful for identifying chromosomes with elevated somatic mutation burden.
    """
    if mutations_df.empty or "chrom" not in mutations_df.columns:
        print("No chromosomal data to plot.")
        return

    # Define standard chromosome order
    chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    present = [c for c in chrom_order if c in mutations_df["chrom"].values]

    counts = mutations_df["chrom"].value_counts().reindex(present, fill_value=0)

    fig, ax = plt.subplots(figsize=(14, 5))
    counts.plot(kind="bar", color=sns.color_palette("muted")[0], edgecolor="white", ax=ax)

    ax.set_title("Chromosomal Distribution of High-Impact Mutations", fontsize=13, fontweight="bold")
    ax.set_xlabel("Chromosome", fontsize=11)
    ax.set_ylabel("Mutation Count", fontsize=11)
    ax.tick_params(axis="x", rotation=0)

    for container in ax.containers:
        ax.bar_label(container, padding=3, fontsize=8)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved: {save_path}")


def plot_cadd_distribution(mutations_df: pd.DataFrame, save_path: str = "figures/cadd_distribution.png") -> None:
    """
    Histogram and KDE of CADD PHRED scores for annotated mutations.
    CADD >= 20 indicates top 1% most deleterious variants genome-wide.
    CADD >= 30 indicates top 0.1%.
    """
    cadd_available = mutations_df.dropna(subset=["cadd_score"])

    if cadd_available.empty:
        print("No CADD scores available to plot.")
        return

    fig, ax = plt.subplots(figsize=(9, 5))

    sns.histplot(
        cadd_available["cadd_score"],
        bins=30,
        kde=True,
        color=sns.color_palette("muted")[2],
        edgecolor="white",
        ax=ax
    )

    # Reference lines
    ax.axvline(20, color="orange", linestyle="--", linewidth=1.2, label="CADD ≥ 20 (top 1%)")
    ax.axvline(30, color="red",    linestyle="--", linewidth=1.2, label="CADD ≥ 30 (top 0.1%)")
    ax.legend(fontsize=9)

    ax.set_title("CADD PHRED Score Distribution of High-Impact Mutations", fontsize=13, fontweight="bold")
    ax.set_xlabel("CADD PHRED Score", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved: {save_path}")


def plot_pathway_enrichment(
    genes: list,
    gene_sets: list = ["KEGG_2021_Human", "Reactome_2022"],
    top_n: int = 15,
    save_path: str = "figures/pathway_enrichment.png"
) -> None:
    """
    Horizontal bar plot of enriched pathways from mutated genes.
    Uses -log10(adjusted P-value) so longer bars = more significant.
    Queries KEGG and Reactome by default.

    Parameters
    ----------
    genes     : list of Hugo gene symbols
    gene_sets : Enrichr library names to query
    top_n     : number of top pathways to display per library
    save_path : output file path
    """
    if not genes:
        print("No genes provided for pathway enrichment.")
        return

    try:
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=gene_sets,
            organism="human",
            outdir=None,          # suppress file output
            verbose=False
        )
    except Exception as e:
        print(f"Enrichr query failed: {e}")
        return

    results = enr.results.copy()

    if results.empty:
        print("No enrichment results returned.")
        return

    # Filter and compute -log10(adjusted P-value)
    results = results[results["Adjusted P-value"] < 0.05].copy()
    if results.empty:
        print("No pathways significant at FDR < 0.05.")
        return

    results["-log10(FDR)"] = -np.log10(results["Adjusted P-value"].clip(lower=1e-300))
    results = results.sort_values("-log10(FDR)", ascending=False)

    # Plot top_n per gene set
    gene_set_list = results["Gene_set"].unique()
    n_sets = len(gene_set_list)
    fig, axes = plt.subplots(1, n_sets, figsize=(9 * n_sets, 7), squeeze=False)

    for idx, gs in enumerate(gene_set_list):
        ax = axes[0][idx]
        subset = results[results["Gene_set"] == gs].head(top_n)

        # Colour by overlap ratio (genes hit / pathway size)
        subset = subset.copy()
        subset["Overlap_ratio"] = subset["Overlap"].apply(
            lambda x: int(x.split("/")[0]) / int(x.split("/")[1]) if "/" in str(x) else 0
        )

        bars = ax.barh(
            subset["Term"],
            subset["-log10(FDR)"],
            color=plt.cm.YlOrRd(subset["Overlap_ratio"] / subset["Overlap_ratio"].max()),
            edgecolor="white"
        )

        # Significance threshold line
        ax.axvline(-np.log10(0.05), color="grey", linestyle="--",
                   linewidth=1, label="FDR = 0.05")

        ax.set_xlabel("-log₁₀(Adjusted P-value)", fontsize=11)
        ax.set_title(f"Top {top_n} Enriched Pathways\n{gs}", fontsize=12, fontweight="bold")
        ax.invert_yaxis()  # most significant at top
        ax.legend(fontsize=9)

        # Colour bar for overlap ratio
        sm = plt.cm.ScalarMappable(
            cmap="YlOrRd",
            norm=plt.Normalize(vmin=0, vmax=subset["Overlap_ratio"].max())
        )
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, pad=0.02)
        cbar.set_label("Overlap Ratio\n(mutated / pathway size)", fontsize=9)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved: {save_path}")


def print_summary(mutations_df: pd.DataFrame) -> None:
    """Print a concise summary of the high-impact mutation set."""
    if mutations_df.empty:
        print("No high-impact mutations found.")
        return

    print("\n" + "=" * 50)
    print("HIGH-IMPACT MUTATION SUMMARY")
    print("=" * 50)
    print(f"Total high-impact mutations : {len(mutations_df)}")
    print(f"Unique genes affected       : {mutations_df['gene'].nunique()}")
    print(f"Chromosomes involved        : {sorted(mutations_df['chrom'].unique())}")
    print(f"\nConsequence breakdown:")
    print(mutations_df["impact"].value_counts().to_string())

    if mutations_df["cadd_score"].notna().any():
        cadd = mutations_df["cadd_score"].dropna()
        print(f"\nCADD score stats (n={len(cadd)}):")
        print(f"  Mean   : {cadd.mean():.1f}")
        print(f"  Median : {cadd.median():.1f}")
        print(f"  ≥ 20   : {(cadd >= 20).sum()} variants")
        print(f"  ≥ 30   : {(cadd >= 30).sum()} variants")
    print("=" * 50 + "\n")


# Main execution
# Load and filter high-impact mutations (update paths/filters as needed)
# Example: mutations_df = pd.read_csv("high_impact_mutations.csv")
mutations_df = pd.DataFrame()  # Replace with actual data loading

if not mutations_df.empty:

    # 1. Summary statistics
    print_summary(mutations_df)

    # 2. Consequence type distribution
    plot_mutation_types(mutations_df)

    # 3. Chromosomal distribution
    plot_chromosomal_distribution(mutations_df)

    # 4. CADD score distribution (only if scores were loaded)
    plot_cadd_distribution(mutations_df)

    # 5. Pathway enrichment on mutated genes
    genes = mutations_df["gene"].replace("Unknown", pd.NA).dropna().unique().tolist()
    plot_pathway_enrichment(genes)

else:
    print("No high-impact mutations found. Check VCF annotations (CSQ field) and filtering thresholds.")