import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def load_kraken_report(file_path):
    "load kraken report file into a dataframe"
    return pd.read_csv(file_path, sep="\t", header=None, names=["Percentage", "Reads", "Taxon"])


def load_metaphlan_output(file_path):
    "load metaphlan output file into a dataframe"
    return pd.read_csv(file_path, sep="\t", skiprows=1, names=["Taxon", "Relative_Abundance"])


def plot_kraken_results(df):
    "plot Kraken results as a bar chart"
    df = df[df["Percentage"] > 1]
    plt.figure(figsize=(10, 6))
    plt.barh(df["Taxon"], df["Percentage"], color="skyblue")
    plt.xlabel("Percentage of Reads")
    plt.ylabel("Taxon")
    plt.title("Taxonomic Classification (Kraken)")
    plt.tight_layout()
    plt.show()


def plot_kraken_results_pie(df):
    "plot Kraken results as a pie chart"
    df = df[df["Percentage"] > 1]
    plt.figure(figsize=(10, 6))
    plt.pie(df["Percentage"], labels=df["Taxon"], autopct="%1.1f%%")
    plt.title("Taxonomic Classification (Kraken)")
    plt.tight_layout()
    plt.show()


def plot_kraken_results_heatmap(df):
    "plot Kraken results as a heatmap"
    df = df[df["Percentage"] > 1]
    df["Taxon"] = df["Taxon"].apply(lambda x: x.split("|")[-1])
    df = df.pivot(index="Taxon", columns="Reads", values="Percentage")
    plt.figure(figsize=(10, 6))
    sns.heatmap(df, cmap="YlGnBu")
    plt.title("Taxonomic Classification (Kraken)")
    plt.tight_layout()
    plt.show()


def plot_kraken_results_stacked_bar(df):
    "plot Kraken results as a stacked bar chart"
    df = df[df["Percentage"] > 1]
    df["Taxon"] = df["Taxon"].apply(lambda x: x.split("|")[-1])
    df = df.pivot(index="Reads", columns="Taxon", values="Percentage")
    plt.figure(figsize=(10, 6))
    df.plot(kind="bar", stacked=True, cmap="tab20")
    plt.xlabel("Reads")
    plt.ylabel("Percentage")
    plt.title("Taxonomic Classification (Kraken)")
    plt.tight_layout()
    plt.show()


def plot_kraken_results_bubble_chart(df):
    "plot Kraken results as a bubble chart"
    df = df[df["Percentage"] > 1]
    df["Taxon"] = df["Taxon"].apply(lambda x: x.split("|")[-1])
    df = df.pivot(index="Reads", columns="Taxon", values="Percentage")
    plt.figure(figsize=(10, 6))
    for col in df.columns:
        plt.scatter(df.index, [col] * len(df), s=df[col]
                    * 1000, alpha=0.5, label=col)
    plt.xlabel("Reads")
    plt.ylabel("Taxon")
    plt.title("Taxonomic Classification (Kraken)")
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_metaphlan_results(df):
    "plot Metaphlan results as a bar chart"
    df = df[df["Relative_Abundance"] > 1]
    plt.figure(figsize=(10, 6))
    plt.barh(df["Taxon"], df["Relative_Abundance"], color="lightcoral")
    plt.xlabel("Relative Abundance")
    plt.ylabel("Taxon")
    plt.title("Taxonomic Classification (Metaphlan)")
    plt.tight_layout()
    plt.show()


def plot_metaphlan_results_pie(df):
    "plot Metaphlan results as a pie chart"
    df = df[df["Relative_Abundance"] > 1]
    plt.figure(figsize=(10, 6))
    plt.pie(df["Relative_Abundance"], labels=df["Taxon"], autopct="%1.1f%%")
    plt.title("Taxonomic Classification (Metaphlan)")
    plt.tight_layout()
    plt.show()


def plot_metaphlan_results_heatmap(df):
    "plot Metaphlan results as a heatmap"
    df = df[df["Relative_Abundance"] > 1]
    df["Taxon"] = df["Taxon"].apply(lambda x: x.split("|")[-1])
    df = df.pivot(index="Taxon", columns="Sample", values="Relative_Abundance")
    plt.figure(figsize=(10, 6))
    sns.heatmap(df, cmap="YlGnBu")
    plt.title("Taxonomic Classification (Metaphlan)")
    plt.tight_layout()
    plt.show()


def plot_metaphlan_results_stacked_bar(df):
    "plot Metaphlan results as a stacked bar chart"
    df = df[df["Relative_Abundance"] > 1]
    df["Taxon"] = df["Taxon"].apply(lambda x: x.split("|")[-1])
    df = df.pivot(index="Sample", columns="Taxon", values="Relative_Abundance")
    plt.figure(figsize=(10, 6))
    df.plot(kind="bar", stacked=True, cmap="tab20")
    plt.xlabel("Sample")
    plt.ylabel("Relative Abundance")
    plt.title("Taxonomic Classification (Metaphlan)")
    plt.tight_layout()
    plt.show()


def plot_metaphlan_results_bubble_chart(df):
    "plot Metaphlan results as a bubble chart"
    df = df[df["Relative_Abundance"] > 1]
    df["Taxon"] = df["Taxon"].apply(lambda x: x.split("|")[-1])
    df = df.pivot(index="Sample", columns="Taxon", values="Relative_Abundance")
    plt.figure(figsize=(10, 6))
    for col in df.columns:
        plt.scatter(df.index, [col] * len(df), s=df[col]
                    * 1000, alpha=0.5, label=col)
    plt.xlabel("Sample")
    plt.ylabel("Taxon")
    plt.title("Taxonomic Classification (Metaphlan)")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Load and plot Kraken results
    kraken_df = load_kraken_report("kraken_report.txt")
    plot_kraken_results(kraken_df)
    plot_kraken_results_pie(kraken_df)
    plot_kraken_results_heatmap(kraken_df)
    plot_kraken_results_stacked_bar(kraken_df)
    plot_kraken_results_bubble_chart(kraken_df)

    # Load and plot MetaPhlAn results
    metaphlan_df = load_metaphlan_output("metaphlan_output.txt")
    plot_metaphlan_results(metaphlan_df)
    plot_metaphlan_results_heatmap(metaphlan_df)
    plot_metaphlan_results_stacked_bar(metaphlan_df)
    plot_metaphlan_results_bubble_chart(metaphlan_df)
