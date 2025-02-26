import pandas as pd
import matplotlib.pyplot as plt

def load_rna_seq_data(file_path):
    """Load RNA-Seq data and validate required columns."""
    REQUIRED_COLUMNS = {'gene', 'log2FoldChange', 'pvalue' , 'padj'}
    try:
        df = pd.read_csv(file_path)
        if missing := REQUIRED_COLUMNS - set(df.columns):
            raise ValueError(f"Missing columns: {missing}")
        return df
    except FileNotFoundError:
        raise FileNotFoundError(f"File {file_path} not found.")
    except pd.errors.EmptyDataError:
        raise ValueError("The CSV file is empty.")

def filter_genes(df, fold_change_thresh=2.0, p_value_thresh=0.05, use_padj=True):
    """Filter genes using fold change and adjusted/raw p-value."""
    p_col = 'padj' if use_padj else 'pvalue'
    filtered = df[
        (df['log2FoldChange'].abs() >= fold_change_thresh) & 
        (df[p_col] <= p_value_thresh)
    ]
    return filtered.sort_values('log2FoldChange', ascending=False)

def plot_volcano(df, top_n=20):
    """Plot top up/down-regulated genes with absolute fold change."""
    if df.empty:
        print("No data to plot.")
        return
    
    # Sort by absolute fold change and take top N
    df = df.iloc[df['log2FoldChange'].abs().argsort()[::-1]].head(top_n)
    
    plt.figure(figsize=(12, 8))
    colors = ['red' if fc > 0 else 'blue' for fc in df['log2FoldChange']]
    plt.barh(df['gene'], df['log2FoldChange'], color=colors)
    plt.xlabel("log2 Fold Change")
    plt.ylabel("Gene")
    plt.title(f"Top {top_n} Differentially Expressed Genes (Red=Up, Blue=Down)")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    file_path = input("Enter CSV path: ")
    df = load_rna_seq_data(file_path)
    
    # Get user-defined thresholds
    fc_thresh = float(input("Fold change threshold (default=2.0): ") or 2.0)
    p_thresh = float(input("P-value threshold (default=0.05): ") or 0.05)
    use_padj = input("Use adjusted p-value? (y/n): ").lower() == 'y'
    
    filtered = filter_genes(df, fc_thresh, p_thresh, use_padj)
    
    if filtered.empty:
        print("No genes passed filters.")
    else:
        print(f"Found {len(filtered)} genes:\n", filtered.head())
        plot_volcano(filtered)