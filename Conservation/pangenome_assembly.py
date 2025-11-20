import subprocess
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# --- Configuration ---
# Replace these with your actual genome assembly files (FASTA format)
GENOME_FILES = [
    "assembly1.fasta",
    "assembly2.fasta",
    "assembly3.fasta",
    # Add more assemblies as needed
]
REFERENCE_GENOME = "reference.fasta"
PANGENOME_GRAPH_OUTPUT = "pangenome_graph.gfa"
SV_OUTPUT_TABLE = "structural_variants.tsv"
MINIGRAPH_PATH = "minigraph" # Assumes minigraph is in your PATH

def run_minigraph_build(ref_genome, assemblies, output_gfa):
    """
    Step 1: Build a Pangenome Graph using Minigraph.
    
    This function executes the Minigraph command to build a graph from multiple
    assemblies, anchored on a reference genome.
    """
    print(f"--- 1. Building Pangenome Graph: {output_gfa} ---")
    
    # Construct the Minigraph command: minigraph <ref.fa> <asm1.fa> <asm2.fa> ... > graph.gfa
    minigraph_command = [
        MINIGRAPH_PATH,
        ref_genome
    ] + assemblies
    
    print(f"Executing: {' '.join(minigraph_command)}")
    
    try:
        # Execute the command and direct stdout to the output GFA file
        with open(output_gfa, "w") as outfile:
            subprocess.run(
                minigraph_command,
                check=True,  # Raise an exception for non-zero return codes
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True
            )
        print("Minigraph execution complete. Graph saved.")
        return True
    except FileNotFoundError:
        print(f"Error: Minigraph not found. Please ensure '{MINIGRAPH_PATH}' is in your PATH.")
        print("Installation required: https://github.com/lh3/minigraph")
        return False
    except subprocess.CalledProcessError as e:
        print(f"Error during Minigraph execution. Return code: {e.returncode}")
        print(f"Stderr:\n{e.stderr}")
        return False

def parse_gfa_for_sv_counts(gfa_file, output_tsv):
    """
    Step 2: Parse the GFA file to extract and count structural variations (SVs).
    
    In a Minigraph-generated GFA, structural variations often manifest as
    alternative paths or complex sequences (paths). We'll perform a simple
    counting based on the 'S' (Sequence/Node) records.
    
    NOTE: A proper SV analysis requires mapping reads back to the graph or
    using specialized graph-based variant callers, which is complex. This is 
    a simplified illustrative example.
    """
    print(f"--- 2. Parsing GFA for simple Node/Segment statistics ---")
    node_lengths = []
    
    try:
        with open(gfa_file, 'r') as f:
            for line in f:
                if line.startswith('S'):
                    parts = line.strip().split('\t')
                    # Format: S <segment_id> <sequence> ...
                    sequence = parts[2]
                    node_lengths.append(len(sequence))
        
        # Create a DataFrame for visualization
        df = pd.DataFrame(node_lengths, columns=['Node_Length'])
        df['Node_Type'] = df['Node_Length'].apply(lambda x: 'Small Variant' if x < 50 else 'Core Sequence')
        
        # Save simple statistics (optional)
        df.to_csv(output_tsv, sep='\t', index=False)
        print(f"Extracted {len(df)} nodes. Statistics saved to {output_tsv}")
        return df
    except FileNotFoundError:
        print(f"Error: GFA file not found at {gfa_file}")
        return None

def visualize_sv_distribution(df):
    """
    Step 3: Visualize the distribution of SV/Node sizes using Seaborn.
    
    This function creates a plot to show the frequency of different node/segment
    lengths in the pangenome graph.
    """
    if df is None or df.empty:
        print("No data to visualize.")
        return

    print("--- 3. Visualizing Node Length Distribution with Seaborn ---")
    plt.figure(figsize=(10, 6))
    
    # Use log scale for length to better visualize small vs. large nodes
    # Use a histogram or a K.D.E. plot for distribution
    sns.histplot(
        data=df,
        x='Node_Length',
        log_scale=True,
        bins=50,
        kde=True,
        color='darkblue'
    )
    
    plt.title('Distribution of Pangenome Graph Node Lengths (Log Scale)')
    plt.xlabel('Node Length (bp, Log Scale)')
    plt.ylabel('Frequency (Count)')
    plt.grid(axis='y', alpha=0.5)
    plt.savefig("node_length_distribution.png")
    print("Visualization saved as 'node_length_distribution.png'")
    # 


def main():
    if not os.path.exists(REFERENCE_GENOME):
        print(f"Error: Reference genome file '{REFERENCE_GENOME}' not found.")
        print("Please ensure your data files are correctly named.")
        return

    # 1. Build the Pangenome Graph
    success = run_minigraph_build(REFERENCE_GENOME, GENOME_FILES, PANGENOME_GRAPH_OUTPUT)
    
    if success:
        # 2. Parse the GFA output for statistics
        stats_df = parse_gfa_for_sv_counts(PANGENOME_GRAPH_OUTPUT, SV_OUTPUT_TABLE)
        
        # 3. Visualize the results
        visualize_sv_distribution(stats_df)
    
    print("Pangenome analysis workflow finished.")

if __name__ == "__main__":
    main()