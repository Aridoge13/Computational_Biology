import os
import subprocess
import argparse
from Bio import Phylo

# Step 1: Multiple Sequence Alignment with MAFFT
def run_mafft(input_fasta, output_alignment):
    """Run MAFFT to perform multiple sequence alignment."""
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input FASTA file {input_fasta} not found.")
    
    # Run MAFFT
    command = f"mafft --auto {input_fasta} > {output_alignment}"
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Alignment saved to {output_alignment}")
    except subprocess.CalledProcessError as e:
        print(f"Error during MAFFT alignment: {e}")
        raise

# Step 2: Build Phylogenetic Tree with RAxML
def run_raxml(aligned_fasta, output_tree_name, output_dir):
    """Run RAxML to build a phylogenetic tree."""
    if not os.path.exists(aligned_fasta):
        raise FileNotFoundError(f"Aligned FASTA file {aligned_fasta} not found.")
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Run RAxML with the specified output directory
    command = f"raxmlHPC -s {aligned_fasta} -n {output_tree_name} -m GTRGAMMA -p 12345 -w {output_dir}"
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Phylogenetic tree saved to {output_dir}/RAxML_bestTree.{output_tree_name}")
    except subprocess.CalledProcessError as e:
        print(f"Error during RAxML tree building: {e}")
        raise

# Step 3: Visualize Phylogenetic Tree with Biopython
def visualize_tree(tree_file):
    """Visualize the phylogenetic tree using Biopython."""
    if not os.path.exists(tree_file):
        raise FileNotFoundError(f"Tree file {tree_file} not found.")
    
    # Load the tree
    tree = Phylo.read(tree_file, "newick")
    
    # Draw the tree
    Phylo.draw(tree)

# Main function
def main(input_fasta, output_alignment, output_tree_name, output_dir):
    # Step 1: Perform multiple sequence alignment
    print("Running MAFFT for multiple sequence alignment...")
    run_mafft(input_fasta, output_alignment)
    
    # Step 2: Build phylogenetic tree
    print("Running RAxML for phylogenetic tree building...")
    run_raxml(output_alignment, output_tree_name, output_dir)
    
    # Step 3: Visualize the tree
    print("Visualizing the phylogenetic tree...")
    tree_file = f"{output_dir}/RAxML_bestTree.{output_tree_name}"
    visualize_tree(tree_file)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Phylogenetic Analysis Pipeline")
    parser.add_argument("--input", required=True, help="Path to the input FASTA file")
    parser.add_argument("--alignment", required=True, help="Path to save the alignment output")
    parser.add_argument("--tree", required=True, help="Name of the output tree file (without extension)")
    parser.add_argument("--output_dir", required=True, help="Directory to save RAxML output files")
    args = parser.parse_args()

    # Run the pipeline
    main(args.input, args.alignment, args.tree, args.output_dir)