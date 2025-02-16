from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import defaultdict

def count_codons(seq):
    """Count occurrences of each codon in a DNA sequence."""
    codon_counts = defaultdict(int)
    seq = seq.upper()
    for i in range(0, len(seq) - 2, 3):  # Iterate in codon (triplet) steps
        codon = seq[i:i+3]
        if len(codon) == 3:
            codon_counts[codon] += 1
    return codon_counts

def process_fasta(file_path):
    """Read sequences from a FASTA file and analyze codon usage."""
    total_codon_counts = defaultdict(int)
    for record in SeqIO.parse(file_path, "fasta"):
        codon_counts = count_codons(str(record.seq))
        for codon, count in codon_counts.items():
            total_codon_counts[codon] += count
    return total_codon_counts

def plot_codon_usage(codon_counts):
    """Plot codon usage as a bar chart."""
    codons, counts = zip(*sorted(codon_counts.items()))  # Sort for consistent order
    plt.figure(figsize=(12, 6))
    plt.bar(codons, counts, color='skyblue', edgecolor='black')
    plt.xlabel("Codons")
    plt.ylabel("Frequency")
    plt.title("Codon Usage Analysis")
    plt.xticks(rotation=90)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

if __name__ == "__main__":
    fasta_file = input("Enter the path to the FASTA file: ")
    codon_usage = process_fasta(fasta_file)
    plot_codon_usage(codon_usage)
