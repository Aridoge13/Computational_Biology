import re
from Bio import SeqIO

def find_motifs(sequence, motif_pattern):
    """Find all occurrences of a motif in a protein sequence."""
    return [(match.start(), match.end()) for match in re.finditer(motif_pattern, str(sequence))]

def process_fasta(file_path, motif_pattern):
    """Scan multiple protein sequences in a FASTA file for motifs."""
    for record in SeqIO.parse(file_path, "fasta"):
        if matches := find_motifs(record.seq, motif_pattern):
            print(f">{record.id} | Motifs found: {len(matches)}")
            for start, end in matches:
                print(f" - Position: {start+1}-{end} | Motif: {record.seq[start:end]}")
        else:
            print(f">{record.id} | No motifs found.")

if __name__ == "__main__":
    fasta_file = input("Enter the path to the FASTA file: ")
    motif_pattern = input("Enter the motif pattern (regex): ")
    process_fasta(fasta_file, motif_pattern)

