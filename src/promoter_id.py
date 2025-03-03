from Bio import SeqIO
import re

def extract_promoter_region(fasta_file, upstream_length=1000):
    """Extracts promoter regions (upstream sequences) from a FASTA file."""
    return {record.id: record.seq[:upstream_length] for record in SeqIO.parse(fasta_file, "fasta")}

def find_promoter_motifs(sequence):
    """Finds common promoter motifs like TATA box and GC-rich regions."""
    motifs = {
        "TATA Box": re.search(r"TATA[AT]A[AT]", str(sequence)),
        "GC-Rich": re.search(r"G{3,}C{3,}", str(sequence))
    }
    return {motif: match.group(0) if match else "Not Found" for motif, match in motifs.items()}

def main():
    fasta_file = input("Enter the path to the FASTA file: ")
    promoters = extract_promoter_region(fasta_file)

    for gene_id, seq in promoters.items():
        motifs = find_promoter_motifs(seq)
        print(f"\n> {gene_id} - Promoter Region:")
        print(f"Sequence: {seq[:50]}...")  # Display only the first 50 bases
        print(f"Motifs: {motifs}")

if __name__ == "__main__":
    main()
