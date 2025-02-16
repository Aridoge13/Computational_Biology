from Bio import SeqIO

def gc_content(seq):
    """Calculate GC content of a DNA sequence."""
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0

def process_fasta(file_path):
    """Read a FASTA file and compute GC content for each sequence."""
    for record in SeqIO.parse(file_path, "fasta"):
        gc = gc_content(record.seq)
        print(f">{record.id}\nGC Content: {gc:.2f}%")

print ("GC content of DNA sequence")

if __name__ == "__main__":
    fasta_file = input("Enter the path to the FASTA file:")
    process_fasta(fasta_file)

print("Thank you for using the GC content calculator!")