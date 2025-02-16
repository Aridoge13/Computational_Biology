from Bio import SeqIO, Align
from Bio.Align import PairwiseAligner

def gc_content(seq):
    """Calculate GC content of a DNA sequence."""
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0

def process_fasta(file_path):
    """Read a FASTA file and compute GC content for each sequence."""
    sequences = list(SeqIO.parse(file_path, "fasta"))
    for record in sequences:
        gc = gc_content(record.seq)
        print(f">{record.id}\nGC Content: {gc:.2f}%")
    return sequences

def align_sequences(seq1, seq2):
    """Perform pairwise sequence alignment."""
    aligner = PairwiseAligner()
    aligner.mode = "local"
    alignments = aligner.align(seq1, seq2)
    # print(f"Number of alignments: {len(alignments)}") 
    for alignment in alignments:
        print(alignment)
        break
    
if __name__ == "__main__":
    fasta_files = input("Enter the paths to the FASTA files separated by a space:").split()
    sequences = [seq for file in fasta_files for seq in process_fasta(file)]
    
    if len(sequences) >= 2:
        align_sequences(sequences[0].seq, sequences[1].seq)
