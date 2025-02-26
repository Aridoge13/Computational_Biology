from Bio import SeqIO
from Bio.Seq import Seq

def translate_dna(dna_seq, frame=1):
    """Translate a DNA sequence into an amino acid sequence with a specified reading frame.

    Args:
        dna_seq (str): DNA sequence (consisting of A, T, C, G).
        frame (int): Reading frame (1, 2, or 3).

    Returns:
        str: Translated amino acid sequence.
    """
    adjusted_seq = dna_seq[frame - 1:]  # Adjust sequence based on reading frame
    seq_obj = Seq(adjusted_seq)
    protein_seq = seq_obj.translate(to_stop=False)  # Translate entire sequence without stopping at stop codons
    return str(protein_seq)

def translate_all_frames(dna_seq):
    """Translate the DNA sequence in all three forward reading frames."""
    return {frame: translate_dna(dna_seq, frame) for frame in range(1, 4)}

def process_fasta(file_path):
    """Read a FASTA file and translate each sequence."""
    try:
        records = list(SeqIO.parse(file_path, "fasta"))
        if not records:
            print("Error: No sequences found in the provided FASTA file.")
            return

        for record in records:
            print(f"\n> {record.id}")
            frame = int(input("Choose reading frame (1, 2, or 3): ") or 1)
            if frame not in [1, 2, 3]:
                print("Invalid frame selected. Using default frame 1.")
                frame = 1

            protein_seq = translate_dna(str(record.seq), frame)
            print(f"Protein sequence (Frame {frame}): {protein_seq}")

            show_all = input("Translate in all frames? (y/n): ").strip().lower()
            if show_all == "y":
                translations = translate_all_frames(str(record.seq))
                for f, seq in translations.items():
                    print(f"Frame {f}: {seq}")

    except FileNotFoundError:
        print("Error: File not found. Please check the file path.")
    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    print("DNA to Protein Translator")
    fasta_path = input("Enter the path to the FASTA file: ").strip()
    process_fasta(fasta_path)

if __name__ == "__main__":
    main()

