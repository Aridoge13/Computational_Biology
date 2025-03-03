from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

def analyze_protein(seq):
    """Computes molecular weight, isoelectric point, and hydrophobicity."""
    protein = ProteinAnalysis(str(seq))
    properties = {
        "Molecular Weight": protein.molecular_weight(),
        "Isoelectric Point (pI)": protein.isoelectric_point(),
        "Hydrophobicity": protein.gravy()
    }
    return properties

def process_fasta(fasta_file):
    """Reads a FASTA file and processes each protein sequence."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        props = analyze_protein(record.seq)
        print(f"\n> {record.id} - Protein Properties:")
        for prop, value in props.items():
            print(f"{prop}: {value:.2f}")

def main():
    fasta_file = input("Enter the path to the FASTA file: ")
    process_fasta(fasta_file)

if __name__ == "__main__":
    main()
