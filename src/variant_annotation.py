import pysam
from Bio import SeqIO
import argparse

# Step 1: Load the reference genome
def load_reference_genome(ref_fasta):
    """Load the reference genome from a FASTA file."""
    reference = SeqIO.read(ref_fasta, "fasta")
    return str(reference.seq)

# Step 2: Identify variants
def call_variants(bam_file, ref_genome, output_vcf):
    """Call variants (SNPs and indels) from a BAM file."""
    bam = pysam.AlignmentFile(bam_file, "rb")
    variants = []

    for pileupcolumn in bam.pileup():
        pos = pileupcolumn.pos  # 0-based position
        ref_base = ref_genome[pos].upper()

        bases = {}
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                read_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                bases[read_base] = bases.get(read_base, 0) + 1

        for base, count in bases.items():
            if base != ref_base and count > 5:  # Simple threshold for variant calling
                variants.append((pos + 1, ref_base, base))  # 1-based position for VCF

    # Write variants to a VCF file
    with open(output_vcf, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for pos, ref, alt in variants:
            vcf.write(f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

    return variants

# Step 3: Annotate variants
def annotate_variants(variants, ref_genome):
    """Annotate variants using Biopython."""
    annotations = []
    for pos, ref, alt in variants:
        codon_start = (pos - 1) // 3 * 3  # Find the start of the codon
        codon = ref_genome[codon_start:codon_start + 3]
        new_codon = list(codon)
        new_codon[(pos - 1) % 3] = alt
        new_codon = "".join(new_codon)

        # Translate codons to amino acids
        from Bio.Seq import Seq
        ref_aa = str(Seq(codon).translate())
        new_aa = str(Seq(new_codon).translate())

        # Determine mutation type
        if ref_aa == new_aa:
            mutation_type = "synonymous"
        elif new_aa == "*":
            mutation_type = "nonsense"
        else:
            mutation_type = "missense"

        annotations.append((pos, ref, alt, ref_aa, new_aa, mutation_type))

    return annotations

# Main function
def main(bam_file, ref_fasta, output_vcf):
    # Load reference genome
    ref_genome = load_reference_genome(ref_fasta)

    # Call variants
    variants = call_variants(bam_file, ref_genome, output_vcf)

    # Annotate variants
    annotations = annotate_variants(variants, ref_genome)

    # Print annotated variants
    print("POS\tREF\tALT\tREF_AA\tALT_AA\tMUTATION_TYPE")
    for ann in annotations:
        print("\t".join(map(str, ann)))

# Parse command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Variant Calling and Annotation")
    parser.add_argument("--bam", required=True, help="Path to the input BAM file")
    parser.add_argument("--ref", required=True, help="Path to the reference genome FASTA file")
    parser.add_argument("--output", required=True, help="Path to the output VCF file")
    args = parser.parse_args()

    # Run the pipeline
    main(args.bam, args.ref, args.output)