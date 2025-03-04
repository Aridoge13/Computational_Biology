# Computational Biology Beginner Projects

## Description
This project provides a set of Python scripts for analyzing nucleotide sequences, focusing on codon usage, GC content calculation, and sequence alignment. It is designed to assist in bioinformatics research and genomic, transcriptomic and proteomic data analysis.
To access the programs, please open the folder named "src". Thank you. 

## Features
- codon_analysis.py for codon analysis
    - Analyze codon usage in DNA sequences from a FASTA file and visualize the results as a bar chart.
    - The script reads sequences from a FASTA file, counts the occurrences of each codon, and generates a bar chart showing codon usage frequencies.
    - How to use:
        - Install dependencies:
            - Ensure the following Python libraries are installed:
                - Biopython: pip install biopython
                - Matplotlib: pip install matplotlib

        - Run the Script:
            - Run the script and provide the path to the FASTA file when prompted:
            python codon_analysis.py
            - Input: When prompted, enter the path to the FASTA file containing DNA sequences.
    
    - Output:
        - A bar chart displaying the frequency of each codon in the sequences.
        - The chart is displayed using Matplotlib.

- gc_calc_pwa.py for calculating GC content and global alignment of 2 or more nucleotide sequences
    - Computes GC content of multiple nucleotide sequences.
    - Performs global sequence alignment for two or more DNA sequences.
    - Provides visual representation of GC content distribution.

- gc_calc.py for calculating GC content of a single DNA sequence
    - Calculates the GC content of a single DNA sequence.
    - Outputs GC percentage for sequence analysis.

- protein_seq_motif.py for detecting protein motifs in a given sequence
    - Detects the site of the desired protein motif
    - Detects the number of the desired protein motifs present in the protein sequence

- rna_analysis.py for basic RNA-Seq Data Analysis
    - Basic gene expression studies. 

- seq_align_score_calc.py
    - Implements a simple Needleman-Wunsch or Smith-Waterman algorithm for pairwise sequence alignment.
    - Comparing scoring matrices like PAM and BLOSUM.

- dna-protein.py for DNA to Protein translation
    - Detects the codons in the DNA to predict the primary structure of the protein that is translated from the DNA sequence.
    - Can give the result for 3 reading frames.

- promoter_id.py
    - Extracts the promoter region (upstream sequence) of a gene.
    - Analyzes GC content and possible promoter motifs (like the TATA box).

- protein_property.py
    - Computes molecular weight, isoelectric point (pI), and hydrophobicity for proteins.

- phylogenetic_analysis.py
    - Build phylogenetic trees from multiple sequence alignments and analyze evolutionary relationships.
    - Ensure the input FASTA file contains at least 4 sequences.
    - The output directory (--output_dir) will be created automatically if it doesn't exist.
    - How to use:
        - Install Dependencies:
        - Ensure the following tools are installed:
            - MAFFT: sudo apt install mafft
            - RAxML: sudo apt install raxml
            - Biopython: pip install biopython
        - Run the Script:
        - Use the following command to run the pipeline:
            - python phylogenetic_analysis.py --input <input_fasta> --alignment <output_alignment> --tree <tree_name> --output_dir <output_directory>
        - Arguments:
            - --input: Path to the input FASTA file containing sequences.
            - --alignment: Path to save the aligned sequences (output from MAFFT).
            - --tree: Name of the output tree file (without extension).
            - --output_dir: Directory to save RAxML output files.
        - Output:
            - Alignment File: Saved to the path specified by --alignment.
            - Phylogenetic Tree: Saved to <output_directory>/RAxML_bestTree.<tree_name>.
            - Visualization: The phylogenetic tree is displayed using Biopython.

- variant_calling.py
    - Identify and annotate genetic variants (SNPs and indels) from a BAM file using a reference genome.
    - The script assumes the BAM file is sorted and indexed. If not, use samtools sort and samtools index to prepare the BAM file.
    - The reference genome FASTA file should be in standard FASTA format.
        - How to use:
            - Install Dependencies:
            Ensure the following Python libraries are installed:
            - pysam: pip install pysam
            - Biopython: pip install biopython
    - Run the Script:
        - Use the following command to run the pipeline: python variant_calling.py --bam <input_bam> --ref <reference_fasta> --output <output_vcf>
    - Arguments:
        - --bam: Path to the input BAM file containing aligned reads.
        - --ref: Path to the reference genome FASTA file.
        - --output: Path to save the output VCF file containing called variants.

    - Output:
        - VCF File: Saved to the path specified by --output. Contains the called variants
        - Annotations: Printed to the terminal, showing the variant position, reference allele, alternate allele, reference amino acid, alternate amino acid, and mutation type (synonymous, missense, or nonsense).


## Installation
Instructions on how to install and set up your project.

```bash
# Clone the repository
git clone https://github.com/Aridoge13/Computational_Biology.git


# Install dependencies
pip install