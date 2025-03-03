# Computational Biology Beginner Projects

## Description
This project provides a set of Python scripts for analyzing nucleotide sequences, focusing on codon usage, GC content calculation, and sequence alignment. It is designed to assist in bioinformatics research and genomic, transcriptomic and proteomic data analysis.
To access the programs, please open the folder named "src". Thank you. 

## Features
- codon_analysis.py for codon analysis
    - Reads DNA sequences from a FASTA file.
    - Counts codon usage frequency.
    - Visualizes codon distribution using a bar chart.

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

## Installation
Instructions on how to install and set up your project.

```bash
# Clone the repository
git clone https://github.com/Aridoge13/Computational_Biology.git


# Install dependencies
pip install