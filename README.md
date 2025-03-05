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
    - Calculate GC content for sequences in a FASTA file and perform pairwise sequence alignment between two sequences.
    - The script reads sequences from one or more FASTA files, computes the GC content for each sequence, and aligns the first two sequences using local pairwise alignment.
    
    - How to use:
        - Install Dependencies:
            - Ensure the following Python library is installed:
                - Biopython: pip install biopython
        
        - Run the Script:
            - Run the script and provide the paths to the FASTA files when prompted: python gc_calc_pwa.py
            - Input: When prompted, enter the paths to the FASTA files separated by a space (e.g., file1.fasta file2.fasta).
        
        - Output: 
            - GC Content: The GC content of each sequence is printed to the terminal.
            - Pairwise Alignment: If at least two sequences are provided, the script performs a local pairwise alignment between the first two sequences and prints the alignment.


- gc_calc.py for calculating GC content of a single DNA sequence
    - Calculate the GC content of DNA sequences in a FASTA file.
    - The script reads sequences from a FASTA file and computes the GC content for each sequence
    
    - How to use:
        - Install Dependencies:
            - Ensure the following Python library is installed:
                - Biopython: pip install biopython
        
        - Run the Script:
            - Run the script and provide the path to the FASTA file when prompted: python gc_calc.py
            - Input: When prompted, enter the path to the FASTA file (e.g., file.fasta).
        
        - Output:
            - GC Content: The GC content of each sequence is printed to the terminal.
            - A friendly message is displayed at the end of the program: "Thank you for using the GC content calculator!"


- protein_seq_motif.py for detecting protein motifs in a given sequence
    - Scan protein sequences in a FASTA file for occurrences of a specific motif.
    - The script uses regular expressions (regex) to identify motifs in protein sequences and reports their positions and sequences.
    
    - How to use:
        - Install Dependencies:
            - Ensure the following Python library is installed:
                - Biopython: pip install biopython
        
        - Run the Script:
            - Run the script and provide the path to the FASTA file and the motif pattern when prompted: python protein_seq_motif.py
            - Input:
                - When prompted, enter the path to the FASTA file (e.g., file.fasta).
                - Enter the motif pattern as a regular expression (e.g., [AG].{2}G for motifs starting with A or G, followed by any two characters, and ending with G).
        - Output: For each sequence in the FASTA file, the script prints:
            - The sequence ID.
            - The number of motifs found.
            - The position and sequence of each motif (if any).
            - If no motifs are found, it reports No motifs found.


- rna_analysis.py for basic RNA-Seq Data Analysis
    - Perform basic RNA-Seq data analysis, including filtering differentially expressed genes and visualizing results.
    - The script reads RNA-Seq data from a CSV file, filters genes based on fold change and p-value thresholds, and generates a bar plot of the top differentially expressed genes.
    
    - How to use:
        - Install Dependencies:
            - Ensure the following Python libraries are installed:
                - Pandas: pip install pandas
                - Matplotlib: pip install matplotlib
        
        - Run the Script:
            - Run the script and provide the path to the CSV file when prompted: python rna_analysis.py 
            - Input:
            - When prompted, enter the path to the CSV file containing RNA-Seq data.
            - Provide the following thresholds (or use defaults):
                - Fold change threshold (default: 2.0).
                - P-value threshold (default: 0.05).
                - Whether to use adjusted p-values (y for yes, n for no).
        
        - Output:
            - Filtered Genes: The script prints the number of genes passing the filters and displays the top results.
            - Volcano Plot: A bar plot of the top differentially expressed genes is generated, with up-regulated genes in red and down-regulated genes in blue.



- seq_align_score_calc.py
    - Perform global and local sequence alignment using the Needleman-Wunsch and Smith-Waterman algorithms.
    - The script reads sequences from FASTA files, aligns them using predefined scoring matrices (BLOSUM62 for global alignment and PAM250 for local alignment), and outputs the alignment and scores.
    
    - How to use:
        - Install Dependencies:
            - Ensure the following Python library is installed:
                - NumPy: pip install numpy
        - Run the Script:
            - Run the script and provide the paths to the FASTA files when prompted: python seq_align_score.py
            - Input: When prompted, enter the paths to the FASTA files separated by a space (e.g., file1.fasta file2.fasta).
        - Output:
            - Global Alignment:
                - Uses the Needleman-Wunsch algorithm with the BLOSUM62 scoring matrix.
                - Prints the aligned sequences and the alignment score.
            - Local Alignment:
                - Uses the Smith-Waterman algorithm with the PAM250 scoring matrix.
                - Prints the maximum local alignment score.



- dna-protein.py for DNA to Protein translation
    - Translate DNA sequences into protein sequences using specified reading frames.
    - The script reads DNA sequences from a FASTA file, translates them into amino acid sequences, and allows the user to choose a specific reading frame or translate in all three forward frames.
    
    - How to use:
        - Install Dependencies:
            - Ensure the following Python library is installed:
                - Biopython: pip install biopython
        - Run the Script:
            - Run the script and provide the path to the FASTA file when prompted: python dna-protein.py
        - Input:
            - When prompted, enter the path to the FASTA file containing DNA sequences.
            - Choose a reading frame (1, 2, or 3) for translation.
            - Optionally, translate the sequence in all three forward frames.
        
        - Output:
            - Protein Sequence: The translated protein sequence for the selected reading frame is displayed.
            - All Frames: If requested, the script displays the translated sequences for all three forward reading frames.



- promoter_id.py
    - Extract promoter regions from a FASTA file and identify common promoter motifs (e.g., TATA box and GC-rich regions).
    - The script reads sequences from a FASTA file, extracts the upstream promoter regions, and searches for motifs in these regions.
    
    - How to use:
        - Install Dependencies:
            - Ensure the following Python library is installed:
                - Biopython: pip install biopython
        - Run the Script:
            - Run the script and provide the path to the FASTA file when prompted: python promoter_id.py
        - Input: When prompted, enter the path to the FASTA file containing DNA sequences.
        - Output:
            - Promoter Regions: For each sequence, the script extracts the upstream promoter region (default: 1000 bases).
            - Motifs: The script identifies and reports the presence of common promoter motifs (TATA box and GC-rich regions).
            - The first 50 bases of the promoter region and the motifs found are printed for each sequence.



- protein_property.py
    - Compute molecular weight, isoelectric point (pI), and hydrophobicity for protein sequences in a FASTA file.
    - The script reads protein sequences from a FASTA file and calculates key physicochemical properties using Biopython's ProteinAnalysis.

    - How to use:
        - Install Dependencies:
            - Ensure the following Python library is installed:
                - Biopython: pip install biopython
        - Run the Script:
            - Run the script and provide the path to the FASTA file when prompted: python protein_property.py
        - Input:
            - When prompted, enter the path to the FASTA file containing protein sequences.
        - Output:
            - Protein Properties: For each sequence, the script calculates and prints:
                - Molecular weight (in Daltons).
                - Isoelectric point (pI).
                - Hydrophobicity (GRAVY score).



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