import numpy as np
from itertools import product

# Simplified scoring matrices (subset for demonstration)
BLOSUM62 = {
    ('A', 'A'): 4, ('A', 'R'): -1, ('A', 'N'): -2,
    ('R', 'A'): -1, ('R', 'R'): 5, ('R', 'N'): 0,
    ('N', 'A'): -2, ('N', 'R'): 0, ('N', 'N'): 6
}

PAM250 = {
    ('A', 'A'): 2, ('A', 'R'): -2, ('A', 'N'): 0,
    ('R', 'A'): -2, ('R', 'R'): 6, ('R', 'N'): 0,
    ('N', 'A'): 0, ('N', 'R'): 0, ('N', 'N'): 2
}

def create_matrix(rows, cols):
    return np.zeros((rows, cols))

def needleman_wunsch(seq1, seq2, matrix, gap_penalty=-4):
    m, n = len(seq1), len(seq2)
    score = create_matrix(m+1, n+1)
    
    # Initialize gaps
    for i in range(1, m+1):
        score[i][0] = gap_penalty * i
    for j in range(1, n+1):
        score[0][j] = gap_penalty * j

    # Fill matrix
    for i, j in product(range(1, m+1), range(1, n+1)):
        match = score[i-1][j-1] + matrix.get((seq1[i-1], seq2[j-1]), -3)
        delete = score[i-1][j] + gap_penalty
        insert = score[i][j-1] + gap_penalty
        score[i][j] = max(match, delete, insert)
    
    return score

def smith_waterman(seq1, seq2, matrix, gap_penalty=-4):
    m, n = len(seq1), len(seq2)
    score = create_matrix(m+1, n+1)
    max_score = 0
    
    for i, j in product(range(1, m+1), range(1, n+1)):
        match = max(0, score[i-1][j-1] + matrix.get((seq1[i-1], seq2[j-1]), -3))
        delete = max(0, score[i-1][j] + gap_penalty)
        insert = max(0, score[i][j-1] + gap_penalty)
        score[i][j] = max(match, delete, insert)
        if score[i][j] > max_score:
            max_score = score[i][j]
    
    return score, max_score

def traceback(score, seq1, seq2, matrix, gap_penalty=-4):
    align1, align2 = "", ""
    i, j = len(seq1), len(seq2)
    
    while i > 0 and j > 0:
        current = score[i][j]
        diagonal = score[i-1][j-1]
        up = score[i][j-1]
        left = score[i-1][j]
        
        if current == diagonal + matrix.get((seq1[i-1], seq2[j-1]), -3):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif current == left + gap_penalty:
            align1 = seq1[i-1] + align1
            align2 = f"-{align2}"
            i -= 1
        else:
            align1 = f"-{align1}"
            align2 = f"{seq2[j-1]}{align2}"
            j -= 1
    
    return align1, align2

def process_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ""
        for line in file:
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                    sequence = ""
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences

if __name__ == "__main__":
    fasta_files = input("Enter the paths to the FASTA files separated by a space:").split()
    sequences = [seq for file in fasta_files for seq in process_fasta(file)]

    seq1 = sequences[0]
    seq2 = sequences[1]

    print("Global Alignment (Needleman-Wunsch) with BLOSUM62:")
    score_matrix = needleman_wunsch(seq1, seq2, BLOSUM62)
    a1, a2 = traceback(score_matrix, seq1, seq2, BLOSUM62)
    print(f"Alignment:\n{a1}\n{a2}\nScore: {score_matrix[-1][-1]}\n")

    print("Local Alignment (Smith-Waterman) with PAM250:")
    sw_matrix, max_score = smith_waterman(seq1, seq2, PAM250)
    print(f"Max local score: {max_score}")