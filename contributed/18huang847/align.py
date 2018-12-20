#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''A program for calculating global/local/suboptimal alignment.
Input 2 sequences(DNA and protein) and  substitution_matrix,
use Needleman-Wunsch Smith-Waterman and Waterman-Eggert algorithm.
Convert DNA into protein, align the two protein sequnce by Amino acid.
Scoring by substitution matrix, which include scores with each two Amimo acid,
and gap penalty, which differ in global/local/suboptimal alignment problems.
Backtrace the route of score,and finally,
output scoring table && alignment result shown by conventional pairs.

Author: Bo Huang Intoduction to Bioinformatics, Biol47800, Fall 2018
'''

import os
import sys
import numpy as np

# examine the input file number
# if input arg number is invalid, give correct input format and exit with -1
if len(sys.argv) != 4:
    print("Usage: python bio1_tmp.py DNA_fasta_file protein_fasta_file substitution_matrix_file")
    sys.exit(-1)

#global variables
AA_index = [] # list all AA with order in substitution_matrix order
gap = 0 # gap value from substitution_matrix
N = 10 # maximun n value in suboptimal alignment

# reform fasta file format sequence into a line and store into a string
def fasta2seq(fasta):
    seq = ""
    for line in fasta:
        # fasta file starts with ">" line, ignore the line, read sequence only
        if line.startswith('>'):
            continue
        else:
            seq += line.strip() # sequence in fasta file may be stored in many line, put them together
    return seq

# translate dna sequence into protein sequence
def DNA2protein(seq):
    protein = ""
    table = { # table of Nucleic acid and Amino acid
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    for i in range(0, len(seq), 3): # 3 Nucleic acid correspond to 1 Amino acid
        protein += table[seq[i:i+3]]
    return protein

# tranform the substituetion matric file to list and get gap value
def getProteinSubstituetionMatrix(lines):
    global AA_index
    global gap
    substitution_matrix = []
    for i, line in enumerate(lines):
        linesep = line.split()
        # get AA index from the 1st line
        if 0 == i:
            for j in range(23): # total 23 Amino acid
                AA_index.append(linesep[j])
        # get gap value from last line
        elif 24 == i:
            gap = int(linesep[1])
        else:
            tmp = []
            # get score for AA pairs
            for j in range(1,24):
                tmp.append(int(linesep[j]))
            substitution_matrix.append(tmp)
    return substitution_matrix

# read data from file
def prepareData():
    # read DNA fasta file
    with open(sys.argv[1]) as fin:
        DNA_fasta = fin.readlines()
    # read protein fasta file
    with open(sys.argv[2]) as fin:
        protein_fasta = fin.readlines()
    # read Substituetion Matrix file
    with open(sys.argv[3]) as fin:
        substitution_matrix_file = fin.readlines()
    return fasta2seq(protein_fasta), DNA2protein(fasta2seq(DNA_fasta)), getProteinSubstituetionMatrix(substitution_matrix_file)

# global alignment
def globalAlignment(seq1, seq2, mat):
    print("Global Alignment:")
    global AA_index
    global gap
    best_score = np.zeros((len(seq1), len(seq2)), dtype = np.int) # scoring tabel initalization
    source  = np.zeros((len(seq1), len(seq2)), dtype = np.int) # direction of score source each pair
    # Traversing seq1
    for i, aa1 in enumerate(seq1):
        # Traversing seq2
        for j, aa2 in enumerate(seq2):
            matscore = mat[AA_index.index(aa1)][AA_index.index(aa2)] # read score for aa1 and aa2 from substitution matrix
            # calculate the score from top
            if 0 == i: # head of seq1
                score1 = (j + 1) * gap
            else:
                score1 = best_score[i-1][j] + gap
            # calculate the score from left
            if 0 == j: # head of seq2
                score2 = (i + 1) * gap
            else:
                score2 = best_score[i][j-1] + gap
            # calculate the score from substitution matrix
            if 0 == i: # head of seq1
                score4 = matscore + j * gap
            elif 0 == j: # head of seq2
                score4 = matscore + i * gap
            else: # top left + substitution matrix score
                score4 = matscore + best_score[i-1][j-1]
            best_score[i][j] = max(score1, score2, score4) # get best score from score1 score2 score3
            # record the score source
            if score4 == best_score[i][j]:
                source[i][j] = 4
            elif score2 == best_score[i][j]:
                source[i][j] = 2
            else:
                source[i][j] = 1
    # deal with edge source
    for i in range(len(seq1)):
        source[i][0] = 1
    for j in range(len(seq2)):
        source[0][j] = 2
    print(best_score)
    i = len(seq1) - 1
    j = len(seq2) - 1
    route = [] # store route of alignment
    # backtrace
    while i > 0 or j > 0:
        route.append(source[i][j])
        if source[i][j] == 4:
            i -= 1
            j -= 1
        elif source[i][j] == 2:
            j -= 1
        else:
            i -= 1
    # translate backtrace result into conventional pairs
    result = [] # store final result
    i = len(seq1) - 1
    j = len(seq2) - 1
    for direction in route:
        if direction == 2: # left
            result.append(('.', seq2[j]))
            j -= 1
        elif direction == 1: # top
            result.append((seq1[i], '.'))
            i -= 1
        else: # top left
            result.append((seq1[i], seq2[j]))
            i -= 1
            j -= 1
    # result above does not include the first pair, add it!
    if i == 0 and j == 0:
        result.append((seq1[i], seq2[j]))
    elif i == 0 and j != 0:
        result.append(('.', seq2[j]))
    else:
        result.append((seq1[i], '.'))
    result.reverse() # reverse the result to sequence's order
    # print result at command line
    for i in range(len(result)):
        print(str(result[i][0]), end='')
    print()
    for i in range(len(result)):
        print(str(result[i][1]), end='')
    print()
    return

# local alignment
def localAlignment(seq1, seq2, mat):
    print("Local Alignment:")
    global AA_index
    global gap
    best_score = np.zeros((len(seq1), len(seq2)), dtype = np.int) # scoring tabel initalization
    source  = np.zeros((len(seq1), len(seq2)), dtype = np.int) # direction of score source each pair
    # Traversing seq1
    for i, aa1 in enumerate(seq1):
        # Traversing seq2
        for j, aa2 in enumerate(seq2):
            matscore = mat[AA_index.index(aa1)][AA_index.index(aa2)] # read score for aa1 and aa2 from substitution matrix
            # calculate the score from top
            if 0 == i:
                score1 = gap
            else:
                score1 = best_score[i-1][j] + gap
            # calculate the score from left
            if 0 == j:
                score2 = gap
            else:
                score2 = best_score[i][j-1] + gap
            # calculate the score from substitution matrix
            if 0 == i: # head of seq1
                score4 = matscore + j * gap
            elif 0 == j: # head of seq2
                score4 = matscore + i * gap
            else: # top left + substitution matrix score
                score4 = matscore + best_score[i-1][j-1]
            best_score[i][j] = max(score1, score2, score4, 0) # get best score from score1 score2 score3 and 0
            #record the score source
            if score4 == best_score[i][j]:
                source[i][j] = 4
            elif score2 == best_score[i][j]:
                source[i][j] = 2
            elif score1 == best_score[i][j]:
                source[i][j] = 1
            else: # max = 0
                source[i][j] = 0
    # deal with edge source
    for i in range(len(seq1)):
        source[i][0] = 0
    for j in range(len(seq2)):
        source[0][j] = 0
    print(best_score)
    # backtrace
    highest_locs = np.where(best_score == np.max(best_score)) # get position of highest score
    i, j = highest_locs[0][0], highest_locs[1][0]
    route = [] # store route of alignment
    while i >= 0 and j >= 0 and source[i][j] > 0: # if source == 0, stop backtrace
        route.append(source[i][j])
        if source[i][j] == 4:
            i -= 1
            j -= 1
        elif source[i][j] == 2:
            j -= 1
        else:
            i -= 1
    # translate backtrace result into conventional pairs
    result = [] # store final result
    highest_locs = np.where(best_score == np.max(best_score))
    i, j = highest_locs[0][0], highest_locs[1][0]
    for direction in route:
        if direction == 2: # left
            result.append(('.', seq2[j]))
            j -= 1
        elif direction == 1: # top
            result.append((seq1[i], '.'))
            i -= 1
        elif direction == 4: # top left
            result.append((seq1[i], seq2[j]))
            i -= 1
            j -= 1
    # only 1 pair match doesn't has route, add it!
    if len(route) == 0:
        result.append((seq1[i], seq2[j]))
    result.reverse() # reverse the result to sequence's order
    # print result at command line
    for i in range(len(result)):
        print(str(result[i][0]), end='')
    print()
    for i in range(len(result)):
        print(str(result[i][1]), end='')
    print()
    return

#suboptimal alignment
def suboptimalAlignment(seq1, seq2, mat):
    print("Suboptimal Alignment:")
    global AA_index
    global gap
    global N
    blacklist = [] # store pairs of previous alignment
    # do local alignment for N times
    while N > 0:
        print("%d:" % N)
        best_score = np.zeros((len(seq1), len(seq2)), dtype = np.int) # scoring tabel initalization
        source  = np.zeros((len(seq1), len(seq2)), dtype = np.int) # direction of score source each pair
        # traversing seq1
        for i, aa1 in enumerate(seq1):
            # traversing seq2
            for j, aa2 in enumerate(seq2):
                # if (i, j) is in blacklist(previous alignment), give 0 for score and route, continue to next for loop
                if (i,j) in blacklist:
                    best_score[i][j] = 0
                    source[i][j] = 0
                else:
                    matscore = mat[AA_index.index(aa1)][AA_index.index(aa2)] # read score for aa1 and aa2 from substitution matrix
                    # calculate the score from left
                    if 0 == i:
                        score1 = gap
                    else:
                        score1 = best_score[i-1][j] + gap
                    # calculate the score from top
                    if 0 == j:
                        score2 = gap
                    else:
                        score2 = best_score[i][j-1] + gap
                    # calculate the score from substitution matrix
                    if 0 == i: # head of seq1
                        score4 = matscore + j * gap
                    elif 0 == j: # head of seq2
                        score4 = matscore + i * gap
                    else: # top left + substitution matrix score
                        score4 = matscore + best_score[i-1][j-1]
                    best_score[i][j] = max(score1, score2, score4, 0) # get best score from score1 score2 score3 and 0
                    #record the score source
                    if score4 == best_score[i][j]:
                        source[i][j] = 4
                    elif score2 == best_score[i][j]:
                        source[i][j] = 2
                    elif score1 == best_score[i][j]:
                        source[i][j] = 1
                    else:
                        source[i][j] = 0
        # deal with edge situation
        for i in range(len(seq1)):
            source[i][0] = 0
        for j in range(len(seq2)):
            source[0][j] = 0
        # backtrace
        highest_locs = np.where(best_score == np.max(best_score)) # get position of highest score
        i, j = highest_locs[0][0], highest_locs[1][0]
        # if scoring table is all 0, stop while
        if np.max(best_score) == 0:
            break
        route = [] # store route of alignment
        while i >= 0 and j >= 0 and source[i][j] > 0: # if source == 0, stop backtrace
            route.append(source[i][j])
            blacklist.append((i, j)) # add this pair to blacklist
            if source[i][j] == 4:
                i -= 1
                j -= 1
            elif source[i][j] == 2:
                j -= 1
            elif source[i][j] == 1:
                i -= 1
        # translate backtrace result into conventional pairs
        result = [] # store the final result
        highest_locs = np.where(best_score == np.max(best_score))
        i, j = highest_locs[0][0], highest_locs[1][0]
        for direction in route:
            if direction == 2: # left
                result.append(('.', seq2[j]))
                j -= 1
            elif direction == 1: # top
                result.append((seq1[i], '.'))
                i -= 1
            elif direction == 4: # top left
                result.append((seq1[i], seq2[j]))
                i -= 1
                j -= 1
        # only 1 pair match doesn't has route, add it!
        if len(route) == 0:
            result.append((seq1[i], seq2[j]))
            blacklist.append((i, j))
        result.reverse() # reverse the result to sequence's order
        # print result at command line
        print(best_score)
        for i in range(len(result)):
            print(str(result[i][0]), end='')
        print()
        for i in range(len(result)):
            print(str(result[i][1]), end='')
        print()
        N -= 1 # remain times - 1
    return

# equal to main fuction
if __name__ == "__main__":
    # prepare data from sequences and substitution matrix
    seq1, seq2, mat = prepareData()
    # calculate global alignment
    globalAlignment(seq2, seq1, mat)
    print("--------------------------------------------------")
    # calculate local alignment
    localAlignment(seq2, seq1, mat)
    print("--------------------------------------------------")
    # calculate suboptimal alignment
    suboptimalAlignment(seq2, seq1, mat)