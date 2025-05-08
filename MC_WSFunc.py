#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 22:43:04 2025

@author: mattc
"""

import numpy as np


def scale(x, N=2):
    if N <=1:
        return(x)
    result = []
    for i in x:
        result.append(i**N)
    total = sum(result)
    result = [i/total for i in result]
    return(result)


def custom_score(a, b, match=1, mismatch=-1):
    if a == b:
        return match
    if 'N' in (a, b):
        return 0
    if (a in 'CU' and b == 'Y') or (b in 'CU' and a == 'Y'):
        return 0.25
    if (a in 'AG' and b == 'R') or (b in 'AG' and a == 'R'):
        return 0.25
    return mismatch

def SWA(s1, s2, gap_open=-2, gap_extend=-1):
    m, n = len(s1), len(s2)

    # Initialize matrices
    M = np.zeros((m+1, n+1))
    Ix = np.zeros((m+1, n+1))
    Iy = np.zeros((m+1, n+1))

    max_score = 0
    max_pos = None

    for i in range(1, m+1):
        for j in range(1, n+1):
            score = custom_score(s1[i-1], s2[j-1])

            M[i, j] = max(0, M[i-1, j-1], Ix[i-1, j-1], Iy[i-1, j-1]) + score
            Ix[i, j] = max(0, M[i-1, j] + gap_open + gap_extend, Ix[i-1, j] + gap_extend)
            Iy[i, j] = max(0, M[i, j-1] + gap_open + gap_extend, Iy[i, j-1] + gap_extend)

            local_max = max(M[i, j], Ix[i, j], Iy[i, j])
            if local_max >= max_score:
                max_score = local_max
                max_pos = (i, j)

    # Traceback
    align1, align2 = "", ""
    i, j = max_pos
    matrix = np.argmax([M[i, j], Ix[i, j], Iy[i, j]])

    while i > 0 and j > 0:
        if matrix == 0:
            if M[i, j] == 0:
                break
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            score = custom_score(s1[i-1], s2[j-1])
            if M[i, j] == M[i-1, j-1] + score:
                matrix = 0
            elif M[i, j] == Ix[i-1, j-1] + score:
                matrix = 1
            else:
                matrix = 2
            i -= 1
            j -= 1
        elif matrix == 1:
            if Ix[i, j] == 0:
                break
            align1 = s1[i-1] + align1
            align2 = "-" + align2
            if Ix[i, j] == Ix[i-1, j] + gap_extend:
                matrix = 1
            else:
                matrix = 0
            i -= 1
        else:
            if Iy[i, j] == 0:
                break
            align1 = "-" + align1
            align2 = s2[j-1] + align2
            if Iy[i, j] == Iy[i, j-1] + gap_extend:
                matrix = 2
            else:
                matrix = 0
            j -= 1

    return max_score, align1, align2, i, j

import random

def random_nucleotide_sequence(n):
    return ''.join(random.choices('ACGT', k=n))

def fix_gaps(s1,s2):
    seq = ''
    for x in range(len(s2)):
        if s2[x] == '-':
            continue
        else:
            seq = seq + s1[x]
    return(seq)

def consensus_match_ratios(sequences, consensus):
    """
    Given a list of strings (`sequences`) and a consensus string (`consensus`),
    returns a list of ratios indicating the fraction of sequences matching the
    consensus at each position. Also prints a list of indices sorted by most
    mismatched to least mismatched.

    Parameters:
        sequences (List[str]): List of strings, all of equal length.
        consensus (str): A string of the same length as the sequences.

    Returns:
        List[float]: List of ratios of matches per position.
    """
    if not sequences or not consensus:
        raise ValueError("Sequences and consensus cannot be empty.")

    length = len(consensus)
    if any(len(seq) != length for seq in sequences):
        raise ValueError("All sequences and the consensus must be the same length.")

    ratios = []
    mismatches = []
    
    n_seq = len(sequences)
    for i in range(length):
        if consensus[i] == 'N':
            ratios.append(1)
            mismatches.append((i, 0))
            continue
        elif consensus[i] == 'Y':
            matches = sum(seq[i] == 'C' for seq in sequences) + sum(seq[i] == 'T' for seq in sequences) + sum(seq[i] == 'Y' for seq in sequences)
        elif consensus[i] == 'R':
            matches = sum(seq[i] == 'A' for seq in sequences) + sum(seq[i] == 'G' for seq in sequences) + sum(seq[i] == 'R' for seq in sequences)
        else:
            matches = sum(seq[i] == consensus[i] for seq in sequences)
        ratio = matches/n_seq
        ratios.append(ratio)
        mismatches.append((i, 1 - ratio))

    # Sort indices by mismatch rate descending
    sorted_indices_by_mis = [index for index, _ in sorted(mismatches, key=lambda x: x[1], reverse=True)]
    print("Indices sorted by mismatch frequency (most to least):", sorted_indices_by_mis)

    return(ratios,sorted_indices_by_mis)

def find_all_indices(text, char):
    indices = []
    for index, letter in enumerate(text):
        if letter == char:
            indices.append(index)
    return indices

def find_inserts(sequences, consensus): 
    ratios = [0]*len(consensus)
    for x in sequences:
        inds = find_all_indices(x, '-')
        for z in range(len(inds)):
            inds[z] = inds[z] - z
        for z in inds:
            if (0 <z<len(ratios)):
                ratios[z] = ratios[z] +1
    ratios = [y/len(sequences) for y in ratios]
    import matplotlib.pyplot as plt
    x_labels = list(consensus)
    x_positions = range(len(consensus))

    plt.figure(figsize=(10, 4))
    plt.bar(x_positions, ratios)
    plt.xticks(x_positions, x_labels)
    plt.ylim(0, 1)
    plt.xlabel("Consensus Sequence Position")
    plt.ylabel("Insertion Ratio")
    plt.title("Insertion Ratio at Each Position Compared to Consensus")
    plt.tight_layout()
    plt.show()
    return(ratios)

def make_inserts(sequences, ratios):
    probs = [x/sum(ratios) for x in ratios]
    last = probs[0]
    for i in range(1, len(probs)):
        if probs[i] > 0:
            probs[i] = probs[i] +last
            last = probs[i]
    seqs2 = []        
    print(probs)
    for x in sequences:
        n = np.random.random()
        for idx, i in enumerate(probs):
            if n <= i:
                seqs2.append(x[:idx] +'A' + x[idx:])
                seqs2.append(x[:idx] +'C' + x[idx:])
                seqs2.append(x[:idx] +'G' + x[idx:])
                seqs2.append(x[:idx] +'U' + x[idx:])
                break
    seqs2 = {'sequence':seqs2, 'mismatches':[-1]*len(seqs2)}
    import pandas as pd
    return(pd.DataFrame(seqs2))   
    
def plot_match_ratios(consensus, ratios, text=''):
    """
    Plots a bar graph of match ratios with consensus characters as x-axis labels.

    Parameters:
        consensus (str): The consensus sequence.
        ratios (List[float]): List of match ratios at each position.
    """
    import matplotlib.pyplot as plt
    x_labels = list(consensus)
    x_positions = range(len(consensus))

    plt.figure(figsize=(10, 4))
    plt.bar(x_positions, ratios)
    plt.xticks(x_positions, x_labels)
    plt.ylim(0, 1)
    plt.xlabel("Consensus Sequence Position")
    plt.ylabel("Match Ratio")
    plt.title("Match Ratio at Each Position Compared to Consensus for" + text)
    plt.tight_layout()
    plt.show()


def recurs(tree, score_tree, score_max, current_score = 0, current_seq ='', max_len = 10):
    sequences = []
    seq_scores = []
    for x in range(4):
        if current_score + score_tree[0][x] <= score_max:
           sequences.append(current_seq +  tree[0][x])
           seq_scores.append(current_score + score_tree[0][x])
    if (len(tree) == 1) or (len(current_seq)==max_len-1):
        return(sequences,seq_scores)
    else:
        final_sequences = []
        final_scores = []
        for x in range(len(sequences)):
            tempseq, tempscores = recurs(tree[1::], score_tree[1::], score_max, seq_scores[x], sequences[x])
            final_sequences = final_sequences + list(tempseq)
            final_scores = final_scores + list(tempscores)
        return(final_sequences, final_scores)
           
                
#if I dont choose indices, it is (69 choose N/2)(4^(N/2))


def generate_sequences(consensus, indices, score_max=10, maxlen = 10):
    """
    generates fake sequences from the consensus, within a Smith Waterman score of score_max
    It substitutes out most frequently mutated SNPs first, going in order of the list of indices
    maxlen is the maximum number of nucleotides to substitute

    """
    tree = [['A', 'G','C', 'U']]*len(indices)
    score_tree = []
    for x in indices:
        if consensus[x] == 'Y':
            score_tree.append([2,2,.75,.75])
        elif consensus[x] == 'R':
            score_tree.append([.75, .75, 2,2])
        elif consensus[x] not in ['A', 'G','C', 'U']:
            score_tree.append([1,1,1,1])
        elif consensus[x] == 'A':
            score_tree.append([0,2,2,2])
        elif consensus[x] == 'G':
            score_tree.append([2,0,2,2])
        elif consensus[x] == 'C':
            score_tree.append([2,2,0,2])
        elif consensus[x] == 'U':
            score_tree.append([2,2,2,0])
    seqs, scores = recurs(tree, score_tree, score_max, 0, '', maxlen)
    data = {"sequence": seqs, "scores":scores}
    import pandas as pd
    data = pd.DataFrame(data)
    del seqs, scores
    for x in range(len(data['sequence'])):
        temp = list(consensus)
        seq= data['sequence'][x]
        for y in range(len(seq)):
            temp[indices[y]] = seq[y]
        data['sequence'][x] = ''.join(temp)
    return(data)

def expand_possible(frame):
    result = {"sequence":[], 'scores':[]}
    sequences = list(frame['sequence'])
    scores = list(frame['scores'])
    del frame
    for x in range(len(sequences)):
        replacements = sequences[x].replace('A','')
        replacements = replacements.replace('C','')
        replacements = replacements.replace('G','')
        replacements = replacements.replace('U','')
        if len(replacements)==0:
            continue
        temp = [sequences[x]]
        for y in replacements:
            temp2 = []
            for z in temp:
                if y == 'N':
                    temp2.append(z.replace('N', 'A',1))
                    temp2.append(z.replace('N', 'C',1))
                    temp2.append(z.replace('N', 'G',1))
                    temp2.append(z.replace('N', 'U',1))
                if y == 'R':
                    temp2.append(z.replace('R', 'A',1))
                    temp2.append(z.replace('R', 'G',1))
                if y == 'Y':
                    temp2.append(z.replace('Y', 'C',1))
                    temp2.append(z.replace('Y', 'U',1))
            temp = temp2
        result['sequence'] = result['sequence'] + temp
        result['scores'] = result['scores'] + [scores[x]]*len(temp)
    import pandas as pd
    result = pd.DataFrame(result)
    return(result)

def chunk_filter(seq_list, filepath):
    first_chunk = True
    import pandas as pd
    for chunk in pd.read_csv(filepath, chunksize=20000):
        frame = expand_possible(chunk)
        l1 = len(frame)
        frame = frame[~frame["sequence"].isin(seq_list)]
        print('filtered ' + str(l1-len(frame)) + ' sequences')
        del l1
        frame.to_csv(filepath[:-4] + '_filtered.csv', mode='a', index=False, header=first_chunk)
        first_chunk = False

class HMM:
    def __init__(self, ratios, consensus, scaling = 1):
        if len(ratios) != len(consensus):
            print('faulty input, consensus and match ratios are not same size')
            return(False)
        self.transition=1
        self.consensus = consensus
        self.start_score = len(consensus) - consensus.count('N') - consensus.count('Y') - consensus.count('R')
        self.obs_matrix = [[.25,.5,.75, 1]]*len(ratios)
        for x in range(len(self.obs_matrix)):
            if consensus[x] == 'N':
                self.obs_matrix[x] = [.25,.5,.75, 1]
                continue
            elif consensus[x] == 'R':
                self.obs_matrix[x][0] = ratios[x]/2
                self.obs_matrix[x][2] = ratios[x]/2
                self.obs_matrix[x][1] = (1-ratios[x])/2
            elif consensus[x] == 'Y':
                self.obs_matrix[x][1] = ratios[x]/2
                self.obs_matrix[x][0] = (1-ratios[x])/2
                self.obs_matrix[x][2] = (1-ratios[x])/2
            else:
                i = 'ACGU'.index(consensus[x])
                self.obs_matrix[x][i] = ratios[x]
                for j in range(4):
                    if j==i:
                        continue
                    self.obs_matrix[x][j] = (1-ratios[x])/3
            if scaling > 1:
                self.obs_matrix[x] = scale(self.obs_matrix[x], scaling)
            self.obs_matrix[x] = [self.obs_matrix[x][0], sum(self.obs_matrix[x][0:2]), sum(self.obs_matrix[x][0:3]), 1]
            
    def generate(self):
        result = ''
        s = 'ACGU'
        mismatches = 0
        for x in range(len(self.obs_matrix)):
            n = np.random.random()
            for j in range(4):
                if n <= self.obs_matrix[x][j]:
                    result = result + s[j]
                    if s[j] != self.consensus[x]:
                        mismatches = mismatches - custom_score(s[j], self.consensus[x])
                        #print(custom_score(s[j], self.consensus[x]))
                    break
        return(result, mismatches)  
    
    def generate_n(self,n=1):
        results =[]
        scores = []
        for z in range(n):
            result, score = self.generate()
            results.append(result)
            scores.append(score)
        import pandas as pd
        return(pd.DataFrame({'sequence':results, 'mismatches':scores}))
    
    def __str__(self):
        print(self.consensus)
        print(self.obs_matrix)
    

if __name__ == "__main__":
    # Example usage:
    seq = random_nucleotide_sequence(10)
    print(seq)
    
    
    # Example usage
    seq1 = "UUUUGAAUCCUGGCUCAGGACGAACGCUGGCGGCGUGCUUAACACAUGCAAGUCGAACGAUGAAGGCCGUGCUUGCACGGCUGGAUUAGUGGCGAACGGGUGAGUAACACGUGAGUAACCUGCCCUUCACUCUGGGAUAACCUCGGGAAAUCGGGGCUAAUACCGGAUAUGAGUGUCCACUGCAUGGUGGAUGCUGGAAAGUUUUUCGGUGAAGGAUGGACUCGCGGCCUAUCAGUUUGUUGGUGAGGUAAUGGCUCACCAAGACGAUGACGGGUAGCCGGCCUGAGAGGGCGACCGGCCACACUGGGACUGAGACACGGCCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCACAAUGGGCGAAAGCCUGAUGCAGCGACGCCGCGUGAGGGAUGACGGCCUUCGGGUUGUAAACCUCUUUCAGUAGGGAAGAAGCGAAAGUGACGGUACCUGCAGAAGAAGCGCCGGCUAACUACGUGCCAGCAGCCGCGGUAAUACGUAGGGCGCAAGCGUUGUCCGGAAUUAUUGGGCGUAAAGAGCUUGUAGGUGGCUUGUCGCGUCUGCCGUGAAAAUCCAGGGCUUAACUCUGGACGUGCGGUGGGUACGGGCAGGCUAGAGUGUGGUAGGGGAGACUGGAACUCCUGGUGUAGCGGUGAAAUGCGCAGAUAUCAGGAAGAACACCGAUGGCGAAGGCAGGUCUCUGGGCCAUUACUGACACUGAGAAGCGAAAGCAUGGGGAGCGAACAGGAUUAGAUACCCUGGUAGUCCAUGCCGUAAACGUUGGGCACUAGGUGUGGGCGACAUUCCACGUUGUCUGCGCCGUAGCUAACGCAUUAAGUGCCCCGCCUGGGGAGUACGGCCGCAAGGCUAAAACUCAAAGGAAUUGACGGGGGCCCGCACAAGCGGCGGAGCAUGCGGAUUAAUUCGAUGCAACGCGAAGAACCUUACCAAGGCUUGACAUGCACCGGACAGCUGCAGAGAUGUGGCUUUCUUUGGACUGGUGCACAGGUGGUGC"
    seq3 = 'GAGUUUGAUCCUGGCUCAGGACGAACGCUGGCGGCGUGCUUAACACAUGCAAGUCGAACGAUGAAGCCUUUCGGGGUGGAUUAGUGGCGAACGGGUGAGUAACACGUGGGCAAUCUGCCCUUCACUCUGGGACAAGCCCUGGAAACGGGGUCUAAUACCGGAUAAUACUUUCUCCUGCAUGGGGGAGGGUUGAAAGCUCCGGCGGUGAAGGAUGAGCCCGCGGCCUAUCAGCUUGUUGGUGGGGUGAUGGCCUACCAAGGCGACGACGGGUAGCCGGCCUGAGAGGGCGACCGGCCACACUGGGACUGAGACACGGCCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCACAAUGGGCGGAAGCCUGAUGCAGCGACGCCGCGUGAGGGAUGACGGCCUUCGGGUUGUAAACCUCUUUCAGCAGGGAAGAAGCGUGAGUGACGGUACCUGCAGAAGAAGCGCCGGCUAACUACGUGCCAGCAGCCGCGGUAAUACGUAGGGCGCAAGCGUUGUCCGGAAUUAUUGGGCGUAAAGAGCUCGUAGGCGGCUUGUCGCGUCGGAUGUGAAAGCCCGGGGCUUAACCCCGGGUCUGCAUUCGAUACGGGCAGGCUAGAGUGUGGUAGGGGAGAUCGGAAUUCCUGGUGUAGCGGUGAAAUGCGCAGAUAUCAGGAGGAACACCGGUGGCGAAGGCGGAUCUCUGGGCCAUUACUGACGCUGAGGAGCGAAAGCGUGGGGAGCGAACAGGAUUAGAUACCCUGGUAGUCCACGCCGUAAACGUUGGGAACUAGGUGUUGGCGACAUUCCACGUCGUCGGUGCCGCAGCUAACGCAUUAAGUUCCCCGCCUGGGGAGUACGGCCGCAAGGCUAAAACUCAAAGGAAUUGACGGGGGCCCGCACAAGCAGCGGAGCAUGUGGCUUAAUUCGACGCAACGCGAAGAACCUUACCAAGGCUUGACAUACACCGGAAACGGCCAGAGAUGGUCGCCCCCUUGUGGUCGGUGUACAGGUGGUGCAUGGCUGUCGUCAGCUCGUGUCGUGAGAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUUGUUCUGUGUUGCCAGCAUGCCUUUCGGGGUGAUGGGGACUCACAGGAGACUGCCGGGGUCAACUCGGAGGAAGGUGGGGACGACGUCAAGUCAUCAUGCCCCUUAUGUCUUGGGCUGCACACGUGCUACAAUGGCCGGUACAAUGAGCUGCGAUGCCGUGAGGCGGAGCGAAUCUCAAAAAGCCGGUCUCAGUUCGGAUUGGGGUCUGCAACUCGACCCCAUGAAGUCGGAGUUGCUAGUAAUCGCAGAUCAGCAUUGCUGCGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACGUCACGAAAGUCGGUAACACCCGAAGCCGGUGGCCCAACCCCCUUGUGGGGAGGGAGCUGUCGAAGGUGGGACUGGCGAUUGGGACGAAGUCGUAACAAGGUAGCC'
    seq2 = "AGGUGNUGCAUGGYYGYCGUCAGCUCGUGYCGUGAGUGUUGGGUUAAGUCCCRYAACGAGCGCAACCCU"
    seq4 = "AGGUGAUGCAUGGUUGUCGUCAGCUCGUGUCGUGAGUGUUGGGUUAAGUCCCAUAACGAGCGCAACCCU"
    score, aligned_seq1, aligned_seq2, i, j = SWA(seq1, seq2)
    print("Alignment Score:", score)
    print(aligned_seq1)
    print(aligned_seq2)
    
    score2, aligned_seq3, aligned_seq4 = SWA(random_nucleotide_sequence(463), random_nucleotide_sequence(462))
    print("Alignment Score:", score2)
    print(aligned_seq3)
    print(aligned_seq4)
    
    sequences = ["ATCG", "ATGG", "ATCC", "ATCG", "GGGG"]
    consensus = "ATCG"
    ratios = consensus_match_ratios(sequences, consensus)
    plot_match_ratios(consensus, ratios)