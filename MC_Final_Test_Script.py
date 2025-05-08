#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: mattc
"""

PATH = '/YOUR/PATH/HERE/'

import pandas as pd
import sys
sys.path.append(PATH)
from Bio import SeqIO
from MC_hash_table import DNAHashTable
ht = DNAHashTable(2**25)



preprocessed = pd.read_csv(PATH + 'All_preprocessed_unique.csv')
for x in range(len(preprocessed)): 
    temp = preprocessed['original sequence'][x]
    temp = temp.replace('U', 'T')
    ht.insert_69mer(temp, x)  # store index of a 69-mer
del x, temp


possible = pd.read_csv(PATH + 'possible_sequences.csv')
for x in range(len(possible)): 
    temp = possible['sequence'][x]
    temp = temp.replace('U', 'T')
    ht.insert_69mer(temp, -2)  # store index of a 69-mer
del x, possible, temp

sequences = list(preprocessed['original sequence'])
hits = {'read':[], 'sequence':[], 'species':[]}
c=0
#this was the first one I tried:
#'/home/mattc/Desktop/comppgen_project/SRR32304900_1.fq', "fastq"
for record in SeqIO.parse(PATH +'SRR32477976_test.fasta', 'fasta'):
    seq = str(record.seq)
    if len(seq) <69:
        continue
    for y in range(len(seq)-68):
        result = ht.lookup_69mer(seq[y:(y+69)])
        if result == -2:
            hits['read'].append(record.id)
            hits['sequence'].append(seq[y:(y+69)])
            hits['species'].append('MC predicted')
            print('found predicted!')
            break
        elif result < 0:
            continue
        else:
            temp = seq[y:(y+69)]
            temp = temp.replace('T', 'U')
            if temp in sequences:
                print('found real!')
                hits['read'].append(record.id)
                hits['sequence'].append(seq[y:(y+69)])
                hits['species'].append(result)
                del temp
                break
    if c%100000==0:
        print(c)
    c+=1
del record, seq, result,c,y, sequences


hits = pd.DataFrame(hits)



seqs = list(preprocessed['original sequence'])
for x in range(len(hits)):
    t = hits['sequence'][x].replace('T', 'U')
    if t == (seqs[hits['species'][x]]):
        continue
    elif t in seqs:
        hits['species'][x] = preprocessed[preprocessed['original sequence'] == t].index
        print(hits['species'][x])
    else:
        print('not found')
del x, t,seqs

for x in range(len(hits)):
    i =  hits['species'][x]
    hits['species'][x] = preprocessed['Species'][i]
del i,x
counts = hits['species'].value_counts()
print(counts)

with pd.ExcelWriter(PATH +'SRR32477976_test_hits.xlsx') as writer:
    counts.to_excel(writer, sheet_name = "Summary", index = True)
    hits.to_excel(writer, sheet_name = "All hits", index = False)
