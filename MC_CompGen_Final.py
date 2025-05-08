#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 22:32:22 2025

@author: mattc
"""



def write_csv(path, data):
    text = ''
    for x in range(len(data)):
        text = text + str(data[x]) + ','
    if len(text) > 1:
        text = text[:-1]
    try:
        f = open(path + 'preprocessed.csv', mode="r")
        f = open(path + 'preprocessed.csv', mode="a")
        f.write(text +"\n")
        f.close()   
    except:
        f= open(path + 'preprocessed.csv', mode="w")
        f.write("Species,Original 16s rRNA length,Score,Homologue,Aligned homologue,Aligned consensus,homologue index,consensus index,original sequence\n")
        f.write(text+"\n")
        f.close()
    del f

#consensus 7a and 7b
consensus = "AGGUGNUGCAUGGYYGYCGUCAGCUCGUGYCGUGAGUGUUGGGUUAAGUCCCRYAACGAGCGCAACCCU"
#max score is 63.5 due to R and Y being .25 and N is 0

#consensus 6a and 6b
#consensus = 'AAANTYAAANRAATWGRCGGGGRCCCGCACAAGATGTGGTTTAATTCGA'

count = 0

#This block of code is for if I want to resume where I left off
import pandas as pd
try:
    f = pd.read_csv('/home/mattc/Desktop/comppgen_project/Actinomycetota/Actinomycetota_preprocessed.csv')
    n = len(f)-1
    del f
except:
    n = 0
    
import sys
sys.path.append('/home/mattc/Desktop/comppgen_project/')
import MC_WSFunc as mc
import time

from Bio import SeqIO
for record in SeqIO.parse('/home/mattc/Desktop/comppgen_project/Actinomycetota/Actinomycetota_16S.fasta', "fasta"):
    if count <= n:
        count+=1
        continue
    print("Header:", record.id)
    #print('length: ' +str(len(record.seq)))
    #print("Sequence:", record.seq)
    score, aligned_seq1, aligned_seq2, i ,j = mc.SWA(record.seq, consensus)
    print(score)
    if '-' in aligned_seq2:
        homologue =  mc.fix_gaps(aligned_seq1, aligned_seq2)
    else:
        homologue = aligned_seq1
    homologue = '-'*j + homologue + '-'*(len(consensus)-j-len(homologue))
    #print(len(homologue))
    if len(homologue) != 69:
        print(aligned_seq1)
        print(aligned_seq2)
        print(homologue)
        print(i)
        print(j)
        break
    original = str(record.seq)
    if len(original) > (i+69):
        original = original[i:(i+69)]
    else:
        original = original[i::]
    write_csv('/home/mattc/Desktop/comppgen_project/Actinomycetota/Actinomycetota_', [record.description, len(record.seq), score, homologue, aligned_seq1, aligned_seq2,i,j, original])
    count+=1
    if count%50000==0:
        print(count)
        time.sleep(420)
del score, aligned_seq1, aligned_seq2, n, count, record, homologue, i, j, original

count=0
for record in SeqIO.parse('/home/mattc/Desktop/comppgen_project/cyanobacteriota/cyanobacteriota_16srRNA.fasta', "fasta"):
    print("Header:", record.id)
    #print('length: ' +str(len(record.seq)))
    #print("Sequence:", record.seq)
    score, aligned_seq1, aligned_seq2, i ,j = mc.SWA(record.seq, consensus)
    print(score)
    if '-' in aligned_seq2:
        homologue =  mc.fix_gaps(aligned_seq1, aligned_seq2)
    else:
        homologue = aligned_seq1
    homologue = '-'*j + homologue + '-'*(len(consensus)-j-len(homologue))
    #print(len(homologue))
    if len(homologue) != 69:
        print(aligned_seq1)
        print(aligned_seq2)
        print(homologue)
        print(i)
        print(j)
        break
    original = str(record.seq)
    if len(original) > (i+69):
        original = original[i:(i+69)]
    else:
        original = original[i::]
    write_csv('/home/mattc/Desktop/comppgen_project/cyanobacteriota/cyanobacteriota_', [record.description, len(record.seq), score, homologue, aligned_seq1, aligned_seq2,i,j, original])
    count+=1
    if count%50000==0:
        print(count)
        time.sleep(420)
del score, aligned_seq1, aligned_seq2, count, record, homologue, i, j, original

count=0
for record in SeqIO.parse('/home/mattc/Desktop/comppgen_project/Verrucomicrobiota/Verrucomicrobiota.fasta', "fasta"):
    print("Header:", record.id)
    #print('length: ' +str(len(record.seq)))
    #print("Sequence:", record.seq)
    score, aligned_seq1, aligned_seq2, i ,j = mc.SWA(record.seq, consensus)
    print(score)
    if '-' in aligned_seq2:
        homologue =  mc.fix_gaps(aligned_seq1, aligned_seq2)
    else:
        homologue = aligned_seq1
    homologue = '-'*j + homologue + '-'*(len(consensus)-j-len(homologue))
    #print(len(homologue))
    if len(homologue) != 69:
        print(aligned_seq1)
        print(aligned_seq2)
        print(homologue)
        print(i)
        print(j)
        break
    original = str(record.seq)
    if len(original) > (i+69):
        original = original[i:(i+69)]
    else:
        original = original[i::]
    write_csv('/home/mattc/Desktop/comppgen_project/Verrucomicrobiota/Verrucomicrobiota_', [record.description, len(record.seq), score, homologue, aligned_seq1, aligned_seq2,i,j, original])
    count+=1
    if count%50000==0:
        print(count)
        time.sleep(420)
del score, aligned_seq1, aligned_seq2, count, record, homologue, i, j, original

#mc.SWA(mc.random_nucleotide_sequence(900), consensus)

import MC_WSFunc as mc
import pandas as pd
consensus = "AGGUGNUGCAUGGYYGYCGUCAGCUCGUGYCGUGAGUGUUGGGUUAAGUCCCRYAACGAGCGCAACCCU"

preprocessed = pd.read_csv('/home/mattc/Desktop/comppgen_project/Actinomycetota/Actinomycetota_preprocessed.csv')
preprocessed = preprocessed[~preprocessed['original sequence'].duplicated(keep='first')]
preprocessed.to_csv('/home/mattc/Desktop/comppgen_project/Actinomycetota/Actinomycetota_preprocessed_unique.csv', index = False)

preprocessed = pd.read_csv('/home/mattc/Desktop/comppgen_project/cyanobacteriota/Cyanobacteriota_preprocessed.csv')
preprocessed = preprocessed[~preprocessed['original sequence'].duplicated(keep='first')]
preprocessed.to_csv('/home/mattc/Desktop/comppgen_project/cyanobacteriota/Cyanobacteriota_preprocessed_unique.csv', index = False)

preprocessed = pd.read_csv('/home/mattc/Desktop/comppgen_project/Verrucomicrobiota/Verrucomicrobiota_preprocessed.csv')
preprocessed = preprocessed[~preprocessed['original sequence'].duplicated(keep='first')]
preprocessed.to_csv('/home/mattc/Desktop/comppgen_project/Verrucomicrobiota/Verrucomicrobiota_preprocessed_unique.csv', index = False)

preprocessed = pd.read_csv('/home/mattc/Desktop/comppgen_project/Bacteroidota/Bacteroidota_preprocessed.csv')
preprocessed = preprocessed[~preprocessed['original sequence'].duplicated(keep='first')]
preprocessed.to_csv('/home/mattc/Desktop/comppgen_project/Bacteroidota/Bacteroidota_preprocessed_unique.csv', index = False)
del preprocessed

df1 = pd.read_csv('/home/mattc/Desktop/comppgen_project/Actinomycetota/Actinomycetota_preprocessed_unique.csv')
df1['Source Phylum'] = 'Actinomycetota'
df2 = pd.read_csv('/home/mattc/Desktop/comppgen_project/cyanobacteriota/Cyanobacteriota_preprocessed_unique.csv')
df2['Source Phylum'] = 'Cyanobacteriota'
df3 = pd.read_csv('/home/mattc/Desktop/comppgen_project/Verrucomicrobiota/Verrucomicrobiota_preprocessed_unique.csv')
df3['Source Phylum'] = 'Verrucomicrobiota'
df4 = pd.read_csv('/home/mattc/Desktop/comppgen_project/Bacteroidota/Bacteroidota_preprocessed_unique.csv')
df4['Source Phylum'] = 'Bacteroidota'
preprocessed = pd.concat([df4, df1, df2, df3])
preprocessed = preprocessed[~preprocessed['original sequence'].duplicated(keep='first')]
del df1,df2,df3, df4
preprocessed.to_csv('/home/mattc/Desktop/comppgen_project/All_preprocessed_unique.csv', index = False)

sources = preprocessed['Source Phylum'].value_counts()
import matplotlib.pyplot as plt
fig1, ax1 = plt.subplots()
total = len(preprocessed)
def my_fmt(x):
    """
    from here:
    https://stackoverflow.com/questions/59644751/show-both-value-and-percentage-on-a-pie-chart

    """
    print(round(x,1))
    return '{:.1f}%\n({:.0f})'.format(round(x,1), total*x/100)
ax1.pie(list(sources), labels=list(sources.index), autopct=my_fmt,
        shadow=False, startangle=90)
# Equal aspect ratio ensures that pie is drawn as a circle.
ax1.axis('equal')  
# Adding a title
plt.title('Source phyla of unique 16s rRNA consensus alignments, before filtering')
# Displaying the chart
plt.show()
del ax1, fig1

plt.hist(preprocessed["Score"], bins=50, color='skyblue', edgecolor='black')
plt.title('Smith Waterman Scores of aligned 16s rRNA genes to consensus sequence, before filtering')
plt.xlabel('Score')
plt.ylabel('Frequency')
plt.show()

import seaborn as sns
sns.boxplot(x='Source Phylum', y='Score', data=preprocessed)
plt.xticks(fontsize=9) 
plt.title('Smith Waterman Scores of aligned 16s rRNA genes to consensus sequence, before filtering')
plt.show()

preprocessed = preprocessed[preprocessed["Score"] >14] 
preprocessed.to_csv('/home/mattc/Desktop/comppgen_project/All_preprocessed_unique.csv', index = False)


plt.hist(preprocessed["Score"], bins=50, color='skyblue', edgecolor='black')
plt.title('14 <= Smith Waterman Scores of aligned 16s rRNA genes to consensus sequence')
plt.xlabel('Score')
plt.ylabel('Frequency')
plt.show()

c = 0
colors = ['blue', 'green', 'red', 'yellow']
sources = sources.sort_values(ascending=False)
for x in list(sources.index):
    subset = preprocessed[preprocessed['Source Phylum']==x]
    #plt.hist(subset['Score'], bins=30, alpha=0.5, label=x)
    sns.histplot(data=subset['Score'], color=colors[c], alpha=1, kde=True, label=x)
    c+=1
del x, subset,c,colors, sources
# Add labels and title
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('14 <= Smith Waterman Score of consensus sequence by phylum')
# Add legend
plt.legend(loc='upper left')
# Show the plot
plt.show()

sns.boxplot(x='Source Phylum', y='Score', data=preprocessed)
plt.xticks(fontsize=9) 
plt.title('14 <= Smith Waterman Scores of aligned 16s rRNA genes to consensus sequence')
plt.show()

sources = preprocessed['Source Phylum'].value_counts()
fig1, ax1 = plt.subplots()
total = len(preprocessed)
def my_fmt(x):
    """
    from here:
    https://stackoverflow.com/questions/59644751/show-both-value-and-percentage-on-a-pie-chart

    """
    print(round(x,1))
    return '{:.1f}%\n({:.0f})'.format(round(x,1), total*x/100)
ax1.pie(list(sources), labels=list(sources.index), autopct=my_fmt,
        shadow=False, startangle=90)
# Equal aspect ratio ensures that pie is drawn as a circle.
ax1.axis('equal')  
# Adding a title
plt.title('Source phyla of unique 16s rRNA consensus alignments')
# Displaying the chart
plt.show()
del ax1, fig1

seq_list = list(set(preprocessed["original sequence"]))
ratios, indices = mc.consensus_match_ratios(list(preprocessed['Homologue']), consensus)
mc.plot_match_ratios(consensus, ratios, ' filtered aligned 16s rRNA genes')

insert_ratios = mc.find_inserts(list(preprocessed["Aligned consensus"]), consensus)

preprocessed2 = preprocessed.iloc[0:0]
for x in range(len(preprocessed)): 
    temp = preprocessed.iloc[[x]]
    seq = list(temp['original sequence'])[0]
    if (not isinstance(seq, str)):
        continue
    if 'N' in seq:
        for z in ['A', 'U', 'G', 'C']:
            seq2 = seq.replace('N', z)
            temp['original sequence'] = [seq2]
            preprocessed2 = pd.concat([preprocessed2, temp], ignore_index=True)
    else:
        preprocessed2 = pd.concat([preprocessed2, temp], ignore_index=True)
del preprocessed, indices, seq, seq2, x, z, temp
preprocessed2 = preprocessed2[~preprocessed2['original sequence'].duplicated(keep='first')]
preprocessed2.to_csv('/home/mattc/Desktop/comppgen_project/All_preprocessed_unique.csv', index = False)
del preprocessed2


hmm = mc.HMM(ratios, consensus, scaling = 1.2)
possible = hmm.generate_n(400000)
possible = possible[~possible['sequence'].duplicated(keep='first')]

hmm = mc.HMM(ratios, consensus, scaling = 1.05)
possible2 = hmm.generate_n(400000)
possible2 = possible2[~possible2['sequence'].duplicated(keep='first')]

hmm = mc.HMM(ratios, consensus, scaling = 1.01)
possible3 = hmm.generate_n(400000)
possible3 = possible3[~possible3['sequence'].duplicated(keep='first')]

hmm = mc.HMM(ratios, consensus, scaling = 1.005)
possible4 = hmm.generate_n(400000)
possible4 = possible4[~possible4['sequence'].duplicated(keep='first')]

possible = pd.concat([possible, possible2, possible3, possible4])
del possible2, possible3, possible4,hmm


ratios, indices = mc.consensus_match_ratios(list(possible['sequence']), consensus)
mc.plot_match_ratios(consensus, ratios, ' predicted 16s rRNA genes (no inserts)')
del indices, ratios

possible_inserts = mc.make_inserts(list(possible['sequence']), insert_ratios)
possible = possible[~possible['sequence'].duplicated(keep='first')]

possible_inserts = possible_inserts[~possible_inserts['sequence'].duplicated(keep='first')]

possible = pd.concat([possible, possible_inserts], ignore_index=True)
del possible_inserts, consensus

seqs = list(possible['sequence'])
filtered = []
c=0
for i in range(len(seqs)):
    x = seqs[i]
    flag = False
    if len(x) ==70:
        x = x[:-1]
    if x in seq_list:
        flag = True
        print('filtered!')
    if flag == False:
        filtered.append(i)
    if c%50000==0:
        print(c)
    c+=1

possible = possible.loc[filtered]
del filtered, x,c,flag, seq_list,i

possible = possible[~possible['sequence'].duplicated(keep='first')]
possible.to_csv('/home/mattc/Desktop/comppgen_project/possible_sequences.csv', index = False)
del possible

import pandas as pd
import sys
sys.path.append('/home/mattc/Desktop/comppgen_project/')
import MC_WSFunc as mc
import time
from Bio import SeqIO
from MC_hash_table import DNAHashTable
ht = DNAHashTable(2**25)

preprocessed = pd.read_csv('/home/mattc/Desktop/comppgen_project/All_preprocessed_unique.csv')
for x in range(len(preprocessed)): 
    temp = preprocessed['original sequence'][x]
    temp = temp.replace('U', 'T')
    ht.insert_69mer(temp, x)  # store index of a 69-mer
del x, temp


possible = pd.read_csv('/home/mattc/Desktop/comppgen_project/possible_sequences.csv')
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
for record in SeqIO.parse('/home/mattc/Desktop/comppgen_project/SRR32477976.fasta', 'fasta'):
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

with pd.ExcelWriter('/home/mattc/Desktop/comppgen_project/SRR32477976_hits.xlsx') as writer:
    counts.to_excel(writer, sheet_name = "Summary", index = True)
    hits.to_excel(writer, sheet_name = "All hits", index = False)