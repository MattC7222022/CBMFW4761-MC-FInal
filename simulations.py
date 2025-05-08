#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 00:03:37 2025

@author: mattc
"""


def write_simulations(path, data):
    text = ''
    for x in range(len(data)):
        text = text + str(data[x]) + ','
    if len(text) > 1:
        text = text[:-1]
    try:
        f = open(path + 'simulations.csv', mode="r")
        f = open(path + 'simulations.csv', mode="a")
        f.write(text +"\n")
        f.close()   
    except:
        f= open(path + 'simulations.csv', mode="w")
        f.write("Random nucleotide length,Mean score,Median score,Standard deviation,Ninety fifth percentile\n")
        f.write(text+"\n")
        f.close()
    del f
    
path = '/home/mattc/Desktop/comppgen_project/'
consensus = "AGGUGNUGCAUGGYYGYCGUCAGCUCGUGYCGUGAGUGUUGGGUUAAGUCCCRYAACGAGCGCAACCCU"

import sys
sys.path.append('/home/mattc/Desktop/comppgen_project/')
import MC_WSFunc as mc
import numpy as np
for n in range(400,1925,25):
    scorelist = []
    for z in range(25):
        seq = mc.random_nucleotide_sequence(n)
        score, aligned_seq1, aligned_seq2, i ,j = mc.SWA(seq, consensus)
        scorelist.append(score)
    print("done with " +str(n))
    write_simulations(path, [n, np.mean(scorelist), np.median(scorelist), np.std(scorelist), np.percentile(scorelist,95)])
