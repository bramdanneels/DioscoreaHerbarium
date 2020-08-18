#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Bram Danneels (bram.danneels@ugent.be)

Script to summarize nucleotide misincorporation data from MapDamage2 (https://ginolhac.github.io/mapDamage/)

input: "misincorporation.txt" file(s) from MapDamage2 output
output: table containing the average frequency of specific nucleotide misincorporations (columns), per input sample/file (rows)

usage: python avg_conversions.py inputfile1 inputfile2 inputfile3 ...
"""

from sys import argv
import numpy as np

script, *files = argv

conv_dict = {'A':4, 'C':5, 'G':6, 'T':7, '-':0, 'S':-1}

print_bool = False
for file in files:
    file_dict = {}
    name = file.split('.')[0]
    for line in open(file):
        if not line.startswith('#'):
            if line.startswith('Chr'):
                headers = line.strip().split()
            else:
                info = line.strip().split()
                end, orient, pos = info[1:4]
                if end == '5p' and orient == '+':
                    for mut, count in zip(headers[9:], info[9:]):
                        ref = conv_dict[mut[0]]
                        if ref > 0:
                            if int(info[ref]) > 0:
                                freq = int(count)/int(info[ref])
                            else:
                                freq = 0
                        elif ref < 0:
                            if int(info[8]) > 0:
                                freq = int(count)/int(info[8])
                            else:
                                freq = 0
                        else:
                            freq = count
                        file_dict[mut] = file_dict.get(mut, [])
                        file_dict[mut].append(freq)
    heads = ['name']
    freq_string = [name]
    for mut in sorted(file_dict.keys()):
        freqs = file_dict[mut]
        av_freq = np.mean([float(x) for x in freqs])
        heads.append(mut)
        freq_string.append(av_freq)
    if print_bool == False:
        print(*heads)
        print_bool = True
    print(*freq_string)
                            
                            
                    
                    