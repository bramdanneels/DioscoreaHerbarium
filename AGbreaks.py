#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Bram Danneels (bram.danneels@ugent.be)

Script to summarize purine enrichment statistics from MapDamage2 output

input: "dnacomp.txt" file(s) generated from MapDamage2 (https://ginolhac.github.io/mapDamage/)
output: table with in the columns statistics related to purine enrichment around strand breaks, with each row representing a sample/input file

usage: python AGbreaks.py inputfile1 inputfile2 ...

abbreviations:
A-1/A-5: number of adenosines at position -1/-5 of a read start (proxy for strand break)
G-1/G-5: number of guanines at position -1/-5 of a read start (proxy for strand break)
AG-1/AG-5: number of purines (adenosine + guanine) at position -1/-5 of a read start
relA/G/AG: relative difference of adenosine/guanine/both between position -1 compared to -5
absA/G/AG: absolute difference of adenosine/guanine/both between position -1 compared to -5
"""

from sys import argv

script, *files = argv

print('Sample', 'A-1', 'G-1', 'AG-1', 'A-5', 'G-5', 'AG-5', 'relA', 'relG', 'relAG', 'absA', 'absG', 'absAG')
for file in files:
    name = file.split('.')[0]
    A1counts = []
    G1counts = []
    AG1counts = []
    A5counts = []
    G5counts = []
    AG5counts = []
    for line in open(file):
        if not line.startswith('#'):
            cont, end, std, pos, a, c, g, t, total = line.strip().split()
            if end == '5p' and pos == '-1':
                if int(total) != 0:
                    A = int(a)/int(total)
                    A1counts.append(A)
                    G = int(g)/int(total)
                    G1counts.append(G)
                    AG = (int(a)+int(g))/int(total)
                    AG1counts.append(AG)
                else:
                    A1counts.append(0)
                    G1counts.append(0)
                    AG1counts.append(0)
            if end == '5p' and pos == '-5':
                if int(total) != 0:
                    A = int(a)/int(total)
                    A5counts.append(A)
                    G = int(g)/int(total)
                    G5counts.append(G)
                    AG = (int(a)+int(g))/int(total)
                    AG5counts.append(AG)
                else:
                    A1counts.append(0)
                    G1counts.append(0)
                    AG1counts.append(0)
    A1avg = sum(A1counts)/len(A1counts)
    G1avg = sum(G1counts)/len(G1counts)
    AG1avg = sum(AG1counts)/len(AG1counts)
    A5avg = sum(A5counts)/len(A5counts)
    G5avg = sum(G5counts)/len(G5counts)
    AG5avg = sum(AG5counts)/len(AG5counts)
    print(name, A1avg, G1avg, AG1avg, A5avg, G5avg, AG5avg, A1avg/A5avg, G1avg/G5avg, AG1avg/AG5avg, A1avg-A5avg, G1avg-G5avg, AG1avg-AG5avg)
                