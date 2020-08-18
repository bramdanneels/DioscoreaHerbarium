#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Bram Danneels (bram.danneels@ugent.be)

Script for performing metagenome assembly of trimmed SE reads

usage: python SE_assembly.py trimmed_reads.fastq kraken_db nr_cores

required software:
    R (https://www.r-project.org/)
    SMALT (https://www.sanger.ac.uk/tool/smalt-0/)
    SAMtools (https://github.com/samtools/)
    SPAdes (https://cab.spbu.ru/software/spades/)
    Kraken (https://ccb.jhu.edu/software/kraken/)
    Blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  
input:
    trimmed_reads: fastq file of (trimmed) single-ended Illumina sequencing reads, all gathered in one .fastq file
    kraken_db: location of kraken database for classifying metagenomic contigs
    nr_cores: number of cores for running spades, mapping, and running kraken
output:
    folder SE_assembly containing:
        SPAdes output of the metagenome assembly
        Summary files of coverage and GC-content of metagenomic contigs, necessary for further processing (see SE_assembly_part2.py script)
        Plot where every conting > 500 bp are plotted based on their read coverage and GC content, and colored based on taxonomy      
"""

import sys  
from Bio import SeqIO
from Bio.SeqUtils import GC
import subprocess as sp
import matplotlib.pyplot as plt
import numpy as np

def spades(reads, cores):
    '''
    Performs Spades Assembly on the trimmed, paired reads
    '''
    sp.call('spades.py -s {} -k 21,25,33,37,45 -t {} -o SE_assembly'.format(reads, cores), shell=True)
    return()

def calc_stats(ref, reads, cores):
    '''
    Calculates coverage and GC content of contigs
    '''
    print('Calculating GC content\n')
    contigs = SeqIO.parse(open(ref, 'rU'), 'fasta')
    GCOUT = open('GCstats.txt', 'w')
    print('Group.1\tGC', file=GCOUT)
    for contig in contigs:
        print('\t'.join([str(contig.id), str(GC(contig.seq))]), file=GCOUT)
    GCOUT.close()
    
    print('Mapping reads to contigs...\n')
    sp.call('smalt index -k 13 -s 4 smalt-index {}'.format(ref), shell=True)
    sp.call('smalt map -n {} -f sam -o smalt_mapping.sam smalt-index {}'.format(cores, reads), shell=True)
    sp.call('samtools sort smalt_mapping.sam -O sam -T temp -o mapping.sorted.sam', shell=True)
    sp.call('samtools depth -a mapping.sorted.sam > samtools_depth.txt', shell=True)
    print('Summarizing statistics...\n')
    sp.call('Rscript ~/scripts/assembly/mapping_statsGC.R', shell=True)
    sp.call('rm smalt_mapping.sam', shell=True)
    return()

def run_kraken(contigs, kraken_db, cores):
    '''
    Runs Kraken on the metagenome contigs
    '''
    print('Classifying metagenome contigs...\n')
    sp.call('kraken --preload --threads {} --db {} {} > SE_assembly/meta_kraken_out'.format(cores,kraken_db, contigs), shell=True)
    sp.call('kraken-translate --db {} > SE_assembly/meta_kraken_lab'.format(kraken_db), shell=True)
    print('Finished classifying metagenome contigs!\n')
    return()

def plot_kraken(stats, kraken_lab, cutoff):
    '''
    Create plots for visually inspecting metagenome contigs
    '''
    print('Creating annotated plots...\n')
    label_dict = {}
    for line in open(kraken_lab):
        contig, tax = line.strip().split('\t')
        tax = tax.split(';')
        if 'Burkholderia' in tax:
            label_dict[contig] = 'red'
        elif 'Betaproteobacteria' in tax:
            label_dict[contig] = 'orange'
        elif 'Bacteria' in tax:
            label_dict[contig] = 'yellow'
        elif 'Viridiplantae' in tax:
            label_dict[contig] = 'green'
        elif 'Eukaryota' in tax:
            label_dict[contig] = 'blue'
        else:
            label_dict[contig] = 'grey'
        
    labeldict={'red':'Burkholderia', 'orange':'Betaproteobacteria', 'yellow':'Bacteria', 'green':'Viridiplantae', 'blue':'Eukaryota', 'grey':'Unclassified'}
    dotdict = {'red':([],[]), 'orange':([],[]), 'yellow':([],[]), 'green':([],[]), 'blue':([],[]), 'grey':([],[])}
    for line in open(stats):
        if 'contig' in line:
            pass
        else:
            contig, length, cov, gc = line.strip().split('\t')
            if contig not in label_dict.keys():
                label_dict[contig] = 'grey'
            if int(length) >= cutoff:
                color = label_dict[contig]
                dotdict[color][0].append(float(gc))
                dotdict[color][1].append(np.log10(float(cov)))
    
    for color in ['grey', 'blue', 'green', 'yellow', 'orange', 'red']:
        plt.scatter(dotdict[color][0], dotdict[color][1], color=color, s=10, label=labeldict[color])
    plt.title('%G+C vs Coverage, contigs >= {} nt'.format(str(cutoff)))
    plt.xlabel('% G+C')
    plt.ylabel('Log10(Coverage)')
    plt.legend(fontsize='xx-small')
    plt.savefig('SE_assembly/Cov_vs_GC_kraken.png', dpi=600, format='png')
    plt.close()
    
    dotdict = {'red':([],[]), 'orange':([],[]), 'yellow':([],[]), 'green':([],[]), 'blue':([],[]), 'grey':([],[])}
    for line in open(stats):
        if 'contig' in line:
            pass
        else:
            contig, length, cov, gc = line.strip().split('\t')
            if contig not in label_dict.keys():
                label_dict[contig] = 'grey'
            if int(length) >= cutoff:
                color = label_dict[contig]
                dotdict[color][0].append(int(length))
                dotdict[color][1].append(float(gc))
    
    for color in ['grey', 'blue', 'green', 'yellow', 'orange', 'red']:
        plt.scatter(dotdict[color][0], dotdict[color][1], color=color, s=10, label=labeldict[color])
    plt.title('%G+C vs Contig Length, contigs >= {} nt'.format(str(cutoff)))
    plt.xlabel('Contig Length')
    plt.ylabel('% G+C')
    plt.legend(fontsize='xx-small')
    plt.savefig('SE_assembly/GC_vs_Len_kraken.png', dpi=600, format='png')
    print('Finished creating plots!\n')
    return()
    
script, reads, kraken_db, cores = sys.argv

if int(assembly) != 0:
	print('Running first assembly...\n')
	spades(reads, cores)
	print('Assembly finished\n')
	print('Cleaning up assembly files...\n')
sp.call('mv SE_assembly/contigs.fasta ./',shell=True)
sp.call('rm -r SE_assembly/*', shell=True)
sp.call('mv contigs.fasta SE_assembly/', shell=True)
print('Calculating statistics...\n')
calc_stats('SE_assembly/contigs.fasta', reads, cores)
sp.call('mv contig_statistics.txt SE_assembly/', shell=True)
sp.call('mv contig_coverage.pdf SE_assembly/', shell=True)
sp.call('rm GCstats.txt mapping.* samtools_depth.txt smalt-index*', shell=True)
print('Finished calculating statistics\n')
run_kraken('SE_assembly/contigs.fasta', cores)
plot_kraken('SE_assembly/contig_statistics.txt', 'SE_assembly/meta_kraken_lab', 500)
print('SE assembly finished\n')




