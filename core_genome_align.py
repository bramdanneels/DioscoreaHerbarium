#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@authors: Bram Danneels (bram.danneels@ugent.be) & Aurélien Carlier

Scripts that creates single-copy core genome alignments from orthofinder (https://github.com/davidemms/OrthoFinder) output
It is required that all input gene ID's follow the format XXX_YYY where XXX is a unique gene identifier per input genome, and YYY can be anything

usage: python orthologs, fasta_prot, fasta_nucl

required software:
    muscle (https://www.drive5.com/muscle/)
    t-coffee (https://github.com/cbcrg/tcoffee)
    trimal (http://trimal.cgenomics.org/)

input:
    orthologs: "Orthogroups.txt" file from Orthofinder output
    fasta_prot: protein fasta file including all protein sequences used in the Orthofinder run
    fasta_nucl: CDS fasta file including the CDS sequence of all proteins used in the Orthofinder run
output:
    nucl_alignment.phy: nucleotide alignment of the single-copy core genome
    partitions.txt: a file denoting the boundaries of each separate gene in the alignment, used for example by RaxML for creating phylogenetic trees
"""

from __future__ import print_function
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from sys import argv
import subprocess as sp
from Bio import AlignIO

#usage: python whole_genome_align.py orthomcl_out.txt prots.fa nucl.fa 
#orthomcl_out.txt is a set of one_to_one orthologous genes following the same format as the output of orthomcl v1.2
#prots.fa is a fasta file of protein sequences for all the genomes used in the orthomcl run. The sequence identifiers have to match the identifiers in orthomcl_out.txt
#nucl.fa is as prots.fa but with nucleotide sequences. This file is used as the reference for the back translation of the alignments.

script, orthologs, fasta_prot, fasta_nucl = argv

fasta_prot_dict = SeqIO.index(fasta_prot, "fasta")
fasta_nucl_dict = SeqIO.index(fasta_nucl, "fasta")

def fasta_extract(fasta_dict,IDlist):
    #where IDlist is a list of tuples ("species","ID")
    protlist = []
    #sorting list of IDs per species name
    #IDlist = sorted(IDlist, key= lambda species: species[0])
    for x in IDlist:
        species = x[0]
        seq = fasta_dict[x[1]].seq
        seqrec = SeqRecord(seq, name=species, id=species)
        protlist.append(seqrec)   
    return protlist

print('Calculating maximum number of taxa...')
max_taxa = 0   
for line in open(orthologs):
    if line.strip() != '':
        group, genes = line.strip().split(':')
        genes = genes.split()
        nr_taxa = len(set([gene.split('_')[0] for gene in genes]))
        max_taxa = max(max_taxa, nr_taxa)
print('Max number of taxa = {}'.format(max_taxa))
print()

print('Aligning core genes...')
align_nucllist = []
for line in open(orthologs):
    if line.strip() != '':
        species_list = []
        group, genes = line.strip().split(':')
        genes = genes.split()
        nr_taxa = len(set([gene.split('_')[0] for gene in genes]))  
        if max_taxa == nr_taxa and nr_taxa == len(genes):
            for gene in genes:
                species = gene.split('_')[0]
                species_list.append((species, gene))
            protlist = fasta_extract(fasta_prot_dict, species_list)
            OUT = open('temp_prot.fa', 'w')
            fasta_out = FastaIO.FastaWriter(OUT, wrap=None)
            fasta_out.write_file(protlist)
            OUT.close()
            OUT = open('temp.fa', 'w')
            nucllist = fasta_extract(fasta_nucl_dict, species_list)
            fasta_out = FastaIO.FastaWriter(OUT, wrap=None)
            fasta_out.write_file(nucllist)
            OUT.close()
            OUT = open('temp_nucl.fa', 'w')
            IN = open('temp.fa')
            for line in IN:
                if line.startswith('>'):
                    print(line.strip(), file=OUT)
                else:
                    print(line.strip()[:-3].upper(), file=OUT)
            IN.close()
            OUT.close()
            sp.call('muscle -in temp_prot.fa -out temp.clw -clwstrict -quiet', shell=True)
            sp.call('t_coffee -other_pg seq_reformat -in temp_nucl.fa -in2 temp.clw -action +thread_dna_on_prot_aln -output clustalw > temp.clustal', shell=True)
            sp.call('~/Software/trimal/source/trimal -in temp.clustal -fasta -out temp.nucl.fasta -gt 0.5', shell=True)
            nucl_align = AlignIO.read('temp.nucl.fasta', 'fasta')
            nucl_align.sort()
            align_nucllist.append(nucl_align)

sp.call('rm temp*', shell=True)

print('Writing final alignment and partitions file...')
partitions = open("partitions.txt","w")
nucl_align_all = align_nucllist[0]
part = 1
start = 1
for temp_nucl_align in align_nucllist:
    partitions.write('DNA, p{}={}-{}\n'.format(part, start, start+len(temp_nucl_align[0])-1))
    if part != 1:
        nucl_align_all += temp_nucl_align
    part += 1
    start += len(temp_nucl_align[0])
AlignIO.write(nucl_align_all, 'nucl_alignment.phy', 'phylip')

