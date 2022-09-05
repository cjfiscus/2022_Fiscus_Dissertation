#!/usr/bin/env python
import sys
import re
import os
import statistics
from Bio import SeqIO

#By Dan Koenig 2018
#Script to use Kmer table to calculate median values for a genomic sequence
#Each kmer in the fasta sequence will be assigned the value from the kmer table and then a global median will be taken for the sequences
#Kmers not found in the table will be ignored (this will include all kmers that have ambigous nucleotides if not found in the table)

#command line inputs:
#1: kmer table file
#Input kmer table format is tab delimited two column 1: kmer 2: statistic
#f = gzip.open(sys.argv[1], 'r')
#2: fasta file
#3: length of kmer to analyze (should be consistent with kmer table)
#Returned values:
#1: Chromosome
#2: Region start position
#3: Region end position
#4: median value for the region
#5: number of kmers used for #4
#If no kmers are found that match the kmer table NA will be returned for the median value


#################################################################################################
#################################################################################################
#############################FUNCTIONS###########################################################
#################################################################################################
#################################################################################################
#Function to reverse complement sequences in the input kmer table
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
#Function to process individual fasta
def processthisseq (fastaseq,merlen,keyindex,normd):
    fastaseq = fastaseq.upper()
    sum_score = []
    kmercount = 0
    for i in range(0,len(fastaseq)):
        if len(fastaseq[i:i + merlen]) == merlen and fastaseq[i:i+merlen] in keyindex:
            kmercount += 1
            sum_score.append(float(keyindex[fastaseq[i:i+merlen]]))
    if kmercount >0:
        return [statistics.median(sum_score),kmercount]
    else:
        return ["NA",kmercount]

#################################################################################################
#################################################################################################
#############################COMMANDLINE PARAMETERS##############################################
#################################################################################################
#################################################################################################
tablefile = sys.argv[1]
g = sys.argv[2]
merlen = int(sys.argv[3])

#################################################################################################
#################################################################################################
#############################KMER TABLE PARSE AND STORE##########################################
#################################################################################################
#################################################################################################
#create a dictionary with kmer keys and table values from the input table
#Will reverse complment every kmer and add it to the dictionary with the same value
keydata = {}
fulldata = 0
with open(tablefile,'r') as f:
    for line in f:
        if not line: break
        line = line.rstrip('\n')
        linesplit = line.split()
        fulldata += float(linesplit[1])
        keydata[str(linesplit[0])] = linesplit[1]
        keydata[revcompl(str(linesplit[0]))] = linesplit[1]
f.close
#################################################################################################
#################################################################################################
#############################FASTA PARSE AND SCORE CALCULATIONS##################################
#################################################################################################
#################################################################################################
with open(g, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seqscore = processthisseq(record.seq,merlen,keydata,fulldata)
        if seqscore[1] > 0:
            print (record.id + '\t' + str(round(seqscore[0],4)) + "\t" + str(seqscore[1]))
        else:
            print (record.id+"\t"+"NA"+"\t"+"0")
