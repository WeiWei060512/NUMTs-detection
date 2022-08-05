#!/usr/bin/env python

import Bio
from Bio import AlignIO
from Bio import SeqIO
import glob
from Bio.Align.Applications import ClustalwCommandline
import os
import pandas as pd
import sys
import re

# input aln - multiple alignment file from NUMT sequences
# input numts - target numts contains numtID, numt_start and numt_end
aln,numts = sys.argv[1:]

# make sampleID, numtID and contigID
index1 = re.sub(r".*LP", "LP", aln)
index2 = re.sub(r"\..*", "", index1)
numtID = index2.split('_')[2] + '_' + index2.split('_')[3]
sampleID = index2.split('_')[0] + '_' + index2.split('_')[1]
index3 = re.sub(r".human.*", "", index1)
contigID = re.sub(r".*\.", "", index3)

mynumts = pd.read_csv(numts, names=['numtID','numt_start','numt_end','comment'], sep="\t")
mysequence =  AlignIO.read(aln,'fasta')

myseq = []
for i in mysequence:
	myseq.append(i.seq)

myout_list = []
myseq_zip = list(zip(*myseq))
for i in range(0, (len(myseq_zip))):
	myout_list.append(myseq_zip[i])

mydf1 = pd.DataFrame(myout_list)
mydf1.columns = ['ref','alt']
mydf1['pos'] = mydf1.index + 1 

## look for the starting position
for index, row in mydf1.iterrows():
    if row['alt']!='-':
        start_pos = row['pos']
        break

mydf1.rev = mydf1.sort_values(by='pos', ascending=False)
for index, row in mydf1.rev.iterrows():
    if row['alt']!='-':
        end_pos = row['pos']
        break

mydf2 = mydf1[ (mydf1['pos'] >= start_pos) & (mydf1['pos'] <= end_pos) & (mydf1['pos'] <= 16569)]
mydf2['sampleID'] = sampleID
mydf2['numtID'] = numtID
mydf2['contigID'] = contigID

mydf2['numtLength'] = len(mydf2)
mydf2 = pd.merge(mydf2, mynumts, on='numtID', how='left')

mydf3 = mydf2[~ (mydf2['ref']==mydf2['alt'])]
mydf3_dedup = mydf3.drop_duplicates(subset=['ref','alt','pos'])
mydf3['varCount'] = len(mydf3_dedup)
mydf2['varCount'] = len(mydf3_dedup)

mydf4 = mydf3[ (mydf3['pos'] > mydf3['numt_start']) & (mydf3['pos'] < mydf3['numt_end'])]

mydf2.to_csv(aln + '.numt.tsv', index=0, sep='\t')
mydf3_dedup.to_csv(aln + '.numtVar.tsv', index=0, sep='\t')
mydf4.to_csv(aln + '.numtVarFilterPos.tsv', index=0, sep='\t')

