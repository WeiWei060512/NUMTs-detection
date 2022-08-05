#!/usr/bin/env python

################################################################################
## This script generates the nuclear reference genome file used in enrichment permutation test
################################################################################

import fileinput
import sys, os
import pandas as pd
import glob
import scipy.stats as stats
import numpy as np

### this script creates a fake whole genome bed file for chromosome enrichment analysis

input1 = sys.argv[1] # bedfile of each chromosome from whole genome

df_raw = pd.read_csv(input1, names=['chr','start','end'],sep="\t")
df = df_raw[df_raw.chr.isin(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'])]

refStart_noChr = []
refEnd_noChr = []
for index, row in df.iterrows():
        if row['chr'] == "chr1":
                refStart_noChr.append(row['start'])
                refEnd_noChr.append(row['end'])
        elif row['chr'] == "chr2":
                refStart_noChr.append(row['start'] + 230481121)
                refEnd_noChr.append(row['end'] + 230481121)
        elif row['chr'] == "chr3":
                refStart_noChr.append(row['start'] + 471029362)
                refEnd_noChr.append(row['end'] + 471029362)
        elif row['chr'] == "chr4":
                refStart_noChr.append(row['start'] + 669129506)
                refEnd_noChr.append(row['end'] + 669129506)
        elif row['chr'] == "chr5":
                refStart_noChr.append(row['start'] + 858882173)
                refEnd_noChr.append(row['end'] + 858882173)
        elif row['chr'] == "chr6":
                refStart_noChr.append(row['start'] + 1040147551)
                refEnd_noChr.append(row['end'] + 1040147551)
        elif row['chr'] == "chr7":
                refStart_noChr.append(row['start'] + 1210226075)
                refEnd_noChr.append(row['end'] + 1210226075)
        elif row['chr'] == "chr8":
                refStart_noChr.append(row['start'] + 1369196210)
                refEnd_noChr.append(row['end'] + 1369196210)
        elif row['chr'] == "chr9":
                refStart_noChr.append(row['start'] + 1513964346)
                refEnd_noChr.append(row['end'] + 1513964346)
        elif row['chr'] == "chr10":
                refStart_noChr.append(row['start'] + 1635754899)
                refEnd_noChr.append(row['end'] + 1635754899)
        elif row['chr'] == "chr11":
                refStart_noChr.append(row['start'] + 1769017905)
                refEnd_noChr.append(row['end'] + 1769017905)
        elif row['chr'] == "chr12":
                refStart_noChr.append(row['start'] + 1903551647)
                refEnd_noChr.append(row['end'] + 1903551647)
        elif row['chr'] == "chr13":
                refStart_noChr.append(row['start'] + 2036689468)
                refEnd_noChr.append(row['end'] + 2036689468)
        elif row['chr'] == "chr14":
                refStart_noChr.append(row['start'] + 2134672596)
                refEnd_noChr.append(row['end'] + 2134672596)
        elif row['chr'] == "chr15":
                refStart_noChr.append(row['start'] + 2225240745)
                refEnd_noChr.append(row['end'] + 2225240745)
        elif row['chr'] == "chr16":
                refStart_noChr.append(row['start'] + 2309882073)
                refEnd_noChr.append(row['end'] + 2309882073)
        elif row['chr'] == "chr17":
                refStart_noChr.append(row['start'] + 2391688017)
                refEnd_noChr.append(row['end'] + 2391688017)
        elif row['chr'] == "chr18":
                refStart_noChr.append(row['start'] + 2474608233)
                refEnd_noChr.append(row['end'] + 2474608233)
        elif row['chr'] == "chr19":
                refStart_noChr.append(row['start'] + 2554697883)
                refEnd_noChr.append(row['end'] + 2554697883)
        elif row['chr'] == "chr20":
                refStart_noChr.append(row['start'] + 2613138641)
                refEnd_noChr.append(row['end'] + 2613138641)
        elif row['chr'] == "chr21":
                refStart_noChr.append(row['start'] + 2677082909)
                refEnd_noChr.append(row['end'] + 2677082909)
        elif row['chr'] == "chr22":
                refStart_noChr.append(row['start'] + 2717171532)
                refEnd_noChr.append(row['end'] + 2717171532)
        elif row['chr'] == "chrX":
                refStart_noChr.append(row['start'] + 2756331314)
                refEnd_noChr.append(row['end'] + 2756331314)
        elif row['chr'] == "chrY":
                refStart_noChr.append(row['start'] + 2911224348)
                refEnd_noChr.append(row['end'] + 2911224348)
        else:
                print("Something is wrong!")

df['refStart_noChr'] = refStart_noChr
df['refEnd_noChr'] = refEnd_noChr
df.to_csv(input1 + '.convert.noGap.tsv', index=False,sep = '\t')


