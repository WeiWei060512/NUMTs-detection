#!/usr/bin/env python
################################################################################
## This script does the enrichment permutation test
## Example of run: python enrichment_simulation.py inputFile 100 10 200 outputFile
################################################################################

import fileinput
import sys, os
import pandas as pd
import glob
import scipy.stats as stats
import numpy as np
import random
import pybedtools 
from pybedtools import BedTool

input1,numtsNo,numtsTargetNo,numtsLength,output = sys.argv[1:]

numtsNo = int(numtsNo) #total number of NUMTs
numtsTargetNo = int(numtsTargetNo) #number of NUMTs in the target/testing region 
numtsLength = int(numtsLength) #length of flaking region of NUMTs  

simuRuns = 1001 #number of simulation
numtsTargetPerc = numtsTargetNo/numtsNo
df_ref = pd.read_csv(input1, sep="\t") #bedfile of target/testing regions
df_ref['chr'] = 'chr'
df_ref = df_ref[['chr','refStart_noChr','refEnd_noChr']]
df_ref = pybedtools.BedTool.from_dataframe(df_ref)

freq_list=list()
for i in range(1, simuRuns):
    randomStart = random.sample(range(1, 2937639397), numtsNo) # nuclear genome    
    randomStart = random.sample(range(1, 16570), numtsNo) # for mtDNA genome

    randomNUMTs = pd.DataFrame(randomStart, columns=['randomStart'])
    randomNUMTs['randomStart'] = randomNUMTs['randomStart'].astype(int)
    randomNUMTs['randomEnd'] = randomNUMTs['randomStart'] + numtsLength
    randomNUMTs['chr'] = "chr"
    randomNUMTs = randomNUMTs[['chr','randomStart','randomEnd']]
    randomNUMTs = pybedtools.BedTool.from_dataframe(randomNUMTs)
    overlap_freq = (randomNUMTs + df_ref).count() / numtsNo
    freq_list.append(overlap_freq)

df_pro = pd.DataFrame(freq_list, columns=["Percentage"])
Pvalue_less = 1 - (len(df_pro[df_pro['Percentage'] > numtsTargetPerc ]) / simuRuns)
Pvalue_greater = 1 - (len(df_pro[df_pro['Percentage'] < numtsTargetPerc ]) / simuRuns)

df_pro.columns = ["Pvalue_less(" + str(Pvalue_less) + ") & Pvalue_greater(" + str(Pvalue_greater) + ")"]

print("Pvalue(less):", Pvalue_less)
print("Pvalue(greater):", Pvalue_greater)

df_pro.to_csv(input1 + '.enrichmentOutput.' + output +'.tsv', index=True, sep = '\t')
