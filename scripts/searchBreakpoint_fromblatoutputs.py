#!/usr/bin/env python
################################################################################
## This script takes outputs from searchNumtCluster_fromDiscordantReads.py and NUMT_detection.sh 
## to look for NUMT breakpoints
################################################################################

import fileinput
import sys, os
import numpy as np
import pandas as pd
import glob
import scipy.stats as stats


def f_nu(row):
	mismatchLEN = 3 
	readLEN = 150
	if ((row['strand'] == "+") & (row['Qend'] >= (readLEN-mismatchLEN))):
		pointFrom = "nu_Tstart_Bright"
	elif ((row['strand'] == "+") & (row['Qstart'] <= mismatchLEN)):
		pointFrom = "nu_Tend_Bleft"
	elif row['strand'] == "-":
		pointFrom = "nu_NegStrand"
	else:
		pointFrom = "nu_useLess"

	return pointFrom

def f_mt(row):
	mismatchLEN = 3
	readLEN = 150
	if ((row['strand'] == "+") & (row['Qend'] >= (readLEN-mismatchLEN))):
		pointFrom = "mt_Tstart"	
	elif ((row['strand'] == "+") & (row['Qstart'] <= mismatchLEN)):
		pointFrom = "mt_Tend"
	elif ((row['strand'] == "-") & (row['Qend'] >= (readLEN-mismatchLEN))):
		pointFrom = "mt_Tend"
	elif ((row['strand'] == "-") & (row['Qstart'] <= mismatchLEN)):
		pointFrom = "mt_Tstart"
	else:
		pointFrom = "mt_useLess"
	return pointFrom

####################  blat psl   ############################
## START, END position, flake length 500bp in the input file

INPUT_PSL,SAMPLEID,CHR,START,END,OUTPUT=sys.argv[1:]

mismatchLEN = 3 # number of mismatched bases
readLEN = 150 # wgs read length
START = int(START)
END = int(END)
hg37CHR = CHR.replace("chr","") ## delete when align to hg38 only

df_input = pd.read_csv(INPUT_PSL, skiprows=5, sep="\t", names=["match","misMatch","repMatch","Ns","QgapCount","QgapBases","TgapCount","TgapBases","strand","Qname","Qsize","Qstart","Qend","Tname","Tsize","Tstart","Tend","blockCount","blockSizes","qStarts","tStarts"])
df_input['matchLEN'] = df_input['Tend'] - df_input['Tstart']

df1 = df_input[(df_input['matchLEN'] < (readLEN - 10)) & (df_input['misMatch'] <= mismatchLEN)] 
df2 = df1[(df1['Tend'] >= (readLEN - mismatchLEN)) | (df1['Tend'] <= mismatchLEN )]

df_mt = df2[df2['Tname'].isin(['chrM','MT'])]
df_nu = df2[ (df2['Tname'].isin([CHR, hg37CHR])) & (df2['Tstart'] > START) & (df2['Tend'] < END)]
df_m = pd.concat([df_nu, df_mt])

df_nu['pointGroup'] = df_nu.apply(f_nu, axis=1)
df_mt['pointGroup'] = df_mt.apply(f_mt, axis=1)
df_nu['chr'] = CHR
df_nu['Group'] =  df_nu['pointGroup'].str.replace(r'_T.*B','')
df_mt['chr'] = 'chrM'

## list all nuclear breakpoints ##
df_nuLeft = df_nu[df_nu['pointGroup'] == 'nu_Tend_Bleft']
df_nuLeftG = df_nuLeft.groupby(['pointGroup','Group','chr','Tend','strand']).size().reset_index(name="readsCount")
df_nuRight = df_nu[df_nu['pointGroup'] == 'nu_Tstart_Bright']
df_nuRightG = df_nuRight.groupby(['pointGroup','Group','chr','Tstart','strand']).size().reset_index(name="readsCount")
df_nuBothG = pd.concat([df_nuLeftG, df_nuRightG])

## list nuclear breakpoints also mapped to mtDNA ##
df_nu_mt = df_nu[df_nu['Qname'].isin(list(df_mt['Qname']))]
df_nuLeft_mt = df_nu_mt[df_nu_mt['pointGroup'] == 'nu_Tend_Bleft']
df_nuLeftG_mt = df_nuLeft_mt.groupby(['pointGroup','Group','chr','Tend','strand']).size().reset_index(name="readsCount")
df_nuRight_mt = df_nu_mt[df_nu_mt['pointGroup'] == 'nu_Tstart_Bright']
df_nuRightG_mt = df_nuRight_mt.groupby(['pointGroup','Group','chr','Tstart','strand']).size().reset_index(name="readsCount")
df_nuBothG_mt = pd.concat([df_nuLeftG_mt, df_nuRightG_mt])

### define all mtDNA breakpoints  ### 
df_mtTend = df_mt[df_mt['pointGroup'] == 'mt_Tend']
df_mtTendG = df_mtTend.groupby(['pointGroup','chr','Tend','strand']).size().reset_index(name="readsCount")
df_mtTstart = df_mt[df_mt['pointGroup'] == 'mt_Tstart']
df_mtTstartG = df_mtTstart.groupby(['pointGroup','chr','Tstart','strand']).size().reset_index(name="readsCount")
df_mtBothG = pd.concat([df_mtTendG, df_mtTstartG])
df_mtBothG['Group'] = 'UKn'

## list mtDNA breakpoints also mapped to nuclear ###
df_mtTendLeft = df_mtTend[df_mtTend['Qname'].isin(list(df_nuLeft['Qname']))]
df_mtTendRight = df_mtTend[df_mtTend['Qname'].isin(list(df_nuRight['Qname']))]
df_mtTstartLeft = df_mtTstart[df_mtTstart['Qname'].isin(list(df_nuLeft['Qname']))]
df_mtTstartRight = df_mtTstart[df_mtTstart['Qname'].isin(list(df_nuRight['Qname']))]

df_mtLeftTstartG = df_mtTstartLeft.groupby(['pointGroup','chr','Tstart','strand']).size().reset_index(name="readsCount")
df_mtLeftTstartG['Group'] = "mtLeft" ## "mtLeft": left breakpoint of mtDNA
df_mtLeftTendG = df_mtTendLeft.groupby(['pointGroup','chr','Tend','strand']).size().reset_index(name="readsCount")
df_mtLeftTendG['Group'] = "mtLeft"
df_mtRightTstartG = df_mtTstartRight.groupby(['pointGroup','chr','Tstart','strand']).size().reset_index(name="readsCount")
df_mtRightTstartG['Group'] = "mtRight" ## "mtRight": right breakpoint of mtDNA
df_mtRightTendG = df_mtTendRight.groupby(['pointGroup','chr','Tend','strand']).size().reset_index(name="readsCount")
df_mtRightTendG['Group'] = "mtRight"
df_mtConfG = pd.concat([df_mtLeftTstartG,df_mtLeftTendG,df_mtRightTstartG,df_mtRightTendG])

## output all possible nuclear and mtDNA breakpoints ##
df_numtAllG = pd.concat([df_nuBothG, df_mtBothG]) 
#df_numtAllG['clusterID'] = CLUSTERID
df_numtAllG['sampleID'] = SAMPLEID
df_numtAllG['Tstart'] = df_numtAllG['Tstart'].fillna(-1).astype(int)
df_numtAllG['Tend'] = df_numtAllG['Tend'].fillna(-1).astype(int)

## output confident breakpoints ##
df_numtConfG = pd.concat([df_nuBothG_mt, df_mtConfG])
#df_numtConfG['clusterID'] = CLUSTERID
df_numtConfG['sampleID'] = SAMPLEID
df_numtConfG['Tstart'] = df_numtConfG['Tstart'].fillna(-1).astype(int)
df_numtConfG['Tend'] = df_numtConfG['Tend'].fillna(-1).astype(int)

## write to file ##
df_numtConfG.to_csv(OUTPUT + '.Breakpoints.tsv', sep = '\t',header=False, index=False)

