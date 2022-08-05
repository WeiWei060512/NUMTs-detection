#!/usr/bin/env python

################################################################################
## This script takes outputs from NUMT_detection.sh to generates NUMT clusters
################################################################################

import fileinput
import sys, os
import pandas as pd
import glob
import scipy.stats as stats
import numpy as np

#################### Extract cluster from mtDNA discordant sam files  #############################

def cluster(data, maxgap):
	data.sort()
	groups = [[data[0]]]
	for x in data[1:]:
		if abs(x - groups[-1][-1]) <= maxgap:
            		groups[-1].append(x)
		else:
			groups.append([x])
	return groups


sampleID,wgsBAM,input1 = sys.argv[1:]


sampleID = sampleID.replace(".mt.disc","")

df0 = pd.read_csv(input1, names =['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','SM','RG','NM','BC','OC','ZX','ZY','SA'], sep="\t", engine='python', comment="@")

df = df0[~df0.RNAME.str.contains("Un_|random|\.")]
df1 = df[(df["MAPQ"].astype(int) > 0) & ((df['RNAME'] == "MT") | (df['RNAME'] == "chrM") | (df['RNEXT'] == "MT") | (df['RNEXT'] == "chrM"))]

##### remove reads map to mtDNA #####
df3 = df1[~df1.RNAME.isin(['chrM','MT'])]

##### order by chromosome and pos #####
df3.sort_values(['RNAME','POS'], inplace=True)

##### look for the clutser by mapgap on pos #####
##### extract the clusters with number of elements no less than 5 #####

output1 = pd.DataFrame([])

df_chr = df3.groupby(['RNAME'])
for clusterID, myclusters in df_chr:
    myclusters['POS'] = myclusters['POS'].astype(int)
    sub_cluster = cluster(myclusters['POS'].tolist(), maxgap=500)
    output = pd.DataFrame([])
    for x in sub_cluster:
        if len(x)>=2:
            mycluster = x
            df_cluster = df3[df3.POS.isin(mycluster)]
            df_cluster_pairMT = df[df.QNAME.isin(df_cluster['QNAME'])]
            mt_cluster = cluster(df_cluster['PNEXT'].tolist(), maxgap=500)
            for y in mt_cluster:
                df_cluster_pairMT_out = df_cluster_pairMT[df_cluster_pairMT.PNEXT.isin(y)]
                df_cluster_pairMT_out['subCluster_No'] = len(y)
                df_cluster_pairMT_out['Cluster_No'] = len(x)			
                df_cluster_pairMT_out['Cluster_ID'] = clusterID + "_" + str(min(df_cluster['POS'])) + "_" + \
                                                      str(max(df_cluster['POS'])) + "_MTboth_"  + \
                                                      str(min(df_cluster_pairMT_out['PNEXT'])) + "_" + \
                                                      str(max(df_cluster_pairMT_out['PNEXT']))
		output = output.append(df_cluster_pairMT_out, ignore_index=True)

            output1 = output1.append(output)
	    #output2 = df[df.QNAME.isin(output1['QNAME'])]
            output1["IndividualID"] = sampleID
            if len(output1) != 0:
                output1['clusterID'] = output1['RNAME'].astype(str) + "_" + output1['Cluster_No'].astype(str)
		output2 = output1.groupby(['IndividualID','Cluster_ID','Cluster_No','subCluster_No']).size().to_frame('size').reset_index()
            else: 
                next
#output2 = pd.DataFrame(output2)
#output2.columns=['IndividualID','Cluster_ID','Cluster_No','subCluster_No','size']
output1.to_csv(input1 + '.cluster.tsv', sep = '\t',header=True, index=False)
output2.to_csv(input1 + '.cluster.summary.tsv', sep = '\t',header=False, index=True)

##### generate cluster range for defining the breakpoints ######

#output2['Cluster_ID'].replace("_MT","",inplace=True)
#output2 = output2.drop_duplicates(["Cluster_ID"])
output2['disFile'] = input1
output2['splitFile'] = input1.replace('disc','split')
output2['wgsBAM'] = wgsBAM
df_pos = pd.DataFrame(output2['Cluster_ID'].str.split('_').tolist(),columns = ['chr','start','end','chrM','mtstart','mtend'])
del output2['Cluster_ID']
del output2['subCluster_No']
del output2['size']
output3 = pd.concat([output2, df_pos[['chr','start','end']]], axis=1)
output3 = output3.drop_duplicates(['chr','start','end'])
output3['start'] = output3['start'].astype(int) - 500
output3['end'] = output3['end'].astype(int) + 500 + 150
output3['Cluster_No'] = output3['Cluster_No'].astype(int)

##### output to file #####
output3.to_csv(input1 + '.breakpointINPUT.tsv', sep = '\t',header=False, index=False)

