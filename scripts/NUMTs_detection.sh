#! /bin/bash

################################################################################
## This script detects NUMTs from whole genome sequencing BAM files
## samtools, samblaster and blat need to be installed to run the pipeline
## samtools can be downloaded at http://www.htslib.org/download/
## samblaster can be downloaded at https://github.com/GregoryFaust/samblaster
## blat can be downloaded at http://hgdownload.soe.ucsc.edu/admin/exe/
################################################################################

## load modules on HPC
module load samblaster/0.1.24
module load samtools 
module load Sambamba/0.6.6


INPUT_BAM=$1 # input WGS bam file
OUTPUT_DIR=$2 # output folder path
REF_GRCh38=$3 # human reference genome

CLUSTER_SCRIPT='./searchNumtCluster_fromDiscordantReads.py'
BREAKPOINT_SCRIPT='./searchBreakpoint_fromblatoutputs.py'

SAMPLE_ID1="${INPUT_BAM##*/}"
SAMPLE_ID2=${SAMPLE_ID1%.bam}
OUTPUT="${OUTPUT_DIR}/${SAMPLE_ID2}"

INPUT_DISC="${OUTPUT}.mt.disc.sam"
INPUT_SPLIT="${OUTPUT_wgs}.mt.split.sam"

samtools view -@ 16 -m 10G -h -F 2 $INPUT_BAM | grep -e @ -e MT -e chrM | samtools sort -@ 16 -m 10G -n  | samtools view -h | samblaster --ignoreUnmated -e -d $INPUT_DISC -s $INPUT_SPLIT -o /dev/null

python $CLUSTER_SCRIPT ${SAMPLE_ID2} ${INPUT_BAM} ${INPUT_DISC}

echo "look for cluster's done, move to look for breakpoints"

filelines=`cat ${INPUT_DISC}.breakpointINPUT.tsv`

for line in $filelines ; do
    echo $line
    INPUT_Dis="$(echo $line | cut -d, -f3)"
    INPUT_Split="$(echo $line | cut -d, -f4)"
    INPUT_WGS="$(echo $line | cut -f5)"
    CHR="$(echo $line | cut -f6)"
    START="$(echo $line | cut -f7)"
    END="$(echo $line | cut -f8)"
    sampleID="$(echo $line | cut -f1)"
    REGION="${CHR}:${START}-${END}"
    OUTPUT="${OUTPUT_DIR}/${sampleID}_${CHR}.${START}.${END}"
    samtools view ${INPUT_WGS} ${REGION} >${OUTPUT}.sam

    awk '$6 !~ /150M|149M|148M|149S|148S/' ${OUTPUT}.sam | cut -f1,10 >${OUTPUT}.fasta
    perl -pi -e 's/^/>/g' ${OUTPUT}.fasta
    perl -pi -e 's/\t/\n/g' ${OUTPUT}.fasta
    blat ${REF_GRCh38}  ${OUTPUT}.fasta  ${OUTPUT}.psl
    python ${BREAKPOINT_SCRIPT} ${OUTPUT}.psl ${sampleID} ${CHR} ${START} ${END} ${OUTPUT}
    rm ${OUTPUT}.fasta
    rm ${OUTPUT}.sam
done

