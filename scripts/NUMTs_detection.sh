#! /bin/bash

module load samblaster/0.1.24
module load samtools 
module load Sambamba/0.6.6


INPUT_BAM=$1
OUTPUT_DIR=$2
REF_GRCh38=$3

CLUSTER_SCRIPT='./searchNumtCluster_fromDiscordantReads.py'
BREAKPOINT_SCRIPT='./searchBreakpoint_fromblatoutputs.py'

SAMPLE_ID1="${INPUT_BAM##*/}"
SAMPLE_ID2=${SAMPLE_ID1%.bam}
OUTPUT="${OUTPUT_DIR}/${SAMPLE_ID2}"

samtools view -@ 16 -m 10G -h -F 2 $INPUT_BAM | grep -e @ -e MT -e chrM | samtools sort -@ 16 -m 10G -n  | samtools view -h | samblaster --ignoreUnmated -e -d ${OUTPUT}.mt.disc.sam -s ${OUTPUT_wgs}.split.sam -o /dev/null

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
    python ${BREAKPOINT_SCRIPT} ${OUTPUT}.psl ${sampleID} ${CHR} ${START} ${END}  ${OUTPUT}
    rm ${OUTPUT}.fasta
    rm ${OUTPUT}.sam
done

