#! /bin/bash

################################################################################
## This script detects methylation sites from long-read sequencing
## nanopolish need to be installed to run the pipeline
## nanopolish can be downloaded at https://github.com/jts/nanopolish
## Example of run: nanopolish_methylationDetection.sh inputDir outputDir ref numtRegion
################################################################################


INPUT_DIR=$1 #inputs directory 
OUTPUT_DIR=$2 #output directory
REF_GRCh38=$3 #reference genome file
regionNUMTnuclear=$4 #NUMT regions e.g chr5:32238477-32438477
SAMPLE_ID="${INPUT_DIR##*/}"

# merge all the fastq file into one file #
fastq_file="${OUTPUT_DIR}/input_fastq*/${SAMPLE_ID}.fastq.gz"
fast5_file="${INPUT_DIR}/fast5_pass/"

cat ${INPUT_DIR}/fastq_pass/*gz >${fastq_file}

# create an index file that links read ids with their signal-level data in the FAST5 files:
# ./nanopolish index -d ${fast5_file} ${fastq_file}

# Align the reads to reference genome
bam_file="${INPUT_DIR}/bam/mapped.bam"
methylationMT_file="${OUTPUT_DIR}/methylation/${SAMPLE_ID}_chrM.methylation.tsv"
methylationNUMTnuclear_file="${OUTPUT_DIR}/methylation/${SAMPLE_ID}_NumtNuclear.methylation.tsv"

regionMT="chrM:1-16569"

./nanopolish call-methylation -t 8 -r ${fastq_file} -b ${bam_file} -g ${REF_GRCh38} -w ${regionNUMTnuclear} > ${methylationNUMTnuclear_file}
./nanopolish call-methylation -t 8 -r ${fastq_file} -b ${bam_file} -g ${REF_GRCh38} -w ${regionMT} > ${methylationMT_file}

