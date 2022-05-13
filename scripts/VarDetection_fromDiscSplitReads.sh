#! /bin/bash

INPUT_LIST="fastaFiles.list" # list of fasta files 

OUTPUT_DIR='outputs'
OUTPUT_DIR_cap3="${OUTPUT_DIR}/assembly"
OUTPUTDIR_psl_human="${INPUT_DIR}_pslHuman"
OUTPUTDIR_psl_chimp="${INPUT_DIR}_pslChimp"
OUTPUTDIR_aln_human="${INPUT_DIR}_alnHuman"
OUTPUTDIR_aln_humanchimp="${INPUT_DIR}_alnHumanChimp"

filelines=`cat $INPUT_LIST`
for line in $filelines ; do
    echo "processing $line ..."
    contigINPUT="$(echo $line | cut -d, -f1)"
    sampleIndex="${contigINPUT##*/}"

    refHuman="./reference/chrM.fa"
    refChimp="./reference/chimpMT.fa"
    refHumanChimp="./reference/Homo_sapiens.GRCh38.Pan_troglodytes.Pan_tro_3.0.dna.chromosome.MT.fa"

    ### assembly ###
    cap3_out="${contigINPUT}.cap"
    ./CAP3/cap3 ${fasta} >${cap3_out}
    
    ### reference alignment - psl file ###
    psl_human="${OUTPUTDIR_psl_human}/${sampleID}_${clusterID}.human.psl"
    blat ${refHuman} ${cap3_out}.idContigs ${psl_human}
 
    ### multiple alignment ###
    cat ${refHuman} ${contigINPUT}  >${contigINPUT}.humanMT.fasta
    cat ${refHumanChimp} ${contigINPUT}  >${contigINPUT}.humanchimpMT.fasta

    alnHuman="${OUTPUTDIR_aln_human}/${sampleIndex}"
    alnHumanChimp="${OUTPUTDIR_aln_humanchimp}/${sampleIndex}"
    clustalo -i ${contigINPUT}.humanMT.fasta  -o ${alnHuman}.humanMTaln.fa --outfmt=fa --force
    clustalo -i ${contigINPUT}.humanchimpMT.fasta  -o ${alnHumanChimp}.humanchimpMTaln.fa --outfmt=fa --force

    ### variant detection ###
    numts_bed="./reference/numts.bed"
    pyHumanChimp="./generateVariantTable.HumanChimp.py"
    pyHuman="./generateVariantTable.Human.py"
    python ${pyHuman} ${alnHumanChimp}.humanMTaln.fa ${numts_bed}
    python ${pyHumanChimp} ${alnHumanChimp}.humanchimpMTaln.fa ${numts_bed}

done
