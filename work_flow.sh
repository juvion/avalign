#!/bin/bash

#define directories
# AVALIGN=/media/expand/avalign
PWD=`pwd`
AVALIGN=${PWD%\tools}
WORKDIR=${AVALIGN}/analysis
REF=${AVALIGN}/ref
TOOLS=${AVALIGN}/tools

DATE=`date +%Y-%m-%d`

export PATH=$PATH:$TOOLS

#--------------------------
#date: 6/2/2015
#work_1a: meta data extract
#--------------------------

function work_1a {
    subdir=$WORKDIR/work_1a
    input=${REF}/multiz100way/alignments/knownGene.exonNuc.100way.fa
    output=${subdir}/knownGene.exonNuc.100way.metadata.out
    error_log=${subdir}/error_${DATE}.log

    if [ ! -d  $subdir ]; then
        mkdir $subdir
    fi
    if [ ! -d ${subdir}/old ]; then
        mkdir ${subdir}/old
    fi

    cd $subdir
    echo "processing $subdir ..."
    if [ -f $output ]; then 
        mv $output ${subdir}/old
    fi

    if [ -f $error_log ]; then 
        mv $error_log ${subdir}/old
    fi
    if [ ! -f $input ]; then 
        echo "Error: input file $input is missing"
    else
        python $TOOLS/metadata.py $input $output > $error_log
    fi

}


#-------------------------------
#date: 6/3/2015
#work_1b: unique ucsc id extract
#-------------------------------

function work_1b {
    subdir=$WORKDIR/work_1b
    input=${REF}/multiz100way/alignments/knownGene.exonAA.100way.fa
    output=${subdir}/knownGene.exonAA.100way_unique_ucsc_id.out
    error_log=${subdir}/error_${DATE}.log

    if [ ! -d  $subdir ]; then
        mkdir $subdir
    fi
    if [ ! -d ${subdir}/old ]; then
        mkdir ${subdir}/old
    fi

    cd $subdir
    echo "processing $subdir ..."
    if [ -f $output ]; then 
        mv $output ${subdir}/old
    fi

    if [ -f $error_log ]; then 
        mv $error_log ${subdir}/old
    fi

    if [ ! -f $input ]; then 
        echo "Error: input file $input is missing"
    else
        python $TOOLS/unique_ucsc_id.py $input $output > $error_log
    fi
}


#-----------------------------------------
#date: 6/5/2015
#work_1c: concatenate the protein sequence
#-----------------------------------------

function work_1c {
    subdir=$WORKDIR/work_1c
    # input=${REF}/multiz100way/alignments/knownGene.exonAA.100way.fa
    # output=${subdir}/knownGene.exonAA.100way_concatenate.out
    input=${REF}/multiz100way/alignments/refGene.exonAA.100way.fa
    output=${subdir}/refGene.exonAA.100way_concatenate.out
    error_log=${subdir}/error_${DATE}.log

    if [ ! -d  $subdir ]; then
        mkdir $subdir
    fi
    if [ ! -d ${subdir}/old ]; then
        mkdir ${subdir}/old
    fi

    cd $subdir
    echo "processing $subdir ..."

    if [ -f $output ]; then 
        mv $output ${subdir}/old
    fi

    if [ -f $error_log ]; then 
        mv $error_log ${subdir}/old
    fi


    if [ ! -f $input ]; then 
        echo "Error: input file $input is missing"
    else
        python $TOOLS/concatenate.py $input $output > $error_log
    fi

}


#-----------------------------------------
#date: 6/10/2015
#work_1d: unique refseq id extract
#-----------------------------------------

function work_1d {
    subdir=$WORKDIR/work_1d
    input=${REF}/multiz100way/alignments/refGene.exonAA.100way.fa
    output=${subdir}/refGene.exonAA.100way_unique_refseq_id.out
    error_log=${subdir}/error_${DATE}.log

    if [ ! -d  $subdir ]; then
        mkdir $subdir
    fi
    if [ ! -d ${subdir}/old ]; then
        mkdir ${subdir}/old
    fi

    cd $subdir
    echo "processing $subdir ..."
    if [ -f $output ]; then 
        mv $output ${subdir}/old
    fi

    if [ -f $error_log ]; then 
        mv $error_log ${subdir}/old
    fi

    if [ ! -f $input ]; then 
        echo "Error: input file $input is missing"
    else
        python $TOOLS/unique_refseq_id.py $input $output > $error_log
    fi
}


#-----------------------------------------
#date: 6/10/2015
#work_1e: match cds ranges
#-----------------------------------------

function work_1e {
    subdir=$WORKDIR/work_1e
    input1=$WORKDIR/work_1d/knownGene.exonAA.100way
    input2=$WORKDIR/work_1b/knownGene_match_ucsc_id-refseq-gene_symbl-cds_range.out
    output=${subdir}/knownGene.exonAA.100way_ucsc_id_match_refGene.out
    error_log=${subdir}/error_${DATE}.log

    if [ ! -d  $subdir ]; then
        mkdir $subdir
    fi
    if [ ! -d ${subdir}/old ]; then
        mkdir ${subdir}/old
    fi

    cd $subdir
    echo "processing $subdir ..."
    if [ -f $output ]; then 
        mv $output ${subdir}/old
    fi

    if [ -f $error_log ]; then 
        mv $error_log ${subdir}/old
    fi

    if [ ! -f $input1 ] || [ ! -f $input2 ]; then
       echo "Error: input files $input1 or $input2 are missing."
    else
        python $TOOLS/match_cds_ranges.py $input1 $input2 $output > $error_log
    fi
}


#---------------------------------------------------
#date: 6/16/2015
#work_1f: check the duplicity of matched ucsc ids.
#---------------------------------------------------

function work_1f {
    subdir=$WORKDIR/work_1f
    input1=${REF}/multiz100way/alignments/knownGene.exonAA.100way.fa
    input2=$WORKDIR/work_1b/knownGene_match_ucsc_id-refseq-gene_symbl-cds_range.out
    error_log=${subdir}/error_${DATE}.log

    if [ ! -d  $subdir ]; then
        mkdir $subdir
    fi
    if [ ! -d ${subdir}/old ]; then
        mkdir ${subdir}/old
    fi

    cd $subdir
    echo "processing $subdir ..."
    if [ -f $output ]; then 
        mv $output ${subdir}/old
    fi

    if [ -f $error_log ]; then 
        mv $error_log ${subdir}/old
    fi
    
    if [ ! -f $input1 ] || [ ! -f $input2 ]; then
       echo "Error: input files $input1 or $input2 are missing."
    else
        python $TOOLS/check_duplicate.py $input1 $input2 > $error_log
    fi
}


#---------------------------------------------------
#date: 6/19/2015
#work_1g: get the range of genes from AVA dumped data.
#---------------------------------------------------

function work_1g {
    subdir=$WORKDIR/work_1g
    input=${REF}/hs_transcript_info.csv
    output=${subdir}/AVA_gene_ranges.in
    error_log=${subdir}/error_${DATE}.log

    if [ ! -d  $subdir ]; then
        mkdir $subdir
    fi
    if [ ! -d ${subdir}/old ]; then
        mkdir ${subdir}/old
    fi

    cd $subdir
    echo "processing $subdir ..."
    if [ -f $output ]; then 
        mv $output ${subdir}/old
    fi

    if [ -f $error_log ]; then 
        mv $error_log ${subdir}/old
    fi
    
    if [ ! -f $input ]; then
       echo "Error: input files $input is missing."
    else
        python $TOOLS/get_range.py $input $output > $error_log
    fi
}



#---------------------------------------------------
#date: 6/30/2015
#work_2a: nucleotides alignment extraction
#---------------------------------------------------
function work_2a {
    subdir=$WORKDIR/work_2a
    input=${REF}/AVA_genes.bed
    error_log=${subdir}/error_${DATE}.log

    if [ ! -d  $subdir ]; then
        mkdir $subdir
    fi
    if [ ! -d ${subdir}/old ]; then
        mkdir ${subdir}/old
    fi

    cd $subdir
    echo "processing $subdir ..."

    if [ -f $error_log ]; then 
        mv $error_log ${subdir}/old
    fi
    
    if [ ! -f $input ]; then
       echo "Error: input files $input is missing."
    else
        python $TOOLS/maf_iterate.py $input $subdir 1 1 25 > $error_log
    fi
}



#-----------------------------
#choose to turn on which work
#-----------------------------
work_flow=( "work_1a:false"
            "work_1b:false"
            "work_1c:false"
            "work_1d:false"
            "work_1e:false" 
            "work_1f:false"
            "work_1g:false"
            "work_2a:true")



#--------------------------
#Excute the switched on work
#--------------------------
for work in "${work_flow[@]}"; do
    echo $work
    Key="${work%%:*}"
    Value="${work##*:}"
    if [[ "$Value" == true ]]; then
        $Key
    fi
done