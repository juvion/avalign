#!/bin/bash

AVALIGN=/media/expand/avalign
WORKDIR=${AVALIGN}/analysis
REF=${AVALIGN}/ref
TOOLS=${AVALIGN}/tools

export PATH=$PATH:$TOOLS


#--------------------------
#date: 6/2/2015
#work_1a: meta data extract
#--------------------------

function work_1a {
    subdir=$WORKDIR/work_1a
    cd $subdir
    echo "processing $subdir ..."
    metadata.py ${REF}/multiz100way/alignments/knownGene.exonNuc.100way.fa ${subdir}/knownGene.exonNuc.100way.metadata.out
}


#-------------------------------
#date: 6/3/2015
#work_1b: unique ucsc id extract
#-------------------------------

function work_1b {
    subdir=$WORKDIR/work_1b
    cd $subdir
    echo "processing $subdir ..."
    unique_ucsc_id.py ${REF}/multiz100way/alignments/knownGene.exonAA.100way.fa ${subdir}/knownGene.exonAA.100way_unique_ucsc_id.out
}


#-----------------------------------------
#date: 6/5/2015
#work_1c: concatenate the protein sequence
#-----------------------------------------

function work_1c {
    subdir=$WORKDIR/work_1c
    cd $subdir
    echo "processing $subdir ..."
    concatenate.py ${REF}/multiz100way/alignments/knownGene.exonAA.100way.fa ${subdir}/knownGene.exonAA.100way_concatenate.out
}

#-----------------------------
#choose to turn on which work
#-----------------------------
work_flow=( "work_1a:false"
            "work_1b:false"
            "work_1c:true" )


#--------------------------
#Excute the turned on work
#--------------------------
for work in "${work_flow[@]}"; do
    Key="${work%%:*}"
    Value="${work##*:}"
    if [[ "$Value" == true ]]; then
        $Key
    fi
done