#!/usr/bin/bash
# dependencies: kycg, testFisher.R bedtools
# source ~/repo/KnowYourCG/sh/testEnrichment.sh
# testEnrichment ~/references/mm10/annotation/cpg/cpg_nocontig.bed.gz /scr1/users/zhouw3/projects/20220609_ExpressionMethylationCorrelation/20220815_Clark/Output/EBcells/AllTestedCPG.bed /mnt/isilon/zhou_lab/projects/20191221_references/mm10/featuresHQ/ /scr1/users/zhouw3/projects/20220609_ExpressionMethylationCorrelation/20220815_Clark/Output/EBcells/NegSig.5.CPGonly.bed /scr1/users/zhouw3/tmp/FeatureAggregationtest3

function set_environment {
    ref=$1                      # assume sorted, three columns
    uni=$2                      # assume sorted, use na to skip
    fea=$3                      # assume sorted, five columns
    qry=$4
    out=$5
    opt=$6
    TMPFDR=${out}.$(basename ${fea}).tmp
}

function bed2cg {
  local ref=$1
  local bed=$2
  local out=$3
  bedtools intersect -a $ref -b $bed -sorted -c |
    cut -f4 | kycg pack - $out
}

function testEnrichment() (     # this spawn a subshell
    set_environment $1 $2 $3 $4 $5 $6
    echo
    echo "================="
    echo "Ref (bed): $ref"
    echo "Universe:  $uni"
    echo "Feature:   $fea"
    echo "Query:     $qry"
    echo "Output:    $out"
    echo "Option:    $opt"
    echo "Temp Dir:  $TMPFDR"
    echo "================="

    mkdir -p $TMPFDR
    rm -rf $TMPFDR/*
    mkdir -p $(dirname $out)

    if [[ $qry != *.cg ]]; then
      echo "Packing query to $TMPFDR/$(basename $qry).cg"
      echo "To save time, please provide .cg file (converted using bed2cg)"
      bed2cg $ref $qry $TMPFDR/$(basename $qry).cg
      qry=$TMPFDR/$(basename $qry).cg
    fi

    if [[ $uni == "na" || $ref == $uni ]]; then
      uni_opt=""
    else
      if [[ $uni != *.cg ]]; then
        echo "Packing universe to $TMPFDR/$(basename $uni).cg"
        echo "To save time, please provide .cg file (converted using bed2cg)"
        bed2cg $ref $uni $TMPFDR/$(basename $uni).cg
        uni=$TMPFDR/$(basename $uni).cg
      fi
      uni_opt="-u $uni"
    fi

    echo "Testing overlaps..."
    if [[ -f $fea ]]; then
      kycg overlap $uni_opt {} $f | testFisher.R stdin >$out
    elif [[ -d $fea ]]; then
      find $fea -name '*.cg' | while read f; do
        kycg overlap $uni_opt $qry $f | awk -v f=$f '{print $0,f;}'
      done | testFisher.R stdin >$out
    fi
    
    echo "All completed."
)
