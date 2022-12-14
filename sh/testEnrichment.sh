#!/usr/bin/bash
# dependencies: kycg, testFisher.R bedtools
# source ~/repo/KnowYourCG/sh/testEnrichment.sh
# testEnrichment /scr1/users/zhouw3/projects/20220609_ExpressionMethylationCorrelation/20220815_Clark/Output/EBcells/NegSig.5.CPGonly.bed /scr1/users/zhouw3/tmp/FeatureAggregationtest3 ~/references/mm10/annotation/cpg/cpg_nocontig.bed.gz /scr1/users/zhouw3/projects/20220609_ExpressionMethylationCorrelation/20220815_Clark/Output/EBcells/AllTestedCPG.bed /mnt/isilon/zhou_lab/projects/20191221_references/mm10/featuresHQ/

function set_environment {
    qry=$1
    ref=$2                      # assume sorted, three columns
    uni=$3                      # assume sorted, use na to skip
    fea=$4                      # assume sorted, five columns
    opt=$5
    TMPFDR=$(basename ${fea}).tmp
}

function bed2cg {
  local ref=$1
  local bed=$2
  bedtools intersect -a $ref -b $bed -sorted -c |
    cut -f4 | kycg3 pack -fa -
}

function testEnrichment() (     # this spawn a subshell
    set_environment $1 $2 $3 $4 $5
    >&2 echo
    >&2 echo "================="
    >&2 echo "Query:     $qry"
    >&2 echo "Ref (bed): $ref"
    >&2 echo "Universe:  $uni"
    >&2 echo "Feature:   $fea"
    >&2 echo "Option:    $opt"
    >&2 echo "Temp Dir:  $TMPFDR"
    >&2 echo "================="

    rm -rf $TMPFDR/*

    if [[ $qry != *.cg && $qry != *.cg.gz ]]; then
      mkdir -p $TMPFDR
      >&2 echo "Packing query to $TMPFDR/$(basename $qry).cg.gz"
      >&2 echo "To save time, please provide .cg.gz file (by using bed2cg)"
      bed2cg $ref $qry | gzip -c >$TMPFDR/$(basename $qry).cg.gz
      qry=$TMPFDR/$(basename $qry).cg.gz
    fi

    if [[ $uni == "na" || $ref == $uni ]]; then
      uni_opt=""
    else
      if [[ $uni != *.cg && $uni != *.cg.gz ]]; then
        mkdir -p $TMPFDR
        >&2 echo "Packing universe to $TMPFDR/$(basename $uni).cg"
        >&2 echo "To save time, please provide .cg file (converted using bed2cg)"
        bed2cg $ref $uni | gzip -c >$TMPFDR/$(basename $uni).cg.gz
        uni=$TMPFDR/$(basename $uni).cg
      fi
      uni_opt="-u $uni"
    fi

    >&2 echo "Testing overlaps..."
    if [[ -f $fea ]]; then
      kycg3 overlap $uni_opt $qry $fea | testFisher.R stdin
    elif [[ -d $fea ]]; then
      find $fea -type f | grep '.cg' | while read f; do
        kycg3 overlap $uni_opt $qry $f | awk -v f=$f '{print $0,f;}'
      done | testFisher.R stdin
    fi
    
    >&2 echo "All completed."
)
