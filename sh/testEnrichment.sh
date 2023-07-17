#!/usr/bin/bash
# dependencies: kycg, testFisher.R bedtools
## one should call this script directly by
## ./testEnrichment.sh /scr1/users/zhouw3/projects/20220609_ExpressionMethylationCorrelation/20220815_Clark/Output/EBcells/NegSig.5.CPGonly.bed mm10

## for full control, we need to source this file
## source ~/repo/KnowYourCG/sh/testEnrichment.sh
## testEnrichment0 /scr1/users/zhouw3/projects/20220609_ExpressionMethylationCorrelation/20220815_Clark/Output/EBcells/NegSig.5.CPGonly.bed ~/references/mm10/annotation/cpg/cpg_nocontig.bed.gz ~/references/mm10/annotation/cpg/cpg_nocontig.bed.gz /mnt/isilon/zhou_lab/projects/20191221_references/mm10/featuresHQ/ >out

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
    cut -f4 | kycg pack -fa -
}

function testEnrichment0() (     # this spawn a subshell
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

    if [[ $qry != *.cg ]]; then
      mkdir -p $TMPFDR
      >&2 echo "Packing query to $TMPFDR/$(basename $qry).cg"
      >&2 echo "To save time, please provide .cg file (by using bed2cg)"
      bed2cg $ref $qry >$TMPFDR/$(basename $qry).cg
      qry=$TMPFDR/$(basename $qry).cg
    fi

    if [[ $uni == "na" || $ref == $uni ]]; then
      uni_opt=""
    else
      if [[ $uni != *.cg ]]; then
        mkdir -p $TMPFDR
        >&2 echo "Packing universe to $TMPFDR/$(basename $uni).cg"
        >&2 echo "To save time, please provide .cg file (converted using bed2cg)"
        bed2cg $ref $uni >$TMPFDR/$(basename $uni).cg
        uni=$TMPFDR/$(basename $uni).cg
      fi
      uni_opt="-u $uni"
    fi

    >&2 echo "Testing overlaps..."
    if [[ -f $fea ]]; then
      kycg overlap $uni_opt $qry $fea | testFisher.R stdin
    elif [[ -d $fea ]]; then
      find $fea -type f | grep '.cg' | while read f; do
        kycg overlap $uni_opt $qry $f | awk -v f=$f '{print $0,f;}'
      done | testFisher.R stdin
    fi
    
    >&2 echo "All completed."
)

if [[ "${BASH_SOURCE[0]}" -ef "$0" ]]; then
    qry=$1
    genome=$2
    testEnrichment0 $qry ~/references/$genome/annotation/cpg/cpg_nocontig.bed.gz na ~/references/$genome/featuresHQ/
fi
