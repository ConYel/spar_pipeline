#!/bin/sh

## combine_partition_results.sh
## alex amlie-wolf
## combine partition information across mRNA, lncRNA, repeats

if [ $# == 2 ]; then
    WDIR=$1
    FPREFIX=$2

    ## call the R script to make the table with partition assignments
    ${RSCRIPT} `dirname $0`/combine_partition_results.R ${FPREFIX} ${WDIR}

    head -1 ${WDIR}/${FPREFIX}_all_partition_summary_per_annot_locus.txt > ${WDIR}/${FPREFIX}_all_partition_summary_per_annot_locus.txt.sorted
    tail -n +2 ${WDIR}/${FPREFIX}_all_partition_summary_per_annot_locus.txt | sort -k1,1 -k2,2n -k3,3n -k6,6 >> ${WDIR}/${FPREFIX}_all_partition_summary_per_annot_locus.txt.sorted
    mv ${WDIR}/${FPREFIX}_all_partition_summary_per_annot_locus.txt.sorted ${WDIR}/${FPREFIX}_all_partition_summary_per_annot_locus.txt

    head -1 ${WDIR}/${FPREFIX}_all_partition_summary_per_unannot_locus.txt > ${WDIR}/${FPREFIX}_all_partition_summary_per_unannot_locus.txt.sorted
    tail -n +2 ${WDIR}/${FPREFIX}_all_partition_summary_per_unannot_locus.txt | sort -k1,1 -k2,2n -k3,3n -k6,6 >> ${WDIR}/${FPREFIX}_all_partition_summary_per_unannot_locus.txt.sorted
    mv ${WDIR}/${FPREFIX}_all_partition_summary_per_unannot_locus.txt.sorted ${WDIR}/${FPREFIX}_all_partition_summary_per_unannot_locus.txt
    
else
    echo "Usage: $0 <Working directory> <file name prefix>"
fi
