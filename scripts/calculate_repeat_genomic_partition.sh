#!/bin/sh

## calculate_repeat_genomic_partition.sh
## alex amlie-wolf
## overlaps SPAR intervals with repeat annotations


if [ $# == 3 ]; then
    WDIR=$1
    FPREFIX=$2
    ANNOTATIONF=$3
    #configFile=$4
    #source "${configFile}"
    
    ANNOT_F="${WDIR}/${FPREFIX}_annot.with_conservation.xls"
    UNANNOT_F="${WDIR}/${FPREFIX}_unannot.with_conservation.xls"

    ## convert the annotated locus file into bed and perform bedtools overlap
    tail -n +2 ${ANNOT_F} | awk -F$'\t' 'BEGIN{OFS=FS} {printf "%s\t%d\t%d\t%s", $1, $2, $3, $25; for(i=4; i<=NF; i++) {if(i!=6 && i!=25) {printf ",%s", $i}}; printf "\t0\t%s\n", $6}' | \
	${BEDTOOLS} intersect -s -sorted -a stdin -b ${ANNOTATIONF} -wao -f 0.51 > ${WDIR}/Genomewide_distribution_patterns_of_small_RNA_loci_repeat_annot_overlaps.txt
	#${BEDTOOLS} intersect -s -sorted -a stdin -b ${ANNOTATIONF} -wao > ${WDIR}/Genomewide_distribution_patterns_of_small_RNA_loci_repeat_annot_overlaps.txt

    ## same for the unannotated
    tail -n +2 ${UNANNOT_F} | awk -F$'\t' 'BEGIN{OFS=FS} {printf "%s\t%d\t%d\t%s", $1, $2, $3, $25; for(i=4; i<=NF; i++) {if(i!=6 && i!=25) {printf ",%s", $i}}; printf "\t0\t%s\n", $6}' | \
	${BEDTOOLS} intersect -s -sorted -a stdin -b ${ANNOTATIONF} -wao -f 0.51 > ${WDIR}/Genomewide_distribution_patterns_of_small_RNA_loci_repeat_unannot_overlaps.txt
	#${BEDTOOLS} intersect -s -sorted -a stdin -b ${ANNOTATIONF} -wao > ${WDIR}/Genomewide_distribution_patterns_of_small_RNA_loci_repeat_unannot_overlaps.txt

    ## call the R script to make the table with partition assignments
    ${RSCRIPT} `dirname $0`/generate_repeat_partition_tables.R ${FPREFIX} ${WDIR}
    
else
    echo "Usage: $0 <Working directory> <file name prefix> <annotation file>"
fi
