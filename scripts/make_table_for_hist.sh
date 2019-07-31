#!/bin/bash
# File: make_table_for_hist.sh
## Konstantinos Con_Yel
## make a table for histogram from soft clipped alignment file Aligned.out.filtered.sorted.bam 

INPUT_F="$1"
nthreads="${2:-4}"

if [ $# -lt 1 ]; then
  echo "USAGE: <INPUT_FOLDER> folder with SPAR results from different samples <nthreads>"
  exit 1
fi

 
for filt_bam in $INPUT_F/*/mapping/Aligned.out.filtered.sorted.bam; do
[ -f "$filt_bam" ] || continue
  filename="${filt_bam%%/mapping*}" #remove trailing path
 echo "samtools view -@ $nthreads  "$filt_bam"   | awk '{print length($10)}'  | sort | uniq -c > ${filename}/results/histogram_table.txt"
 samtools view -@ $nthreads  "$filt_bam"   | awk '{print length($10)}'  | sort | uniq -c > ${filename}/results/histogram_table.txt
  done
