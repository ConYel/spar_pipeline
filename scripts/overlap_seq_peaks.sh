
#/usr/local/bin/bedtools intersect -a scripts/testout/adipose1_star_hc_sorted.bam.collapsed.bed.seq -b ~/SPAR_out/adipose1_star_hc_sorted.bam.annot.final.csv -f 0.8 -s -wo

SEQFILE=$1 # collapsed BED with sequences
PEAKFILE=$2 # called peaks (annotated or unannotated)
minOverlap=${3:-0.9}
if [ $# -lt 2 ]; then
  echo "USAGE: $0 <collapsed.BED.with.seqs> <peaks.BED> [minOverlapPct]"
  exit 1
fi

# output overlaps only
echo -e "#chr\tchrStart\tchrEnd\tLength\tExpression\tStrand\tSequence\t`head -n 1 ${PEAKFILE} | awk '{print substr($0,2)}'`"
${BEDTOOLS} intersect -a ${SEQFILE} -b ${PEAKFILE} -f ${minOverlap} -s -wa -wb | sort -k1,1 -k2,2n -k3,3n -k6,6

