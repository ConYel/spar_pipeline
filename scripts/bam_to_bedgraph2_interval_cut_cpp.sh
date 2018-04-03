set -e

INBAM=$1

if [ ! -s ${INBAM} ]; then
  echo "***ERROR: ${INBAM} does not exist or empty"
  exit 1
fi

BAM=`basename ${INBAM}`
BAMDIR=`dirname ${INBAM}`
BAM2BEDGRAPH=`dirname $0`/bam2bedgraph_interval

# directory where output will be saved
# directory is either provided as input argument or is the same as the input BAM directory
OUTDIR=${2:-${BAMDIR}}
mkdir -p ${OUTDIR}

minReadLength=${3:-15}
maxReadLength=${4:-44}
keep5pClipped=${5:-""}

outPosBedgraph="${OUTDIR}/${BAM}.pos.bedgraph";
outNegBedgraph="${OUTDIR}/${BAM}.neg.bedgraph";
outCollapsedBED="${OUTDIR}/${BAM}.collapsed.bed";
outHeader="${OUTDIR}/${BAM}.chromInfo.txt";
>${outPosBedgraph}
>${outNegBedgraph}
${BAM2BEDGRAPH} ${INBAM} ${minReadLength} ${maxReadLength} ${keep5pClipped} | \
  awk 'BEGIN{ posFile = "'${outPosBedgraph}'";
              negFile = "'${outNegBedgraph}'";
              intervalFile = "'${outCollapsedBED}'";
              headerFile = "'${outHeader}'";
              OFS="\t";
              cmdPos=("sort -k1,1 -k2,2n -k3,3n > " posFile);
              cmdNeg=("sort -k1,1 -k2,2n -k3,3n > " negFile);
              minCov=1;
            }
       {
        if ($1=="chrInfo")
        { 
          chrName=$2; chrLength=$3;
          print chrName, chrLength > headerFile;
          next;
        }
 
        if ($1=="interval")
        {
           # 0-based half-open UCSC format
           # start=0-based; end=1-based
           chr=$2; chrStart=$3; len=$4; chrEnd=chrStart+len;
           expr = $5; strand=$6;
           print chr, chrStart, chrEnd, len, expr, strand > intervalFile;
           next;
        }
        if ($4<minCov) next; 
        if ($6=="+")
          print $1,$2,$3,$4 | cmdPos;
        else
          print $1,$2,$3,$4 | cmdNeg;
       }'

# need to sort because of unnatural order of samtools view 
# sort by chr, chrStart, chrEnd

#if [ -s "${outPosBedgraph}" ]; then
#  sort -k1,1 -k2,2n -k3,3n ${outPosBedgraph} -o ${outPosBedgraph}
#else
#  >"${outPosBedgraph}"
#fi

#if [ -s "${outNegBedgraph}" ]; then
#  sort -k1,1 -k2,2n -k3,3n ${outNegBedgraph} -o ${outNegBedgraph}
#else
#  >"${outNegBedgraph}"
#fi

#printT "BEDGRAPH finished successfuly."
