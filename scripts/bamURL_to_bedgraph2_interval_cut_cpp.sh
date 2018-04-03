set -e

BAMURL=$1 # BAM URL
OUTPREFIX=$2 # output directory + prefix
minReadLength=${3:-15}
maxReadLength=${4:-44}
keep5pClipped=${5:-""}

BAM2BEDGRAPH=`dirname $0`/bam2bedgraph_interval

OUTDIR=${OUTPREFIX##*/}
if [ ! -z ${OUTDIR} ]; then
  mkdir -p ${OUTDIR}
fi

outPosBedgraph="${OUTPREFIX}.pos.bedgraph";
outNegBedgraph="${OUTPREFIX}.neg.bedgraph";
outCollapsedBED="${OUTPREFIX}.collapsed.bed";
outHeader="${OUTPREFIX}.chromInfo.txt";
>${outPosBedgraph}
>${outNegBedgraph}
printT "Loading BAM ${BAMURL}"
#curl -s -L ${BAMURL} | ${BAM2BEDGRAPH} - ${minReadLength} ${maxReadLength} | \
curl -s -L ${BAMURL} | ${BAM2BEDGRAPH} - ${minReadLength} ${maxReadLength} ${keep5pClipped} | \
  awk 'BEGIN{ posFile = "'${outPosBedgraph}'";
              negFile = "'${outNegBedgraph}'";
              headerFile = "'${outHeader}'";
              intervalFile = "'${outCollapsedBED}'";
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
           #if (chr!=chr_prev) print chr;
           #chr_prev = chr; 
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

printT "Loaded BAM successfuly." 
#printT "BEDGRAPH finished successfuly."
