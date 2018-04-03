set -e
#source `dirname $0`/../config.sh

if [ $# -lt 3 ]
then
  echo "USAGE: `basename $0` <input.bedgraph> <chrom.sizes> <out.bigWig"
  exit 1
fi

#chromInfo=`dirname $0`/../annot/chromInfo.txt
INBG=$1 # input bedgraph
chromInfo=$2
OUTBIGWIG=$3 #"${INBG%.*}.bigWig" # output bedgraph

#if [ -f "${INBG}" ] && [ -s "${INBG}" ]; then
if [ -s "${INBG}" ]; then
  ${BGTOBIGWIG} ${INBG} ${chromInfo} ${OUTBIGWIG}
fi
#bedGraphToBigWig
