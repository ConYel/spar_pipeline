set -e

bgfile=$1
chromInfo=$2
outbigwig=$3

if [ $# -lt 3 ]; then
  echo "USAGE: ${0} <input.bedgraph> <chrom.sizes> <out.bigWig>"
  exit 1
fi

#source `dirname $0`/../config.sh
#outbigwig=${bgfile/bedgraph/bigWig}
#chromInfo=`dirname $0`/../annot/chromInfo.txt
# check if bedgraph file exists and is not zero size
if [ -s "${bgfile}" ]; then
  ${BGTOBIGWIG} ${bgfile} ${chromInfo} ${outbigwig}
else
  echo "ERROR: ${bgfile} not found or empty!"
  exit 1
fi
