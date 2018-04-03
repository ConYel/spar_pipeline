set -e

#source `dirname $0`/../config.sh

if [ $# -lt 3 ]
then
  echo "USAGE: `basename $0` .segm"
  exit 1
fi



INSEGM=$1 # input segmentation file
chromInfo=$2 #`dirname $0`/../annot/chromInfo.txt
OUTBIGBED=$3
OUTBED=${OUTBIGBED%.*}.bed

#echo "Processing ${fsegm}"
#OUTBED=${INSEGM%.*}.bed
#OUTBED="${INSEGM}.bed"
#OUTBIGBED="${OUTBED%.*}.bigBed"
# prepare BED file
awk 'BEGIN{OFS="\t"; dashrBaseURL="http://lisanwanglab.org/DASHR/get_bg_smrna.php?"}
         { 
           segmChr=$1;
           segmStart=$2;
           segmEnd=$3;
           segmStrand=$6;
           segmID=$4; segmExpr=$5;
           segmMaxProb5p=$13; segmMaxProb3p=$14;
           segmMaxPos5p=$15; segmMaxPos3p=$16;
           dispStrand="pos"; if (segmStrand=="-") dispStrand="neg";
           dashrLink=("<a href=\"" dashrBaseURL "coord=" segmChr ":" (segmStart+1) "-" segmEnd ":" dispStrand "\" target=\"_blank\">" segmID "</a>")
           print segmChr, segmStart, segmEnd, segmID, 1, segmStrand, dashrLink, segmExpr, segmMaxPos5p, segmMaxPos3p, segmMaxProb5p, segmMaxProb3p
         }' ${INSEGM} > ${OUTBED}

# generate bigBed track
#bedToBigBed -type=bed6+6 -as=bedPeakFields.as -tab ${OUTBED} ${chromInfo} ${OUTBED} 
FIELDSDEF=`dirname $0`/../annot/bedPeakFields.as
#chromInfo=`dirname $0`/../annot/chromInfo.txt
${BEDTOBIGBED} -type=bed6+6 -as=${FIELDSDEF} -tab ${OUTBED} ${chromInfo} ${OUTBIGBED}
echo -e "BigBed track saved to:\n${OUTBIGBED}"

