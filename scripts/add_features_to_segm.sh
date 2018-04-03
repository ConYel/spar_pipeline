# add RPM and quantile columns to segmentation output

INSEGMPOS=$1
INSEGMNEG=$2

exprValIdx=5

# combine segmentations for positive and negative strands
combinedSegm=`dirname ${INSEGMPOS}`/combinedSegm.tmp
>${combinedSegm}
if [ -s ${INSEGMPOS} ]; then
  cat ${INSEGMPOS} >> ${combinedSegm}
fi
if [ -s ${INSEGMNEG} ]; then
  cat ${INSEGMNEG} >> ${combinedSegm}
fi

# sort by expression in ascending order
sort -k${exprValIdx},${exprValIdx}g ${combinedSegm} -o ${combinedSegm}

# totalExpr ~ library size; totalSize = number of segmented loci
aSegmStats=( $(cut -f ${exprValIdx} ${combinedSegm} | awk '{totalExpr+=$1}END{print totalExpr, NR}') )
totalExpr=${aSegmStats[0]}
totalSize=${aSegmStats[1]}
echo "TotalExpr=${totalExpr}"
echo "TotalExpr=${totalSize}"
#> ${combinedSegm}.with_rpm_quantile
awk 'BEGIN{OFS="\t"; FS="\t"; totalExpr='${totalExpr}'+0; totalSize='${totalSize}'+0; exprValIdx='${exprValIdx}'+0; quantileSize = int(totalSize / 10); }
    {
        rpm = $exprValIdx / totalExpr * 1000000;
        #quantile = int(NR / quantileSize);
        #if (quantile>9) quantile = 9;
        #print $0, rpm, quantile;
        percentile = NR / totalSize;
        #print $0, rpm, percentile;
        printf "%s\t%.4f\t%.2f\n", $0, rpm, percentile*100;
    }' ${combinedSegm} | sort -k1,1 -k2,2n -k3,3n -k6,6 | \
awk 'BEGIN{outPos="'${INSEGMPOS}'"; outNeg="'${INSEGMNEG}'";}
    {
      strand = $6;
      if (strand=="+")
         print > outPos;
      else
         print > outNeg;
     }' 


