

SPARPEAKS=$1 # peaks from user data
DASHREXPR=$2 # DASHR expression matrix [loci x tissues]

if [ $# -lt 2 ]; then
  echo "USAGE: $0 <table.peaks> <table.dashr.expr_matrix>"
  exit 1
fi

minOverlap=0.9
intersectExprFile=${SPARPEAKS}.encode.expr.intersect.xls


# get headers for DASHR reference table and SPAR output peak table

peakTableHeader=$(head -n 1 ${SPARPEAKS} | awk '{$0=substr($0,2); print $0;}')
nColumnsPeakTable=$(awk 'BEGIN{FS="\t"}{print NF; exit}' ${SPARPEAKS})

if [[ ${DASHREXPR} == *.gz ]]; then
  nColumnsDashrTable=$( zcat "${DASHREXPR}" | awk 'BEGIN{FS="\t"}{print NF; exit}' )
  dashrExprHeader=$( zcat "${DASHREXPR}" | awk '{$0=substr($0,2); print $0; exit}' )
else
  nColumnsDashrTable=$( awk 'BEGIN{FS="\t"}{print NF; exit}' "${DASHREXPR}" )
  dashrExprHeader=$( head -n 1 "${DASHREXPR}" | awk '{$0=substr($0,2); print $0;}' )
fi


overlapField=$(( nColumnsPeakTable + nColumnsDashrTable + 1 ))
peakOverlapPctField=$(( overlapField+1 ))
dashrOverlapPctField=$(( overlapField+2 ))


overlapFile=${SPARPEAKS}.encode.peaks_overlap_raw
nooverlapFile=${SPARPEAKS}.encode.peaks_no_overlap


echo -e "#${peakTableHeader}\t${dashrExprHeader}\toverlap\tpeakOverlapPct\tENCODEpeakOverlapPct" > ${nooverlapFile}

# overlap
# stranded, reciprocal, min overlapping fraction -f 0.9
${BEDTOOLS} intersect -a ${SPARPEAKS} -b ${DASHREXPR} \
  -s -f ${minOverlap} -r  -wao -sorted | \
    awk 'BEGIN{ OFS="\t";
                dashrPeakStartIdx='${nColumnsPeakTable}'+2;
                overlapFile="'${overlapFile}'";
                nooverlapFile="'${nooverlapFile}'";
                overlapField='${overlapField}'+0;

              }
         {  
             if ($dashrPeakStartIdx == -1) # no overlap
             {
                 # add two overlap fields
                 print $0, 0.0, 0.0 >> nooverlapFile
             }
             else
             {
                 overlapSize=$overlapField
                 dashrPeakSize = $(dashrPeakStartIdx+1)-$dashrPeakStartIdx; # assuming 0-based start
                 sparPeakSize = $3-$2;
                 sparPeakOverlapPct = overlapSize / sparPeakSize;
                 dashrPeakOverlapPct = overlapSize / dashrPeakSize;
                 print $0, sparPeakOverlapPct, dashrPeakOverlapPct > overlapFile
             }
         }'


echo -e "#${peakTableHeader}\t${dashrExprHeader}\toverlap\tpeakOverlapPct\tDASHRpeakOverlapPct" > ${intersectExprFile}.full


#sort -k1,1 -k2,2n -k3,3n -k6,6 -k${peakOverlapPctField},${peakOverlapPctField}gr -k${dashrOverlapPctField},${dashrOverlapPctField}gr ${overlapFile} | sort -k1,1 -k2,2n -k3,3n -k6,6 -s -u >> ${intersectExprFile}.full

# note: output only best overlap
sort -k1,1 -k2,2n -k3,3n -k6,6 -k${peakOverlapPctField},${peakOverlapPctField}gr -k${dashrOverlapPctField},${dashrOverlapPctField}gr ${overlapFile} | awk 'BEGIN{sep="\t"; OFS="\t"}{h=($1 sep $2 sep $3 sep $6); if (a[h]!=1) {print; a[h]=1};}' >> ${intersectExprFile}.full
#awk 'BEGIN{sep="\t"; OFS="\t"; peakOverlapIdx='${peakOverlapPctField}'+0; refOverlapIdx='${dashrOverlapPctField}'}{h=($1 sep $2 sep $3 sep $6); curr_peak_pct=$peakOverlapIdx; curr_ref_pct=$refOverlapIdx; if (h_prev && h!=h_prev) {print line_max; peak_max=0; ref_max=0; line_max=""}; if (h==h_prev) {if (peak_max<curr_peak_pct) peak_max=curr_peak_pct; line_max=$0;}else if (peak_max==curr_peak_pct) {if (ref_max<curr_ref_pct) {ref_max=curr_ref_pct; line_max=$0;}}; h_prev=h;}END{ print line_max; }' ${overlapFile} >> ${intersectExprFile}.full

# shorten the overlap table
# output only peak location/RPM and DASHR columns
#awk 'BEGIN{nColumnsPeaksTable='${nColumnsPeakTable}'+0;}{ printf "%s\t%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$19,$6; for (i=nColumnsPeaksTable+1; i<=NF; ++i) printf "\t%s", $i; printf "\n";  }' ${intersectExprFile}.full > ${intersectExprFile}
cp -p ${intersectExprFile}.full ${intersectExprFile}

# for new, non-DASHR peaks, output default (zero) expression vector
#nPeakFields=$(awk 'BEGIN{FS="\t"}{print NF; exit}' ${SPARPEAKS})
#nDashrFields=$(awk 'BEGIN{FS="\t"}{print NF; exit}' ${DASHRTABLE})


#dashrExprFieldIdx=$(($nPeakFields + 5))
#peakExprFieldIdx=5

