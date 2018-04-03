ANNOTTABLE=$1
UNANNOTTABLE=$2
GENOMEBUILD=${3:-hg19}
#baseURL=${4}

HTMLHEADER="`dirname $0`/../templates/header_html_peak_table.html"
HTMLFOOTER="`dirname $0`/../templates/footer_html_peak_table.html"

cat ${HTMLHEADER}

#echo "<table id=\"peakTable\">"
echo "<table id=\"peakTable\" style=\"display:none\">"

#awk 'BEGIN{ncols=split("1,2,3,4,5,6,19,20,8,11,13,14,18",cols,",");annotStartIdx=21;}{sep="<td>"; if (NR==1) {sep="<th>"; print "<thead>"} printf "<tr>"; gsub(/#/,"",$1); for (i=1; i<=ncols; ++i) printf (sep $cols[i]); for (i=NF-2; i<=NF; ++i) printf (sep $i); for (i=annotStartIdx; i<=annotStartIdx+7; ++i) printf (sep $i); printf "\n"; if (NR==1) {print "</thead>"; print "<tbody>"};}' ${ANNOTTABLE}

#awk 'BEGIN{ncols=split("1,2,3,4,5,6,19,20,8,11,13,14,18",cols,",");annotStartIdx=21;nannotCols=split("NA,-1,-1,NA,OTHER,NA,-1,-1", defaultAnnot,",")}{if (NR==1) next; sep="<td>"; defaultAnnot[6]=$6; printf "<tr>"; for (i=1; i<=ncols; ++i) printf (sep $cols[i]); for (i=NF-2; i<=NF; ++i) printf (sep $i); for (i=1; i<=nannotCols; ++i) printf (sep defaultAnnot[i]); printf "\n";}' ${UNANNOTTABLE}

#ucscURLtemplate="https://genome.ucsc.edu/cgi-bin/hgTracks?db=${GENOMEBUILD}&position=%POSITION%&hgct_customText=track type=bigWig name=peakPos description="" bigDataUrl=https://tesla.pcbi.upenn.edu/~pkuksa/SPAR/SPAR_out/96da941/input.pos.bigWig"


#ucscURLtemplate="http://genome.ucsc.edu/cgi-bin/hgTracks?db=${GENOMEBUILD}&position=%POSITION%&hgt.customText=${baseURL}/ucsc_tracks.txt";

ucscURLtemplate="https://genome.ucsc.edu/cgi-bin/hgTracks?db=${GENOMEBUILD}&position="

#https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1:16861813-16861845&hgct_customText=track type=bigWig name=pos description="bibBedTrackPos" bigDataUrl=https://tesla.pcbi.upenn.edu/~pkuksa/SPAR/SPAR_out/96da941/input.pos.bigWig

#awk 'BEGIN{ucsc="'${ucscURLtemplate}'"; ncols=split("1,2,3,4,5,6,19,20,8,11,13,14,18",cols,",");annotStartIdx=21;}{sep="<td>"; if (NR==1) {sep="<th>"; print "<thead>"} printf "<tr>"; gsub(/#/,"",$1); coord=($1 ":" $2+1 "-" $3); if (NR>1) $4=("<a href=\"" ucsc coord "\" target=_blank>" $4 "</a>"); for (i=1; i<=ncols; ++i) printf (sep $cols[i]); for (i=NF-2; i<=NF; ++i) printf (sep $i); for (i=annotStartIdx; i<=annotStartIdx+7; ++i) printf (sep $i); printf "\n"; if (NR==1) {print "</thead>"; print "<tbody>"};}' ${ANNOTTABLE}
# 06oct2016: changed formatting for annotation ID
#awk 'BEGIN{ucsc="'${ucscURLtemplate}'"; ncols=split("1,2,3,4,5,6,19,20,8,11,13,14,18",cols,",");annotStartIdx=21;}{sep="<td>"; if (NR==1) {sep="<th>"; print "<thead>"} printf "<tr>"; gsub(/#/,"",$1); coord=($1 ":" $2+1 "-" $3); if (NR>1) $4=("<a href=\"" ucsc coord "\" target=_blank>" $4 "</a>"); for (i=1; i<=ncols; ++i) printf (sep $cols[i]); for (i=NF-2; i<=NF; ++i) printf (sep $i); if (NR>1) {annotIDidx=annotStartIdx+3; split($annotIDidx,a,":"); $annotIDidx=(a[1] ":" a[2] "-" a[3] ":" a[4] ":" a[5])}; for (i=annotStartIdx; i<=annotStartIdx+7; ++i) printf (sep $i); printf "\n"; if (NR==1) {print "</thead>"; print "<tbody>"};}' ${ANNOTTABLE}

#manual re-arrangement of columns
#awk 'BEGIN{FS="\t"; ucsc="'${ucscURLtemplate}'"; ncols=split("1,2,3,4,5,6,19,20,8,11,13,14,18",cols,",");annotStartIdx=21;}{sep="<td>"; if (NR==1) {sep="<th>"; print "<thead>"} printf "<tr>"; gsub(/#/,"",$1); coord=($1 ":" $2+1 "-" $3); if (NR>1) $4=("<a href=\"" ucsc coord "\" target=_blank>" $4 "</a>"); for (i=1; i<=ncols; ++i) printf (sep $cols[i]); for (i=annotStartIdx+8; i<=annotStartIdx+10; ++i) printf (sep $i); if (NR>1) {annotIDidx=annotStartIdx+3; split($annotIDidx,a,":"); $annotIDidx=(a[1] ":" a[2] "-" a[3] ":" a[4] ":" a[5])}; for (i=annotStartIdx; i<=annotStartIdx+7; ++i) printf (sep $i); for (i=annotStartIdx+11; i<=NF; ++i) print (sep $i); printf "\n"; if (NR==1) {print "</thead>"; print "<tbody>"};}' ${ANNOTTABLE}

# no re-arrangements
#awk 'BEGIN{FS="\t"; ucsc="'${ucscURLtemplate}'"; ncols=split("1,2,3,4,5,6,19,20,8,11,13,14,18",cols,",");annotStartIdx=21; colNames["repeat_family"]="Repeat?"; colNames["peakProportionOfReadsAtMostCommon3pPosition"]="Same 3p end read %"; colNames["peakProportionOfReadsAtMostCommon5pPosition"]="Same 5p end read %";}{sep="<td>"; if (NR==1) {sep="<th>"; print "<thead>"} printf "<tr>"; gsub(/#/,"",$1); gsub(/^X\./,"",$1); coord=($1 ":" $2+1 "-" $3); if (NR>1) $4=("<a href=\"" ucsc coord "\" target=_blank>" $4 "</a>"); if (NR>1) {annotIDidx=annotStartIdx+3; split($annotIDidx,a,":"); $annotIDidx=(a[1] ":" a[2] "-" a[3] ":" a[4] ":" a[5])}; for (i=1; i<=NF; ++i) { if (NR==1) { sep=("<th class=\"tooltip\" alt=\"" $i "\">"); hdr=$i; if (length(colNames[hdr])>1) $i=colNames[hdr];}; printf (sep $i);}  printf "\n"; if (NR==1) {print "</thead>"; print "<tbody>"};}' ${ANNOTTABLE}
colNames="`dirname $0`/../../SPAR_interactive_table_column_names.txt"
awk 'BEGIN{FS="\t"; ucsc="'${ucscURLtemplate}'"; ncols=split("1,2,3,4,5,6,19,20,8,11,13,14,18",cols,",");annotStartIdx=21; colNamesFile="'${colNames}'"; while(( getline line<colNamesFile) > 0 ) { split(line,a,"\t"); colNames[a[1]]=a[2]; }}{sep="<td>"; if (NR==1) {sep="<th>"; print "<thead>"} printf "<tr>"; gsub(/#/,"",$1); gsub(/^X\./,"",$1); coord=($1 ":" $2+1 "-" $3); if (NR>1) $4=("<a href=\"" ucsc coord "\" target=_blank>" $4 "</a>"); if (NR>1) {annotIDidx=annotStartIdx+3; split($annotIDidx,a,":"); $annotIDidx=(a[1] ":" a[2] "-" a[3] ":" a[4] ":" a[5])}; for (i=1; i<=NF; ++i) { if (NR==1) { sep=("<th class=\"tooltip\" alt=\"" $i "\">"); hdr=$i; if (length(colNames[hdr])>1) $i=colNames[hdr];}; printf (sep $i);}  printf "\n"; if (NR==1) {print "</thead>"; print "<tbody>"};}' ${ANNOTTABLE}

#awk 'BEGIN{ucsc="'${ucscURLtemplate}'"; ncols=split("1,2,3,4,5,6,19,20,8,11,13,14,18",cols,",");annotStartIdx=21;nannotCols=split("NA,-1,-1,NA,OTHER,NA,-1,-1", defaultAnnot,",")}{if (NR==1) next; coord=($1 ":" $2+1 "-" $3); $4=("<a href=\"" ucsc coord "\" target=_blank>" $4 "</a>"); sep="<td>"; defaultAnnot[6]=$6; printf "<tr>"; for (i=1; i<=ncols; ++i) printf (sep $cols[i]); for (i=NF-2; i<=NF; ++i) printf (sep $i); for (i=1; i<=nannotCols; ++i) printf (sep defaultAnnot[i]); printf "\n";}' ${UNANNOTTABLE}

# manual re-arrangement
#awk 'BEGIN{FS="\t"; ucsc="'${ucscURLtemplate}'"; ncols=split("1,2,3,4,5,6,19,20,8,11,13,14,18",cols,",");annotStartIdx=21;nannotCols=split("NA,-1,-1,NA,OTHER,NA,-1,-1", defaultAnnot,",")}{if (NR==1) next; coord=($1 ":" $2+1 "-" $3); $4=("<a href=\"" ucsc coord "\" target=_blank>" $4 "</a>"); sep="<td>"; defaultAnnot[6]=$6; printf "<tr>"; for (i=1; i<=ncols; ++i) printf (sep $cols[i]); for (i=annotStartIdx+8; i<=annotStartIdx+10; ++i) printf (sep $i); for (i=1; i<=nannotCols; ++i) printf (sep defaultAnnot[i]); for (i=annotStartIdx+11; i<=NF; ++i) print (sep $i); printf "\n";}' ${UNANNOTTABLE}

awk 'BEGIN{FS="\t"; ucsc="'${ucscURLtemplate}'"; ncols=split("1,2,3,4,5,6,19,20,8,11,13,14,18",cols,",");annotStartIdx=21;nannotCols=split("NA,-1,-1,NA,OTHER,NA,-1,-1", defaultAnnot,",")}{if (NR==1) next; coord=($1 ":" $2+1 "-" $3); $4=("<a href=\"" ucsc coord "\" target=_blank>" $4 "</a>"); sep="<td>"; defaultAnnot[6]=$6; printf "<tr>"; for (i=1; i<=NF; ++i) printf (sep $i);  printf "\n";}' ${UNANNOTTABLE}




echo "</tbody>"
echo "</table>"

cat ${HTMLFOOTER}

# https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr18:56118320-56118341
