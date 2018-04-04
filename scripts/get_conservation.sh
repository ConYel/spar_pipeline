set -e
# the out.bed will be original bed + mean/min/max/ columns

#http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw

INTABLE=$1 # input table in BED-like format

CONSERVATIONTRACK=$2 # genome-wide conservation bigWig track

# NOTE: this is modified and re-compiled version of bigWigAverageOverBed
# main change: use bigWig fetch for each line instead of per-chromosome method
# the per-chromosome method is *slow* as it essentially constructs chromosome-sized vector from bigWig data
# the default behavior was to use fetching if the number of BED intervals was greater than 3000 and chromosome method otherwise
#BIGWIGAVGOVERBED=bin/userApps/bin/bigWigAverageOverBed

if [ $# -lt 2 ]; then
  echo "USAGE: $0 <table.bed> <conservation.bigWig>"
  exit 1
fi

if [ ! -x "${BIGWIGAVGOVERBED}" ]; then
  echo "ERROR: bigWigAverageOverBed not found or is not executable!"
  exit 1
fi

if [ ! -s "${CONSERVATIONTRACK}" ]; then
  echo "ERROR: ${CONSERVATIONTRACK} conservation track not found!"
  exit 1
fi

if [ ! -s "${INTABLE}" ]; then
  echo "ERROR: ${INTABLE} not found!"
  exit 1
fi


INBED=${INTABLE}.bwAvgOverBed.prepared.bed
OUTTAB=${INBED}.output.tab
OUTBED=${INBED}.output.bed

# prepare BED file for bigWigAverageOverBed
cut -f1-6 ${INTABLE} | \
  awk 'BEGIN{OFS="\t"}
       { 
         if (NR==1) next; # skip header
         $5 = int($5); # make sure the score column is an integer
         strand="Plus";
         if ($6=="-") strand="Minus";
         $4 = ($4 ":" "LOC" NR-1 ":" strand); # make sure names are unique
         print;
       }' > ${INBED}

# run bigWigAverageOverBed
${BIGWIGAVGOVERBED} ${CONSERVATIONTRACK} ${INBED} ${OUTTAB} -bedOut=${OUTBED} -minMax

# add conservation columns to the original table
#paste ${INTABLE} <(echo -e "chr\tchrStart\tchrEnd\tname\tscore\tstrand\tconsMean0\tconsMin\tconsMax" | cat - ${OUTBED} | cut -f7-9)
paste ${INTABLE} <(echo -e "chr\tchrStart\tchrEnd\tname\tscore\tstrand\tconsMean0\tconsMin\tconsMax" | cat - ${OUTBED} | cut -f7-9 | awk 'BEGIN{OFS="\t"}{if (NR==1) { print; next }; for (i=1; i<NF; ++i) printf "%.4f\t", $i; printf "%.4f\n", $NF; }')

