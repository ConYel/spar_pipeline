# output read sequences with expression information
INBED=$1 # collapsed BED (intervals with expression, coordinate-sorted
GENOME=$2 # reference genome sequence (FASTA)
MINEXPR=${3:-10} # minimum expression level
FILTBED=${INBED}.min${MINEXPR}

if [ $# -lt 2 ]; then
  echo "USAGE: $0 <collapsed.BED> <genome.fa> <min_expression>"
  exit 1
fi

awk 'BEGIN{minExpr='${MINEXPR}'+0;}
     {if ($5 > minExpr) print}' ${INBED} > ${FILTBED}

# NOTE: this assumes reads are short enough to be printed on *one* line!
# sequences are from reference genome sequence FWD strand
#paste ${FILTBED} <(awk '{chrStart=$2+1; print $1":"chrStart"-"$3}' ${FILTBED} | xargs ${SAMTOOLS} faidx ${GENOME} | awk '{if ($1!~/^>/) print}') 
####| sort -k5,5nr # sort by expression


# NOTE: this assumes reads are short enough to be printed on *one* line!
# sequences are from reference genome sequence,
# reverse complemented for negative strand
#paste ${FILTBED} <(awk '{chrStart=$2+1; print $1":"chrStart"-"$3}' ${FILTBED} | xargs ${SAMTOOLS} faidx ${GENOME} | awk '{if ($1!~/^>/) print}') | awk 'BEGIN{OFS="\t"; rc["A"]="T"; rc["C"]="G"; rc["G"]="C"; rc["T"]="A"; rc["N"]="N"; rc["a"]="t"; rc["c"]="g"; rc["g"]="c"; rc["t"]="a"; rc["n"]="n"; }{strand=$6; if (strand=="-") { n=split($7,seq,""); revcompseq=""; for (i=n; i>0; --i) {nuc=rc[seq[i]]; if (nuc=="") nuc="N"; revcompseq = (revcompseq nuc);} if (length(revcompseq)!=n) {print "ERROR", $7, revcompseq; exit} ;$7=revcompseq;}; print}'
paste ${FILTBED} <(awk '{chrStart=$2+1; print $1":"chrStart"-"$3}' ${FILTBED} | xargs ${SAMTOOLS} faidx ${GENOME} | awk 'BEGIN{s=""}{if ($1!~/^>/) s=(s $0); else {if (NR==1) next; print s; s="";}}' ) | awk 'BEGIN{OFS="\t"; rc["A"]="T"; rc["C"]="G"; rc["G"]="C"; rc["T"]="A"; rc["N"]="N"; rc["a"]="t"; rc["c"]="g"; rc["g"]="c"; rc["t"]="a"; rc["n"]="n"; }{strand=$6; if (strand=="-") { n=split($7,seq,""); revcompseq=""; for (i=n; i>0; --i) {nuc=rc[seq[i]]; if (nuc=="") nuc="N"; revcompseq = (revcompseq nuc);} if (length(revcompseq)!=n) {print "ERROR", $7, revcompseq; exit} ;$7=revcompseq;}; print}'



