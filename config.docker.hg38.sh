# SPAR config file

export HOMEDIR="/root/miniconda" 
export OUTRESDIR="my_data"

#absolute path to the bin directory
export BINDIR="${HOMEDIR}/bin"

# (default) absolute path to the default SPAR output directory
export SPARDIR="/home/${OUTRESDIR}"

# genome build
export GENOMEBUILD=hg38

# data directory for SPAR (DASHR/ENCODE data, etc)
export SPARDATAPATH="/home/spar_pipeline/annot"

#absolute path to the STAR genome index
export STAR="${BINDIR}/STAR" # STAR
#export STAR="/usr/bin/STAR" # STAR
export genomeDir="/home/${OUTRESDIR}/genome/hg38/star"  # STAR genome index
export GENOMEFA="/home/${OUTRESDIR}/genome/hg38.fa"
export CONSERVATIONTRACK="/home/${OUTRESDIR}/genome/hg38.phastCons100way.bw"

#absolute path to pre-installed STAR, samtools, AWK, etc

export SAMTOOLS="${BINDIR}/samtools" # SAMTOOLS
export BEDTOOLS=${BINDIR}/bedtools

export GAWK="/usr/bin/gawk"

# UCSC tools
export BGTOBIGWIG="${BINDIR}/bedGraphToBigWig"
export BEDTOBIGBED="${BINDIR}/bedToBigBed"
export BIGWIGTOBEDGRAPH="${BINDIR}/bigWigToBedGraph"
export BIGWIGAVGOVERBED="${BINDIR}/bigWigAverageOverBed"

# Adapter trimming
export CUTADAPT="${BINDIR}/cutadapt"

#mapping parameters for STAR
export maxMismatchCnt=1 # maximum number of genomic ~mismatches
export maxMapCnt=100 # maximum number of places a read can map to
export minMappedLength=14 # minimum *mapped* length
export maxReadLength=44 # maximum read length
export minReadLength=14 # minimum read length
export minCoverage=10 # minimum read coverage
export minFoldChange=2 # minimum fold change (read coverage(i)/read_coverage(i-1)) 
                       # for peak detection
export max5pClip=1 # maximum allowed 5' end clipping ('S' in the CIGAR string) 


export RSCRIPT="/usr/bin/Rscript"
#export RSCRIPT="${BINDIR}/Rscript"

function printT
{
  echo "`date +'%b %d %H:%M:%S'` ..... $1"
}
export -f printT
