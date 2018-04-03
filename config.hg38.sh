# SPAR config file

export HOMEDIR=/project/wang/pkuksa 

#absolute path to the bin directory
export BINDIR="${HOMEDIR}/bin"

# (default) absolute path to the default SPAR output directory
export SPARDIR="${HOMEDIR}/SPAR_out"

# genome build
export GENOMEBUILD=hg38

# data directory for SPAR (DASHR/ENCODE data, etc)
export SPARDATAPATH=/project/wang/pkuksa/datasets/SPAR_data

#absolute path to the STAR genome index
export STAR="${BINDIR}/STAR_2.4.0h/bin/Linux_x86_64/STAR" # STAR
export STAR="/appl/bin/STAR" # STAR
export genomeDir="${HOMEDIR}/datasets/${GENOMEBUILD}/star/"  # STAR genome index
export GENOMEFA="${HOMEDIR}/datasets/${GENOMEBUILD}/${GENOMEBUILD}.fa"
export CONSERVATIONTRACK="${HOMEDIR}/datasets/${GENOMEBUILD}/${GENOMEBUILD}.100way.phastCons.bw"

#absolute path to pre-installed STAR, samtools, AWK, etc

export SAMTOOLS="${BINDIR}/samtools-1.2/samtools" # SAMTOOLS
export BEDTOOLS=${BINDIR}/bedtools2/bin/bedtools

export GAWK="${BINDIR}/gawk-4.1.3/gawk"

# UCSC tools
export BGTOBIGWIG="${BINDIR}/bedGraphToBigWig"
export BEDTOBIGBED="${BINDIR}/bedToBigBed"
export BIGWIGTOBEDGRAPH="${BINDIR}/bigWigToBedGraph"
export BIGWIGAVGOVERBED="${BINDIR}/userApps/bin/bigWigAverageOverBed"

# Adapter trimming
export CUTADAPT="${BINDIR}/cutadapt-1.8.1/bin/cutadapt"

#hg 19 chromosome information file
#chromInfo=${SPARPATH}/annot/chromInfo.txt

#mapping parameters for STAR
export maxMismatchCnt=0 # maximum number of genomic mismatches
export maxMapCnt=100 # maximum number of places a read can map to
export minMappedLength=14 # minimum *mapped* length
export maxReadLength=44 # maximum read length
export minReadLength=14 # minimum read length
export minCoverage=10 # minimum read coverage
export minFoldChange=2 # minimum fold change (read coverage(i)/read_coverage(i-1)) 
                       # for peak detection
export max5pClip=1 # maximum allowed 5' end clipping ('S' in the CIGAR string) 


export RSCRIPT="${BINDIR}/R-3.2.3/bin/Rscript"
export RSCRIPT="/appl/R-3.2.3/bin/Rscript"

function printT
{
  echo "`date +'%b %d %H:%M:%S'` ..... $1"
}
export -f printT
