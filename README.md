# SPAR: Small RNA-seq Portal for Analysis of sequencing expeRiments

This repository contains SPAR pipeline scripts.
Provided *SPAR.sh* script automates processing, analysis, visualization of small RNA-seq data in mapped (BAM), genome-wide read coverage (bigWig), and *trimmed* sequencing reads (FASTQ) formats ( bigWig/BAM/trimmmed FASTQ files can be prepared using SPAR preprocessing scripts  [SPAR](https://www.lisanwanglab.org/SPAR)).

## Preparing reference genome (FASTA):

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

chmod a+x twoBitToFa

./twoBitToFa hg19.2bit hg19.fa
```

## Preparing STAR index for the reference genome

```
mkdir -p hg19/star

STAR --runMode genomeGenerate --genomeDir hg19/star --genomeFastaFiles hg19.fa --runThreadN 4
```

## Setting up configuration file

Example configuration files are given in the provided *config.hg19.sh*/config.hg38.sh/config.mm10.sh files for human and mouse genomes.
Please set locations of the programs/tools in the config file to the locations in your system.
The configuration files are also used to specify SPAR analysis parameters.
```
#SPAR config file

export HOMEDIR="${HOME}"

#absolute path to the bin directory
export BINDIR="${HOMEDIR}/bin"

#reference genome
export GENOMEBUILD=hg19

#absolute path to the STAR genome index
export STAR="${BINDIR}/STAR_2.4.0j/bin/Linux_x86_64/STAR" # STAR
export genomeDir="${HOMEDIR}/datasets/${GENOMEBUILD}/star/"  # STAR genome index
export GENOMEFA="${HOMEDIR}/datasets/${GENOMEBUILD}/${GENOMEBUILD}.fa"

#absolute path to pre-installed bedtools, samtools, AWK, etc
export SAMTOOLS="${BINDIR}/samtools-1.2/samtools"
export BEDTOOLS="${BINDIR}/bedtools-2.26/bedtools2/bin/bedtools"

# UCSC tools
export BGTOBIGWIG="${BINDIR}/bedGraphToBigWig"
export BEDTOBIGBED="${BINDIR}/bedToBigBed"

# Adapter trimming
export CUTADAPT="${BINDIR}/cutadapt-1.8.1/bin/cutadapt"

# Analysis parameters for STAR
export maxMismatchCnt=0 # maximum number of genomic mismatches
export maxMapCnt=100 # maximum number of places a read can map to
export minMappedLength=14 # minimum *mapped* length
export maxReadLength=44 # maximum read length
export max5pClip=1 # maximum allowed 5' clip; all read with 5' clipping > max will be discarded
export keep5pClipped=0 # by default all reads clipped at 5' are excluded from analysis
```

### Download/install if necessary any missing programs/tools:

STAR
```
wget https://github.com/alexdobin/STAR/archive/2.5.4b.tar.gz
```

UCSC tools
```
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed
```

Conservation tracks


Human genome (GRCh37/hg19):
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons100way/hg19.100way.phastCons.bw
```

Human genome (GRCh38/hg38):
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.bw
```

Mouse genome (GRCm38/mm10):
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/phastCons60way/mm10.60way.phastCons.bw
```

Please set CONSERVATIONTRACK in the SPAR config file to the absolute path of the conservation bigWig track file.

Samtools
```
apt-get install samtools
```



## Recompiling bam2bedgraph binary (if necessary)

### Linux:
Install htslib
```
wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2
tar xjvf htslib-1.7.tar.bz2
cd htslib-1.7
make
make install
```

Compile bam2bedgraph:
```
g++ -O2 bam2bedgraph.cpp -o bam2bedgraph_interval -lhts
```

### MAC OS:

```
brew install htslib

g++ -O2 bam2bedgraph.cpp -o bam2bedgraph_interval -lhts
```

## Preparing bigWig and BAM files from trimmed FASTQ

To prepare bigWig files and BAM from the example trimmed and gzipped FASTQ file:
```
bash prepare_BAM_and_bigWigs_from_fastq.sh example-data.fastq.gz test_out config.hg19.sh
```

## Running SPAR pipeline

SPAR pipeline (SPAR.sh script) usage:
```
SPAR.sh <reads.fastq.gz|aln.bam|signal.bigWig|URL.bam> <SPAR_output_directory> <SPAR_config.sh>
```

To process small RNA-seq data in *BAM* format:
```
./SPAR.sh input.BAM bam_out config.hg19.sh
```

To process small RNA-seq data in *bigWig* format:
```
./SPAR.sh input.pos.bigWig bigwig_out config.hg19.sh
```

To process web accessible (URL) small RNA-seq data:
```
./SPAR.sh http://dashr2.lisanwanglab.org/DASHRv2/tracks/hg19/DASHR1_GEO_hg19/adipose1_GSE45159.pos.bigWig adipose_out config.hg19.sh
```
