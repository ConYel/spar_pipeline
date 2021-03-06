# SPAR: Small RNA-seq Portal for Analysis of sequencing expeRiments

(Forked from wanglab-upenn bitbucket, https://bitbucket.org/wanglab-upenn/spar_pipeline/src/master/)

Analysis, annotation, visualization of small RNA sequencing experiments

Fast turn-around analyses for small RNA-seq, microRNA-seq, single-cell small RNA, and short total RNA sequencing data, as well as comparison of user-given sequencing data to reference data from DASHR normal human tissues/cell types and human ENCODE cell lines

* Ab initio (annotation-free) small non-coding RNA (sncRNA) discovery and characterization
* Expression comparison of user data with >180 reference human tissues, cell types and cell lines from DASHR and ENCODE
* Annotations for GRCh37/hg19, GRCh38/hg38, mm10

### Supported inputs and protocols:
* Small RNA-sequencing, microRNA-seq, short total RNA sequencing protocols
* Mapped reads (BAM), raw signal tracks (BigWig), and sequencing reads (FASTQ)
* Genome-wide or user-provided target region analysis


This repository contains SPAR pipeline scripts.
Provided *SPAR.sh* script automates processing, analysis, visualization of small RNA-seq data in mapped (BAM), genome-wide read coverage (bigWig), and *trimmed* sequencing reads (FASTQ) formats ( bigWig/BAM/trimmmed FASTQ files can be prepared using SPAR preprocessing scripts  [SPAR](https://www.lisanwanglab.org/SPAR)).



## Citation
"SPAR: Small RNA-seq Portal for Analysis of sequencing expeRiments". Pavel P. Kuksa, Alexandre Amlie-Wolf, Zivadin Katanic, Otto Valladares, Li-San Wang, Yuk Yee Leung.
https://academic.oup.com/nar/article/46/W1/W36/4992647#supplementary-data

## License
SPAR is available for academic and nonprofit use for free ([MIT license](LICENSE.md)).

## Prerequisites

1. Linux/Mac
2. Samtools >=1.2
3. Bedtools >=2.26
4. UCSC tools (bedGraphToBigWig, bigWigToBedGraph, bedToBigBed, bigWigAverageOverBed)
4. R >=3.2.3 and R libraries (parallel, ggplot2, dplyr, plyr, reshape2, RColorBrewer)
5. STAR >=2.4
6. Cutadapt >=1.9

## Installation

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

## Preparing conservation tracks

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

Samtools
```
apt-get install samtools
```

Cutadapt
```
pip install --user --upgrade cutadapt
```
or install from source
```
https://github.com/marcelm/cutadapt
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

## Preparing trimmed FASTQ

To trim adapter sequences from sequencing reads, *trim_adapter.sh* script can be used:
For example, to trim Illumina 3' adapter sequence:
```
bash scripts/trim_adapter.sh reads.fastq.gz -a TGGAATTCTCGGGTGCCAAGG
```
The trimming options are in CUTADAPT format (e.g., -a for 3', -b for 5' adapters, etc).
This will generate trimmed FASTQ (reads_trimmed.fastq.gz).

## Preparing bigWig and BAM files from trimmed FASTQ
SPAR can be run using trimmed FASTQ, bigWig, or BAM as input.

To prepare bigWig files and BAM from the trimmed and gzipped FASTQ file  :
```
bash prepare_BAM_and_bigWigs_from_fastq.sh example-data.fastq.gz test_out config.hg19.sh
```
(see [SPAR_prepare](https://bitbucket.org/wanglab-upenn/spar_prepare) repository and [README](https://bitbucket.org/wanglab-upenn/spar_prepare/src/master/README.md))
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
./SPAR.sh http://tesla.pcbi.upenn.edu/~pkuksa/SPAR/sample/hg19/sample.bam config.hg19.sh
```

### Running in a folder
small RNA adapter trimming with 10 processors
```
for file in my_data/*.fastq.gz;do echo ./spar_prepare/smrna_adapter_cut.sh $file 10;done 
```
optional, copy database files
```
cp hg38.fulltable.no_mRNA_no_lncRNA_piRNAonly.unique_LOC_pirDB.bed spar_pipeline/annot/hg38/hg38.fulltable.no_mRNA_no_lncRNA_piRNAonly.unique_LOC.bed 

cp hg38.fulltable.no_mRNA_no_lncRNA.unique_LOC_pirDB.bed spar_pipeline/annot/hg38/hg38.fulltable.no_mRNA_no_lncRNA.unique_LOC.bed 

```

SPAR pipeline (SPAR.sh script) with 10 processors
```
for file in my_data/*trimmed.fastq.gz;do echo ./spar_pipeline/SPAR.sh $file ./my_data/results_file ./spar_pipeline/config.docker.hg38.sh  10;done 

```

