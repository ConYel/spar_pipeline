#!/bin/bash

## use smrna_adapter_cut.sh from spar_prepare workflow instead
set -e

CUTADAPT=$( command -v cutadapt )
if [ ! -x "${CUTADAPT}" ]; then
  echo "$0 requires cutadapt, but it is not found. Please install cutadapt and/or make sure it is in the path."
  exit 1
fi

if [ $# -lt 1 ]; then
  echo "USAGE: `basename $0` reads.fastq.gz [adapters in cutadapt format: -a ADAPTERSEQ3p, -g ADAPTERSEQ5p]"
  echo -e "Example:\n`basename $0` reads.fastq.gz -a TGGAATTCTCGGGTGCCAAGG"
  exit 1
fi

#>illumina_smrna_1.0_3p
smrna_1_0=TCGTATGCCGTCTTCTGCTTG
#>illumina_smrna_1.5_3p(Illumina_Small_RNA_3p_Adapter_1)
smrna_1_5=ATCTCGTATGCCGTCTTCTGCTTG
#>illumina_truseq
smrna_truseq=TGGAATTCTCGGGTGCCAAGG
smrna_ribo=AGATCGGAAGAGCACACGTCT


INFASTQ=$1
ADAPTERoption=$@
ADAPTERoption=${ADAPTERoption##${INFASTQ}}

validExtension=$( echo "${INFASTQ}" | awk '{if ( tolower($1) ~ /.fastq.gz/ ) print 1; else print 0; exit}' )

if [ "${validExtension}" != "1" ]; then
  echo "Invalid file format. File in *.fastq.gz format is required."
  exit 1;
fi

DATASET=${INFASTQ%.*.*}

# trimming parameters
MIN_TRIMMED_READ_LEN=14
MAX_ADAPTER_ERROR=0.06
MIN_ADAPTER_MATCH=6


TRIMMEDREADS=${DATASET}_trimmed.fastq.gz
UNTRIMMEDREADS=${DATASET}_untrimmed.fastq.gz
TOOSHORTREADS=${DATASET}_tooshort.fastq.gz


${CUTADAPT} ${ADAPTERoption} -e ${MAX_ADAPTER_ERROR} -o ${TRIMMEDREADS} \
--untrimmed-output ${UNTRIMMEDREADS} --too-short-output ${TOOSHORTREADS} \
-O ${MIN_ADAPTER_MATCH} -m ${MIN_TRIMMED_READ_LEN} -n 2 \
--length-tag "length=" ${INFASTQ}

