#!/usr/bin/env bash
set -exv
set -o pipefail

export LC_ALL=C
export LC_numeric="en_US.UTF-8"

if [ $# -lt 3 ]; then
  echo "USAGE: `basename $0` <reads.fastq.gz|aln.bam|signal.bigWig|URL.bam> <SPAR_output_directory> <SPAR_config.sh>"
  exit 1
fi

sparCommandLine="$0 $@"

configFile=${3:-config.sh}

if [ -s ${configFile} ]; then
  source "${configFile}"
else
  echo "ERROR: ${configFile} not found."
  exit 1
fi

genomeBuild=${GENOMEBUILD:-"hg19"}
conservationTrack=${CONSERVATIONTRACK}
doConservation=1 #${doConservation:-1}
doCleanup=1

STARTTIME=$(date +%s)

# input data file
INFILE=$1

# set SPARDIR output directory to specified directory if any, or use
# default directory from config file
SPARDIR=${2:-${SPARDIR}} 

# path do SPAR root
SPARPATH=`dirname $0`

# path to SPAR data, e.g., DASHR/ENCODE data tables
SPARDATAPATH=${SPARDATAPATH:-SPAR_data}


segmHeaderFile=`dirname $0`/templates/header_unannot_table.txt
annotHeaderFile=`dirname $0`/templates/header_annot_table.txt


urlRegex='(https?|ftp|file)://[-A-Za-z0-9\+&@#/%?=~_|!:,.;]*[-A-Za-z0-9\+&@#/%=~_|]'
isURL=0
if [[ ${INFILE} =~ ${urlRegex} ]]; then
  isURL=1
fi



if [ ${isURL} = 1 ]; then

 # check for internet connection
 if ping -c 1 google.com > /dev/null 2>&1; then
   echo "";
 else
   echo "****ERROR: no internet connection";
   exit 1;
 fi

 status=$(curl -s --head -w %{http_code} ${INFILE} -o /dev/null)
 if [ "${status}" = "000" ]; then
   echo -e "****ERROR: network error. status=${status}"
   exit 1
 fi
 if [ ! ${status} -lt 400 ]; then
   echo -e "****ERROR: URL ${INFILE} not found: http_code=${status}"
   exit 1
 fi
elif [ ! -s "${INFILE}" ]; then
  echo -e "*****ERROR: INPUT file\n${INFILE}\nnot found or empty!"
  exit 1
fi

# determine file type
# if BAM, skip mapping, proceed
# if FASTQ, map, then proceed from BAM to BEDGRAPH conversion
# if BIGWIG, proceed from segmentation

FILETYPE=${INFILE##*.}

# detect compressed input
isCompressed=0
if [ "${FILETYPE}" = "gz" ] || [ "${FILETYPE}" = "bz2" ]; then
  isCompressed=1
  FILETYPE=$(echo "${INFILE}" | awk '{ext2=$0; gsub(/^.+\./,"",ext2); ext1=$0; gsub(/\.[^.]+$/,"",ext1); gsub(/^.+\./,"", ext1); print ext1;}')
  ARCHIVETYPE=$(echo "${INFILE}" | awk '{ext2=$0; gsub(/^.+\./,"",ext2); ext1=$0; gsub(/\.[^.]+$/,"",ext1); gsub(/^.+\./,"", ext1); print ext2;}')
fi


lcFILETYPE=$(echo ${FILETYPE} | awk '{print tolower($0)}') # lower-case 
ucFILETYPE=$(echo ${FILETYPE} | awk '{print toupper($0)}') # upper-case 

isBAM=0
isBIGWIG=0
isFASTQ=0

if [ "${lcFILETYPE}" = "bam" ]; then
  isBAM=1
fi

if [ "${lcFILETYPE}" = "bigwig" ]; then
  isBIGWIG=1
fi

if [ "${lcFILETYPE}" = "fastq" ]; then
  isFASTQ=1
fi

if [ ${isURL} = 1 ] && [ ${isBAM} = 1 ]; then
   isURLBAM=1
fi



EXPNAME=`basename "${INFILE}" ".${FILETYPE}"`
if [ "${isCompressed}" = 1 ]; then
  EXPNAME=`basename "${INFILE}" ".${FILETYPE}.${ARCHIVETYPE}"`
fi


if [ "${isFASTQ}" = 1 ]; then
  OUTDIR=${SPARDIR}/${EXPNAME}_m${maxMismatchCnt}_map${maxMapCnt}
  OUTBAM=${OUTDIR}/mapping/Aligned.out.filtered.hardClipped.sorted.bam
  mkdir -p ${OUTDIR}/mapping
elif [ "${isBAM}" = 1 ]; then 
  # if input is mapped BAM file
  if [ -z "$2" ]; then
    # if no output directory is provided, use BAM file directory
    OUTDIR=`dirname ${INFILE}`  #${SPARDIR}/${EXPNAME}
  else
    # if output directory is specified, use it
    OUTDIR=$2
  fi
  OUTBAM=${INFILE}
elif [ "${isBIGWIG}" = 1 ]; then
  if [ -z "$2" ]; then
    # if no output directory is provided, use bigWig file directory
    OUTDIR=`dirname ${INFILE}`  #${SPARDIR}/${EXPNAME}
  else
    # if output directory is specified, use it
    OUTDIR=$2
  fi
  OUTBAM=${INFILE}
fi

# OUTPUTS and OUTPUT DIR STRUCTURE
if [ -z "${OUTDIR}" ]; then
  OUTDIR=SPAR_output
fi
if [ "${OUTDIR}" = "${SPARPATH}" ]; then
  OUTDIR=SPAR_output
fi
outputDirs="tracks figures results logs inputs DASHR_comparison ENCODE_comparison results/byClass"
for d in ${outputDirs}; do
  mkdir -p ${OUTDIR}/${d}
done

mkdir -p ${OUTDIR}/figures/js
cp -p ${SPARPATH}/templates/Chart.PieceLabel.js ${OUTDIR}/figures/js/
mkdir -p ${OUTDIR}/figures/css
cp -p ${SPARPATH}/templates/charts-download.css ${OUTDIR}/figures/css/

# JBrowse peak tracks
# these are bigWigs! since JBrowse did not support bigBed!!
JBPEAKTRACKPLUS=${OUTDIR}/tracks/peaks.pos.jbrowse.bigWig
JBPEAKTRACKMINUS=${OUTDIR}/tracks/peaks.neg.jbrowse.bigWig

MATSEQANNOT=${OUTDIR}/results/mature_seqs_annot.xls
MATSEQUNANNOT=${OUTDIR}/results/mature_seqs_unannot.xls
MATSEQALL=${OUTDIR}/results/mature_seqs_all.xls

PEAKTRACKPLUS=${OUTDIR}/tracks/peaks.pos.bigBed
PEAKTRACKMINUS=${OUTDIR}/tracks/peaks.neg.bigBed
PEAKBEDALL=${OUTDIR}/results/peaks_all.bed
RAWPLUS=${OUTDIR}/tracks/raw.pos.bigWig
RAWPLUSBG=${OUTDIR}/tracks/raw.pos.bedGraph
RAWMINUS=${OUTDIR}/tracks/raw.neg.bigWig
RAWMINUSBG=${OUTDIR}/tracks/raw.neg.bedGraph
PEAKTABLEANNOT=${OUTDIR}/results/peaks_annot.xls
PEAKTABLEUNANNOT=${OUTDIR}/results/peaks_unannot.xls
PEAKTABLEALL=${OUTDIR}/results/peaks_all.xls
GENEEXPRTABLE=${OUTDIR}/results/smRNA_gene_expression.xls
DASHRCOMPANNOT=${OUTDIR}/DASHR_comparison/DASHR_comparison_annot.xls
DASHRCOMPUNANNOT=${OUTDIR}/DASHR_comparison/DASHR_comparison_unannot.xls
DASHRCOMPALL=${OUTDIR}/DASHR_comparison/DASHR_comparison_all.xls
NOTINDASHR=${OUTDIR}/DASHR_comparison/not_in_DASHR.xls
ENCODECOMPANNOT=${OUTDIR}/ENCODE_comparison/ENCODE_comparison_annot.xls
ENCODECOMPUNANNOT=${OUTDIR}/ENCODE_comparison/ENCODE_comparison_unannot.xls
ENCODECOMPALL=${OUTDIR}/ENCODE_comparison/ENCODE_comparison_all.xls
NOTINENCODE=${OUTDIR}/ENCODE_comparison/not_in_ENCODE.xls
HTMLOUTPUT=${OUTDIR}/results.html

resultsTemplate="`dirname $0`/templates/spar_results_template.html"
if [ "${isBAM}" = 1 ]  || [ "${isFASTQ}" = 1 ]; then
  resultsTemplate="`dirname $0`/templates/spar_results_template.with_mat_seqs.html"
fi 
cat ${resultsTemplate} | \
while read -r line; do
    while [[ "$line" =~ ([{][A-Z]*[}]) ]] ; do
        LHS=${BASH_REMATCH[1]}
        RHS="$(eval echo "\"\$${LHS}\"")"
        line=${line//\$${LHS}/${RHS}}
    done
    echo "$line"
done > ${HTMLOUTPUT} 

awk 'BEGIN{d="'${OUTDIR}'/";}{gsub(d,"",$0); print $0;}' ${HTMLOUTPUT} > ${HTMLOUTPUT/.html/.local.html}

keep5pClipped=""
# parse command line arguments
shift 3
while [ "$#" -gt 0 ]; do
  case "$1" in
    --min_read_length=*)
        minReadLength="${1#*=}"
        ;;
    --max_5p_clip=*)
        max5pClip="${1#*=}"
        ;;
    --min_coverage=*)
        minCoverage="${1#*=}"
        ;;
    --max_multi_map=*)
        maxMapCnt="${1#*=}"
        ;;
    --keep5pClipped=*)
        keep5pClipped="${1#*=}"
        ;;
    --help) print_help;;
    *)
      printf "***********************************************************\n"
      printf "* Error: Invalid argument, run --help for valid arguments. *\n"
      printf "***********************************************************\n"
      exit 1
  esac
  shift
done

if [ "${keep5pClipped}" = 1 ]; then
  keep5pClipped="--keep5pClipped"
fi

echo "FILETYPE=${FILETYPE}; isCompressed=${isCompressed}; isBAM=${isBAM}; isURL=${isURL} EXPNAME=${EXPNAME}; OUTDIR=${OUTDIR}; OUTBAM=${OUTBAM}; keep5pClipped=${keep5pClipped}"
mkdir -p ${OUTDIR}


if [ "${isBAM}" = 1 ]; then
  if [ "${isURL}" = 1 ]; then
    OUTBAM=${OUTDIR}/input.bam
  fi
fi



LOGSPAR=${OUTDIR}/`basename ${OUTBAM}`.SPAR.log
RLOG=${OUTDIR}/`basename ${OUTBAM}`.R.log
LOGSPAR=${OUTDIR}/logs/SPAR.log
RLOG=${OUTDIR}/logs/R.log
>${LOGSPAR}
>${RLOG}


function runScript
{
  bash ${SPARPATH}/scripts/$1
}

function printT
{
  echo "`date +'%a, %d %b %Y %H:%M:%S %z'` ..... $1" | tee -a ${LOGSPAR}
}

function printL
{
  echo -e "$1" | tee -a ${LOGSPAR}
}



printL "GENOME BUILD: ${genomeBuild}\n"
printL "Output directory:\n${OUTDIR}\n"
printL "STAR GENOME directory: ${genomeDir}\n"
printL "minRead=${minReadLength}, maxRead=${maxReadLength}, max5pClip=${max5pClip}, minCoverage=${minCoverage}"


printT "SPAR ${ucFILETYPE} run started"

# map reads if FASTQ is given at the input (this will create BAM file)
if [ "${isFASTQ}" = 1 ]; then
  runScript "run_star_smrna2.sh ${INFILE} ${maxMismatchCnt} ${maxMapCnt} ${OUTDIR}/mapping"
fi

# convert BAM to BEDGRAPH
if [ "${isBAM}" = 1 ] || [ "${isFASTQ}" = 1 ]; then
  printT "Converting BAM to bedGraph"
  if [ "${isURL}" = 1 ]; then
       printT "bash ${SPARPATH}/scripts/bamURL_to_bedgraph2_interval_cut_cpp.sh ${INFILE} ${OUTDIR}/input.bam ${minReadLength} ${maxReadLength} ${keep5pClipped}"
       bash ${SPARPATH}/scripts/bamURL_to_bedgraph2_interval_cut_cpp.sh ${INFILE} ${OUTDIR}/input.bam ${minReadLength} ${maxReadLength} ${keep5pClipped}
  else
    printT "bam_to_bedgraph2_interval_cut_cpp.sh ${OUTBAM} ${OUTDIR}/mapping ${minReadLength} ${maxReadLength} ${keep5pClipped}"
    runScript "bam_to_bedgraph2_interval_cut_cpp.sh ${OUTBAM} ${OUTDIR}/mapping ${minReadLength} ${maxReadLength} ${keep5pClipped}"
  fi
 
  if [ "${isBAM}" = 1 ]; then
    if [ "${isURL}" = 1 ]; then
      OUTBAM=${OUTDIR}/input.bam
    else
      OUTBAM=${OUTDIR}/mapping/`basename "${INFILE}"`
    fi
  fi
  chromInfoFile=${OUTBAM}.chromInfo.txt
  echo "${chrInfoFile}"


  printT "Converting bedGraph to bigWig"

  runScript "bedgraph_to_bigwig.sh ${OUTBAM}.pos.bedgraph ${chromInfoFile} ${RAWPLUS}"
  runScript "bedgraph_to_bigwig.sh ${OUTBAM}.neg.bedgraph ${chromInfoFile} ${RAWMINUS}"

fi

  chromInfoFile=`dirname $0`/annot/${GENOMEBUILD}/chromInfo.${GENOMEBUILD}.txt
if [ "${isBIGWIG}" = 1 ]; then
  OUTBAM=${OUTDIR}/`basename "${INFILE}"`
  OUTBAM=${OUTBAM%.*} # remove extension
  OUTBAM=${OUTBAM%.*} # remove strand, if any
  # convert bigwig to bedGraph
  INBIGWIG=${INFILE}
  # need to SET the strand for BIGWIG!!!!
  strand=${INBIGWIG%.*} # remove extension
  strand=${strand##*.}
  if [ "${strand}" != "pos" ] && [ "${strand}" != "neg" ]; then
     # set strand="positive" if strand information is not provided
     strand="pos"
  fi

  # bigWig does not have chromosome information
  chromInfoFile=`dirname $0`/annot/${GENOMEBUILD}/chromInfo.${GENOMEBUILD}.txt

  oppStrand="pos"
  OPPBIGWIG=${INBIGWIG}
  if [ "${strand}" = "pos" ]; then
    OPPBIGWIG=${OPPBIGWIG/.pos./.neg.}
    oppStrand="neg"
  else
    OPPBIGWIG=${OPPBIGWIG/.neg./.pos.}
    oppStrand="pos"
  fi
  
  doOppStrand=0
  if [ -s "${OPPBIGWIG}" ] && [ "${OPPBIGWIG}" != "${INBIGWIG}" ]; then
    doOppStrand=1
  fi


  OUTBEDGRAPH=${OUTBAM}.${strand}.bedgraph
  OUTBEDGRAPHOPP=${OUTBAM}.${oppStrand}.bedgraph
  echo "doOppStrand=${doOppStrand}; OPPBIGWIG=${OPPBIGWIG}"
  if [ "${doOppStrand}" = 0 ]; then
    ${BIGWIGTOBEDGRAPH} ${INBIGWIG} ${OUTBEDGRAPH} # convert bigWig to bedGraph 
    printT "Segmenting BIGWIG"
    runScript "segment_bedgraph_entropy_position2.sh ${OUTBEDGRAPH} ${strand}"

    # add features to segmentation output
    OUTBEDGRAPHPOS=${OUTBAM}.pos.bedgraph.segm
    OUTBEDGRAPHNEG=${OUTBAM}.neg.bedgraph.segm
    bash ${SPARPATH}/scripts/add_features_to_segm.sh ${OUTBEDGRAPHPOS} ${OUTBEDGRAPHNEG}

    numFields=$(awk '{print NF; exit}' ${OUTBEDGRAPH}.segm)
    printT "Creating called peaks track BIGWIG"
    #runScript "create_peak_track.sh ${OUTBEDGRAPH}.segm"
    runScript "create_peak_track.sh ${OUTBEDGRAPH}.segm ${chromInfoFile} ${OUTBEDGRAPH}.segm.bigBed"

    printT "Annotating peaks BIGWIG"
    runScript "annotate_segm31.sh ${OUTBEDGRAPH}.segm"
  else
    echo "${strand}"
    echo "${OUTBEDGRAPH}"
    echo "${OUTBEDGRAPHOPP}"
    # convert bigWig to bedGraph
    ${BIGWIGTOBEDGRAPH} ${INBIGWIG} ${OUTBEDGRAPH} 
    ${BIGWIGTOBEDGRAPH} ${OPPBIGWIG} ${OUTBEDGRAPHOPP} 
   
    if [ "${strand}" = "pos" ]; then
      cp  ${OUTBEDGRAPH} ${RAWPLUSBG}
      cp  ${OUTBEDGRAPHOPP} ${RAWMINUSBG}
    else
      cp  ${OUTBEDGRAPH} ${RAWMINUSBG}
      cp  ${OUTBEDGRAPHOPP} ${RAWPLUSBG}
    fi
    printT "Segmenting BIGWIGs"
    runScript "segment_bedgraph_entropy_position2.sh ${OUTBEDGRAPH} ${strand}"
    runScript "segment_bedgraph_entropy_position2.sh ${OUTBEDGRAPHOPP} ${oppStrand}"
    # add features to segmentation output
    
    OUTBEDGRAPHPOS=${OUTBAM}.pos.bedgraph.segm
    OUTBEDGRAPHNEG=${OUTBAM}.neg.bedgraph.segm
    wc -l ${OUTBAM}.*.bedgraph
    wc -l ${OUTBAM}.*.bedgraph.segm
    bash ${SPARPATH}/scripts/add_features_to_segm.sh ${OUTBEDGRAPHPOS} ${OUTBEDGRAPHNEG}
 

    numFields=$(awk '{print NF; exit}' ${OUTBEDGRAPH}.segm)
    printT "Creating called peaks track BIGWIG"
    echo "${chromInfoFile}"
    runScript "create_peak_track.sh ${OUTBEDGRAPHPOS} ${chromInfoFile} ${PEAKTRACKPLUS}"
    runScript "create_peak_track.sh ${OUTBEDGRAPHNEG} ${chromInfoFile} ${PEAKTRACKMINUS}"
    printT "Annotating peaks BIGWIG"
    runScript "annotate_segm31.sh ${OUTBEDGRAPH}.segm"
    runScript "annotate_segm31.sh ${OUTBEDGRAPHOPP}.segm"
    
  fi

else
  # if BAM or FASTQ

  printT "Segmenting [positive strand]"
  runScript "segment_bedgraph_entropy_position2.sh ${OUTBAM}.pos.bedgraph pos"

  printT "Segmenting [negative strand]"
  runScript "segment_bedgraph_entropy_position2.sh ${OUTBAM}.neg.bedgraph neg"

  # add features to segmentation output
  OUTBEDGRAPHPOS=${OUTBAM}.pos.bedgraph.segm
  OUTBEDGRAPHNEG=${OUTBAM}.neg.bedgraph.segm
  bash ${SPARPATH}/scripts/add_features_to_segm.sh ${OUTBEDGRAPHPOS} ${OUTBEDGRAPHNEG}

  numFields=$(awk '{print NF; exit}' ${OUTBAM}.pos.bedgraph.segm)

  printT "Creating track files (Raw signal)"
  echo -e "${chromInfoFile}\t${OUTBAM}"
  runScript "create_bedgraph_track.sh ${OUTBAM}.pos.bedgraph ${chromInfoFile} ${RAWPLUS}"
  runScript "create_bedgraph_track.sh ${OUTBAM}.neg.bedgraph ${chromInfoFile} ${RAWMINUS}"

  printT "Creating track files (Called peaks)"
  runScript "create_peak_track.sh ${OUTBAM}.pos.bedgraph.segm ${chromInfoFile} ${PEAKTRACKPLUS}"
  runScript "create_peak_track.sh ${OUTBAM}.neg.bedgraph.segm ${chromInfoFile} ${PEAKTRACKMINUS}"

  printT "Annotating [positive strand]"
  runScript "annotate_segm31.sh ${OUTBAM}.pos.bedgraph.segm"
  printT "Annotating [negative strand]"
  runScript "annotate_segm31.sh ${OUTBAM}.neg.bedgraph.segm"

fi

# create bigwig segmentation track for JBrowse
#chromInfoFile
printT "${chromInfo}"
for fBG in ${OUTBAM}*.bedGraph; do
  echo "${fBG}"
  fBW=${fBG/.bedgraph.segm.bedGraph/.segm.bigWig}
  if [ -s "${fBG}" ]; then
    if [[ "${fBG}" =~ ([.]pos[.]) ]]; then
      ${BGTOBIGWIG} ${fBG} ${chromInfoFile} ${JBPEAKTRACKPLUS}
    else
      ${BGTOBIGWIG} ${fBG} ${chromInfoFile} ${JBPEAKTRACKMINUS}
    fi
  fi
done




finalAnnot="${OUTBAM}.annot.final.xls"
finalUnannot="${OUTBAM}.unannot.xls"
finalAnnot="${PEAKTABLEANNOT}"
finalUnannot="${PEAKTABLEUNANNOT}"

cat ${annotHeaderFile} > ${finalAnnot}

positiveannot="${OUTBAM}.pos.bedgraph.segm.annot.final.xls"
if [ -s ${positiveannot} ]; then
  cat ${positiveannot} >> ${finalAnnot}
fi
negativeannot="${OUTBAM}.neg.bedgraph.segm.annot.final.xls"
if [ -s "${negativeannot}" ]; then
  cat ${negativeannot} >> ${finalAnnot}
fi

# SORT final annotation
sort -k1,1 -k2,2n -k3,3n -k6,6 ${finalAnnot} -o ${finalAnnot}

# ALL peak table
segmAllBED=${PEAKBEDALL}
cat ${OUTBAM}*.segm | sort -k1,1 -k2,2n -k3,3n -k6,6 > ${segmAllBED} 
segmAllFile=${PEAKTABLEALL}
cat ${segmHeaderFile} ${segmAllBED} > ${segmAllFile}


cat ${annotHeaderFile} > ${finalUnannot}

positiveunannot="${OUTBAM}.pos.bedgraph.segm.unannotated.bed"
if [ -s "${positiveunannot}"  ]; then
  cat ${positiveunannot} >> ${finalUnannot} #${OUTBAM}.unannot
fi

negativeunannot="${OUTBAM}.neg.bedgraph.segm.unannotated.bed"
if [ -s "${negativeunannot}" ]; then
  cat ${negativeunannot} >> ${finalUnannot} #${OUTBAM}.unannot
fi

# SORT final un-annotation
sort -k1,1 -k2,2n -k3,3n -k6,6 ${finalUnannot} -o ${finalUnannot}

baseCoverage=$(awk '{coverage+=($3-$2)}END{print coverage}' ${finalAnnot} ${finalUnannot})


# build small RNA gene expression table
refGeneTable=`dirname $0`/annot/${genomeBuild}/${genomeBuild}.fulltable.no_mRNA_no_lncRNA.unique_LOC.bed
refGeneExprFile=${OUTBAM}.gene_expression.xls
refGeneExprFile=${GENEEXPRTABLE}
>${refGeneExprFile}
# add gene expression table header
echo -e "#Gene\tGeneClass\tReadCount\tRPM" > ${refGeneExprFile}
librarySize=$(awk '{l+=$5}END{print l}' ${finalAnnot} ${finalUnannot})
awk 'BEGIN{OFS="\t"; sep="\t";offset='${numFields}'+0;librarySize='${librarySize}'+0;}
     {
       if (NR==FNR)
       {
         chr = $1; chrStart = $2; chrEnd = $3; geneName=$4; geneClass=$5;
         strand=$6;
         refGeneTable[geneName] = geneClass;
         refGeneExpr[geneName] = 0.0;
       }
       else
       {
         if (NR==1) {next;} # skip header
         exprVal=$5; # expression value (raw read count)
         exprGeneName=$(4+offset); # expressed gene name
         refGeneExpr[exprGeneName]+=exprVal;
       }
       
     }
     END{ for (g in refGeneTable)
            print g, refGeneTable[g], refGeneExpr[g], refGeneExpr[g]/librarySize*1000000;
     }' ${refGeneTable} ${finalAnnot} >> ${refGeneExprFile}



# compute conservation for both annotated and unannotated loci
if [ "${doConservation}" = 1 ]; then

  printT "Computing conservation..."
  bash ${SPARPATH}/scripts/get_conservation.sh ${finalAnnot} ${conservationTrack} > ${finalAnnot/.xls/.with_conservation.xls} && echo "Annotated set finished." || echo "Conservation computation failed (exit code=$?)"

  bash ${SPARPATH}/scripts/get_conservation.sh ${finalUnannot} ${conservationTrack} > ${finalUnannot/.xls/.with_conservation.xls} && echo "Unannotated set finished." || echo "Conservation computation failed (exit code=$?)"
  printT "Done computing conservation"

fi

  printT "Computing overlapping genomic elements..."
  fprefix="peaks"
  printT "Computing overlapping mRNA elements..."
  echo "calculate_mRNA_genomic_partition.sh ${OUTDIR}/results ${fprefix} ${SPARPATH}/annot/partition_files/${genomeBuild}/all_mRNA_annotations_both_strand_with_names.bed.gz" 
  runScript "calculate_mRNA_genomic_partition.sh ${OUTDIR}/results ${fprefix} ${SPARPATH}/annot/partition_files/${genomeBuild}/all_mRNA_annotations_both_strand_with_names.bed.gz" 
  printT "Computing overlapping lncRNA elements..."
  runScript "calculate_lncRNA_genomic_partition.sh ${OUTDIR}/results ${fprefix} ${SPARPATH}/annot/partition_files/${genomeBuild}/all_lncRNA_annotations_both_strand_with_names.bed.gz" 
  printT "Computing overlapping repeat elements..."
  runScript "calculate_repeat_genomic_partition.sh ${OUTDIR}/results ${fprefix} ${SPARPATH}/annot/partition_files/${genomeBuild}/all_repeat_annotations_both_strand_with_names.bed.gz" 
  printT "Computing overlapping elements combining results ..."
  runScript "combine_partition_results.sh ${OUTDIR}/results ${fprefix}" 
  printT "Done computing overlapping genomic elements."

avgMatLen="-1"
if [ "${isBAM}" = 1 ] || [ "${isFASTQ}" = 1 ]; then

  printT "Running sequence analysis"
  bash ${SPARPATH}/scripts/seq_expr.sh ${OUTBAM}.collapsed.bed ${GENOMEFA} 10 > ${OUTBAM}.mature_seq
  bash ${SPARPATH}/scripts/overlap_seq_peaks.sh ${OUTBAM}.mature_seq ${finalAnnot} 0.9 > ${OUTBAM}.annot.mature_seq.raw.xls || printT "annot overlap failed"
  bash ${SPARPATH}/scripts/overlap_seq_peaks.sh ${OUTBAM}.mature_seq ${finalUnannot} 0.9 > ${OUTBAM}.unannot.mature_seq.raw.xls || printT "unannot overlap failed"
  bash ${SPARPATH}/scripts/filter_mature_seq.sh ${OUTBAM}.annot.mature_seq.raw.xls > ${OUTBAM}.annot.mature_seq.xls
  bash ${SPARPATH}/scripts/filter_mature_seq.sh ${OUTBAM}.unannot.mature_seq.raw.xls > ${OUTBAM}.unannot.mature_seq.xls

  matSeqAllFile=${OUTBAM}.all_mature_seq.xls
  matSeqAllFile=${MATSEQALL}

  tail -n+2 ${OUTBAM}.unannot.mature_seq.xls | cat ${OUTBAM}.annot.mature_seq.xls - | cut -f1-13 > ${matSeqAllFile}

  mv ${OUTBAM}.annot.mature_seq.xls ${MATSEQANNOT}
  mv ${OUTBAM}.unannot.mature_seq.xls ${MATSEQUNANNOT}

  avgMatLen=$(awk '{if (NR==1) next; len+=($3-$2);}END{print len/(NR-1)}' ${matSeqAllFile})

fi


# overlap with DASHR

annotFile=${finalAnnot}
unannotFile=${finalUnannot}
DASHRexpr=${SPARDATAPATH}/annot/DASHRv2/${genomeBuild}/DASHR1_GEO_${genomeBuild}.expressionRPM.txt.fixed
DASHRUnAnnotexpr=${SPARDATAPATH}/annot/DASHRv2/${genomeBuild}/DASHR1_GEO_${genomeBuild}.unannot.expressionRPM.txt.fixed

if [ -s ${DASHRexpr} ]; then
  annotDASHRexprFileOut="${annotFile}.dashr.expr.intersect.xls";
  unannotDASHRexprFileOut="${unannotFile}.dashr.expr.intersect.xls";
  combinedDASHRexprFileOut="${OUTBAM}.combined.dashr.expr.intersect.xls";
  noDASHRoverlapFileOut="${OUTBAM}.peaks.noDASHRoverlap.xls";
  printT "Overlapping with DASHR annotated loci"
  bash ${SPARPATH}/scripts/overlap_with_DASHRv2.sh ${annotFile} ${DASHRexpr} || echo "overlap ${annotFile} with DASHR failed"
  printT "Overlapping with DASHR un-annotated loci"
  bash ${SPARPATH}/scripts/overlap_with_DASHRv2.sh ${unannotFile} ${DASHRUnAnnotexpr} || echo "overlap ${unannotFile} with DASHR failed"
  printT "Combining overlaps with DASHR for annotated and un-annotated loci"
  tail -n +2 ${unannotDASHRexprFileOut} | cat ${annotDASHRexprFileOut} - > ${combinedDASHRexprFileOut} || echo "failed to combine overlaps for annotated and un-annotated loci"
  ls -l ${OUTDIR}/results/*.dashr*
  tail -n +2 ${unannotFile}.dashr.peaks_no_overlap | cat ${annotFile}.dashr.peaks_no_overlap - > ${noDASHRoverlapFileOut} || echo "failed to create no_overlap file"

  mv ${annotDASHRexprFileOut} "${DASHRCOMPANNOT}";
  mv ${unannotDASHRexprFileOut} "${DASHRCOMPUNANNOT}";
  mv ${combinedDASHRexprFileOut} "${DASHRCOMPALL}";
  mv ${noDASHRoverlapFileOut} "${NOTINDASHR}";

fi



# overlap with ENCODE

annotFile=${finalAnnot}
unannotFile=${finalUnannot}
ENCODEtable=${SPARPATH}/annot/ENCODE/A1_B2/ENCODE_tier1and2_annot_PEAKS_pavel.unique.csv
ENCODEUnAnnottable=${SPARPATH}/annot/ENCODE/A1_B2/ENCODE_tier1and2_unannot_PEAKS_pavel.unique.csv
ENCODEexpr=${SPARPATH}/annot/ENCODE/ENCODE_tier1and2_annot_PEAKS_pavel.unique.csv.expression_table.minLength10
ENCODEUnAnnotexpr=${SPARPATH}/annot/ENCODE/ENCODE_tier1and2_annot_PEAKS_pavel.unique.csv.expression_table.minLength10
ENCODEexpr=${SPARDATAPATH}/annot/DASHRv2/${genomeBuild}/ENCODE_dataportal_${genomeBuild}.expressionRPM.txt.fixed
ENCODEUnAnnotexpr=${SPARDATAPATH}/annot/DASHRv2/${genomeBuild}/ENCODE_dataportal_${genomeBuild}.unannot.expressionRPM.txt.fixed
# only loci found in at least 2 tissues
ENCODEUnAnnotexpr=${SPARDATAPATH}/annot/DASHRv2/${genomeBuild}/ENCODE_dataportal_${genomeBuild}.unannot.expressionRPM.txt.fixed.atleast2
if [ -s ${ENCODEexpr} ]; then
  annotENCODEexprFileOut="${annotFile}.encode.expr.intersect.xls";
  unannotENCODEexprFileOut="${unannotFile}.encode.expr.intersect.xls";
  combinedENCODEexprFileOut="${OUTBAM}.combined.encode.expr.intersect.xls";
  noENCODEoverlapFileOut="${OUTBAM}.peaks.noENCODEoverlap.xls";
  printT "Overlapping with ENCODE annotated loci"
  bash ${SPARPATH}/scripts/overlap_with_DASHRv2_encode.sh ${annotFile} ${ENCODEexpr} || echo "overlap ${annotFile} with ENCODE failed"

  printT "Overlapping with ENCODE un-annotated loci"
  bash ${SPARPATH}/scripts/overlap_with_DASHRv2_encode.sh ${unannotFile} ${ENCODEUnAnnotexpr} || echo "overlap ${unannotFile} with ENCODE failed"

  printT "Combining overlaps with ENCODE for annotated and un-annotated loci"
  tail -n +2 ${unannotENCODEexprFileOut} | cat ${annotENCODEexprFileOut} - > ${combinedENCODEexprFileOut} || echo "failed to combine overlaps for annotated and un-annotated loci"
  cat ${OUTDIR}/results/*encode.peaks_no_overlap > ${noENCODEoverlapFileOut} || echo "failed to create no_overlap file"

  mv ${annotENCODEexprFileOut} "${ENCODECOMPANNOT}";
  mv ${unannotENCODEexprFileOut} "${ENCODECOMPUNANNOT}";
  mv ${combinedENCODEexprFileOut} "${ENCODECOMPALL}";
  mv ${noENCODEoverlapFileOut} "${NOTINENCODE}";

fi

printT "DONE."

# Split annotated peak tables by class

for segmFile in ${OUTBAM}*.segm; do
  SEGMFILE=${segmFile}
  break
done
numFields=$(awk '{print NF; exit}' ${SEGMFILE})
annotHeader=$(awk '{gsub(/\t/,","); gsub(/#/,""); print}' ${annotHeaderFile})
printT "SEGMFILE=${SEGMFILE}; numFields=${numFields}"
awk 'BEGIN{ numFields='${numFields}'+0; # number of fields after segm
            #filePrefix=FILENAME; gsub("/.xls$/","",filePrefix)
            filePrefix="'${OUTDIR}'/annot";
            filePrefix="'${OUTDIR}'/results/peaks";
            filePrefix="'${OUTDIR}'/results/byClass/peaks";
            header="#'${annotHeader}'";
            gsub(",","\t",header);
          }
     {
       if (NR==1) next; # skip header
       rnaClass=$(numFields+5);
       outfile = (filePrefix ".byClass." rnaClass ".xls");

       rnaClassCnt[rnaClass]++;
       if (rnaClassCnt[rnaClass]==1)
          print header > outfile
       #print outfile
       
       print > outfile;
     }' ${finalAnnot}

# annotation summary
annotSummary=${OUTBAM}.mapped_reads_annotation_summary.txt
annotSummary=${OUTDIR}/results/annotation_summary.txt
awk 'BEGIN{OFS="\t";totalAnnotPeakCnt=0; totalUnannotPeakCnt=0;totalExprAnnot=0;totalExprUnannot=0;numFields='${numFields}'+0;}{if ($0~/^#/) {next}; if (NR==FNR) {geneID=$(numFields+4); rnaClass=$(numFields+5);  n=split(geneID,a,":"); gene=a[n]; gene=geneID; if (ids[gene]!=1) {ids[gene]=1; perClassGeneCnt[rnaClass]++;totalAnnotGeneCnt++;}; exprVal=$5; exprPerClass[rnaClass]+=exprVal;classCnt[rnaClass]+=1;totalExprAnnot+=exprVal; totalAnnotPeakCnt+=1}else{exprVal=$5; totalExprUnannot+=exprVal;totalUnannotPeakCnt+=1;}}END{totalPeakCnt=totalAnnotPeakCnt+totalUnannotPeakCnt; totalExpr=totalExprAnnot+totalExprUnannot; for (rnaClass in exprPerClass) printf "%s\t%d\t%d\t%d\t%.2f\n", rnaClass, classCnt[rnaClass], perClassGeneCnt[rnaClass]+0, exprPerClass[rnaClass], 100*exprPerClass[rnaClass]/totalExpr; propAnnot=0; if (totalExpr>0) propAnnot=totalExprAnnot/totalExpr; printf "%s\t%d\t%d\t%d\t%.2f\n", "Annotated", totalAnnotPeakCnt, totalAnnotGeneCnt, totalExprAnnot, 100*propAnnot; propUnannot=0; if (totalExpr>0) propUnannot=totalExprUnannot/totalExpr; printf "%s\t%d\t%d\t%d\t%.2f\n", "Unannotated",totalUnannotPeakCnt,totalUnannotPeakCnt,totalExprUnannot,100*propUnannot}' ${finalAnnot} ${finalUnannot} | sort -k1,1 | awk 'BEGIN{OFS="\t"; print "#RNA class","Peaks","Genes", "Reads","Percentage of reads"}{print}' > ${annotSummary} 

chartLabelsFile=${OUTDIR}/results/chart.labels
chartReadsFile=${OUTDIR}/results/chart.reads.data
chartPeaksFile=${OUTDIR}/results/chart.peaks.data
LC_ALL=en_US.UTF-8 awk 'BEGIN{FS="\t";
           chartLabelsFile="'${chartLabelsFile}'";
           chartReadsFile="'${chartReadsFile}'";
           chartPeaksFile="'${chartPeaksFile}'";
     }
     { class[NR]=$1; peaks[NR]=$2; reads[NR]=$4;}
     END{ n = NR;
          for (i=3; i<=n; ++i)
          {
            printf "\"%s\"", class[i] > chartLabelsFile;
            if (i<n) printf "," > chartLabelsFile;
            printf "%s", reads[i] > chartReadsFile;
            if (i<n) printf "," > chartReadsFile;
            printf "%s", peaks[i] > chartPeaksFile;
            if (i<n) printf "," > chartPeaksFile;
            totalPeaks+=peaks[i];
            totalReads+=reads[i];
          }
          printf "\n" > chartLabelsFile
          printf "\n" > chartReadsFile
          printf "\n" > chartPeaksFile
          printf "%\047d", totalPeaks > (chartPeaksFile ".sum") 
          printf "%\047d", totalReads > (chartReadsFile ".sum") 
     }' ${annotSummary}

# reads chart
cat "${chartReadsFile}" "${chartLabelsFile}" "${chartReadsFile}.sum" | \
  awk '{if (FNR==NR) {data[FNR]=$0;} else { gsub("numbers_here", ("[" data[1] "]")); gsub("labels_here", ("[" data[2] "]")); gsub("sum_here", ("\x27" data[3]  "\x27") );  print; }}' - ${SPARPATH}/templates/chart1.download.html > ${OUTDIR}/figures/chart_reads_download.html 

cat "${chartReadsFile}" "${chartLabelsFile}" "${chartReadsFile}.sum" | \
  awk '{if (FNR==NR) {data[FNR]=$0;} else { gsub("numbers_here", ("[" data[1] "]")); gsub("labels_here", ("[" data[2] "]")); gsub("sum_here", ("\x27" data[3]  "\x27") );  print; }}' - ${SPARPATH}/templates/chart1.html > ${OUTDIR}/figures/chart_reads.html 
# peaks chart
cat "${chartPeaksFile}" "${chartLabelsFile}" "${chartPeaksFile}.sum" | \
  awk '{if (FNR==NR) {data[FNR]=$0;} else { gsub("numbers_here", ("[" data[1] "]")); gsub("labels_here", ("[" data[2] "]")); gsub("sum_here", ("\x27" data[3]  "\x27") );  print; }}' - ${SPARPATH}/templates/chart2.download.html > ${OUTDIR}/figures/chart_peaks_download.html 

cat "${chartPeaksFile}" "${chartLabelsFile}" "${chartPeaksFile}.sum" | \
  awk '{if (FNR==NR) {data[FNR]=$0;} else { gsub("numbers_here", ("[" data[1] "]")); gsub("labels_here", ("[" data[2] "]")); gsub("sum_here", ("\x27" data[3]  "\x27") );  print; }}' - ${SPARPATH}/templates/chart2.html > ${OUTDIR}/figures/chart_peaks.html 

if [ 1 -eq 0 ]; then
awk 'BEGIN{FS="\t"; OFS="\t";numFields='${numFields}'+0;}
     {
       if (NR==1) next; # skip header
       annotID=$(numFields+4);
       #print annotID;
       n=split(annotID,a,":");
       geneID=a[n];
       annotClass=$(numFields+5);
       if (ids[geneID]!=1)
       {
          geneCnt[annotClass]+=1;
          ids[geneID]=1;
       }
     }
     END{ for (c in geneCnt)
           print c, geneCnt[c]
     }' ${finalAnnot} | sort -k1,1 | awk 'BEGIN{OFS="\t"; FS="\t"; print "#RNA class", "Genes"}{print; ngenes+=$2}END{print "Total", ngenes}' >> ${annotSummary}
fi

annotLengthSummary=${OUTBAM}.annot.peak.length.stats
>${annotLengthSummary}
#echo -e "\nLength of annotated peaks:\n" >> ${annotLengthSummary}
awk 'BEGIN{FS="\t"; OFS="\t"}
     {
       if (NR==1) next;
       peakStart = $2; # 0-based
       peakEnd = $3; # 1-based according to 0-based, half-open UCSC notation
       peakLength = peakEnd-peakStart;
       lengthCnt[peakLength]+=1;  
       totalPeakCnt+=1;
     }
     END{ for (l in lengthCnt)
          {
            lengthProp = lengthCnt[l]/totalPeakCnt
            print l, lengthCnt[l], lengthProp;
          }
     }' ${finalAnnot} | sort -k1,1n | awk 'BEGIN{FS="\t"; OFS="\t"; print "#PeakLength", "Count", "Fraction"}{print}' > ${annotLengthSummary}

unannotLengthSummary=${finalUnannot}.peak.length.stats
>${unannotLengthSummary}
#echo -e "\nLength of annotated peaks:\n" >> ${annotLengthSummary}
awk 'BEGIN{FS="\t"; OFS="\t"}
     {
       if (NR==1) next;
       peakStart = $2; # 0-based
       peakEnd = $3; # 1-based according to 0-based, half-open UCSC notation
       peakLength = peakEnd-peakStart;
       lengthCnt[peakLength]+=1;  
       totalPeakCnt+=1;
     }
     END{ for (l in lengthCnt)
          {
            lengthProp = lengthCnt[l]/totalPeakCnt
            print l, lengthCnt[l], lengthProp;
          }
     }' ${finalUnannot} | sort -k1,1n | awk 'BEGIN{FS="\t"; OFS="\t"; print "#PeakLength", "Count", "Fraction"}{print}' > ${unannotLengthSummary}





printL "\n===Output==="
printL "Output directory: ${OUTDIR}"
printL "LOG file: ${LOGSPAR}"
if [ "${isFASTQ}" = 1 ]; then
  printL "\nMapping output:"
  printL "${OUTBAM}"
fi

if [ "${isBAM}" = 1 ] || [ "${isFASTQ}" = 1 ]; then
  printL "\nTrack files (Raw signal):"
  ls ${OUTDIR}/tracks/raw*.bigWig | tee -a ${LOGSPAR}
fi
printL "\nTrack files (Called peaks):"
ls ${OUTDIR}/tracks/*.bigBed | tee -a ${LOGSPAR}

printL "\nAnnotation output:"
ls ${finalAnnot} | tee -a ${LOGSPAR}

printL "\nUn-annotated output:" 
ls ${finalUnannot} | tee -a ${LOGSPAR}


printL "\nAnnotation summary:"
ls ${annotSummary}  | tee -a ${LOGSPAR}

if [ "${isFASTQ}" = 1 ]; then 
  printL "\nMapping stats:"
  ls ${OUTDIR}/mapping/MAPSTAT.txt | tee -a ${LOGSPAR}
  printL "\nMapped reads alignment summary:"
  ls ${OUTDIR}/mapping/cigar.stat | tee -a ${LOGSPAR}
fi



printL "\n\n===Run summary==="

if [ "${isFASTQ}" = 1 ]; then 
  printL "FASTQ: ${INFILE}"
  grep -e "FASTQ reads" ${OUTDIR}/mapping/MAPSTAT.txt | awk 'BEGIN{FS="\t"}{printf "Total reads: %d [%.4f%%]\n", $2, $3}' | tee -a ${LOGSPAR}
  grep -e "Reads \[all\]" ${OUTDIR}/mapping/MAPSTAT.txt | awk 'BEGIN{FS="\t"}{printf "Reads after QC: %d [%.4f%%]\n", $2, $3}' | tee -a ${LOGSPAR}
 echo "Check ${OUTDIR}/mapping/MAPSTATS.txt for details"
elif [ "${isBAM}" = 1 ]; then
  printL "BAM: ${INFILE}"
elif [ "${isBIGWIG}" = 1 ]; then
  printL "BIGWIG: ${INFILE}"
fi

numAnnot=$(wc -l ${finalAnnot} | awk '{print $1-1}')
numUnannot=$(wc -l ${finalUnannot} | awk '{print $1-1}')
printL "Annotated loci count: ${numAnnot}"

printL "Un-annotated loci count: ${numUnannot}"


printL "\nAnnotation summary:"
ls ${annotSummary} | tee -a ${LOGSPAR}

printL "`cat ${annotSummary}`"

plot_script=${OUTDIR}/spar_plots.sh
>${plot_script}

cat "${SPARPATH}/scripts/debug_header" > ${plot_script}

<<SKIP_MODULE1
if [  "${isBAM}" = 1 ] || [ "${isFASTQ}" = 1 ]; then
  # insilico cut stats
  module="module1"
  step="insilico_cut"
  outRscript="${OUTDIR}/${module}_${step}.r"
  fprefix=`basename "${INFILE}"`
  fprefix="peaks";
  cat ${SPARPATH}/scripts/R/${module}/${step}/*.r | \
     awk 'BEGIN{wdir="wdir=\"'${OUTDIR}'\""; print wdir; fprefix="fprefix=\"'${fprefix}'\""; print fprefix; maxReadLength="maxReadLength='${maxReadLength}'"; print maxReadLength;}{print}' - > ${outRscript}
  
  
  submoduleDir=${SPARPATH}/scripts/R/${module}/${step}

  makehtmlcmd="make_submodule_html.sh ${submoduleDir} . > ${OUTDIR}/${module}_${step}.html"
  echo "bash ${SPARPATH}/scripts/${makehtmlcmd}" >> ${plot_script}
  echo "${RSCRIPT} ${outRscript} &>> ${RLOG}" >> ${plot_script}
fi 

if [ "${isFASTQ}" = 1 ]; then

  printT "Creating plots (module 1)"
  # mapping stats
  module="module1"
  step="mapping"
  outRscript="${OUTDIR}/${module}_${step}.r"
  submoduleDir=${SPARPATH}/scripts/R/${module}/${step}
  cat ${submoduleDir}/*.r | \
     awk 'BEGIN{wdir="wdir=\"'${OUTDIR}'\""; print wdir;}{print}' - > ${outRscript}
  makehtmlcmd="make_submodule_html.sh ${submoduleDir} . > ${OUTDIR}/${module}_${step}.html"
  echo "bash ${SPARPATH}/scripts/${makehtmlcmd}" >> ${plot_script}
  echo "${RSCRIPT} ${outRscript} &>> ${RLOG}" >> ${plot_script}
fi
SKIP_MODULE1

  module="module3"
  step="annot"
  outRscript="${OUTDIR}/${module}_${step}.r"
  outRscriptPar="${OUTDIR}/${module}_${step}_PAR.r"
  fprefix=`basename "${OUTBAM}"`
  fprefix="peaks"
  submoduleDir=${SPARPATH}/scripts/R/${module} #/${step}
  genomicAnnot="${SPARPATH}/annot/partition_files/${genomeBuild}"
  sparPath="${SPARPATH}"
  cat ${submoduleDir}/*.r | \
     awk 'BEGIN{wdir="wdir=\"'${OUTDIR}'\""; print wdir; fprefix="fprefix=\"'${fprefix}'\""; print fprefix; partition_ref_path="partition_ref_path=\"'${genomicAnnot}'\""; print partition_ref_path;}{print}' - > ${outRscript}

  cat ${submoduleDir}/parallel_template.r | \
     awk 'BEGIN{wdir="wdir=\"'${OUTDIR}'/results\""; print wdir; fprefix="fprefix=\"'${fprefix}'\""; print fprefix; partition_ref_path="partition_ref_path=\"'${genomicAnnot}'\""; print partition_ref_path; spar_path="SPAR_path=\"'${sparPath}'\""; print spar_path;}{print}' - > ${outRscriptPar}

  printT "Creating plots"
  printL "Running R scripts for plotting"
  printL "RSCRIPT=${RSCRIPT}\n${outRscript}\n${submoduleDir}\n"
  rcmd="/home/pkuksa/bin/R-3.2.3/bin/Rscript ${outRscript} 2>&1"
  rcmd="${RSCRIPT} ${outRscript} 2>&1"
  rcmd="${RSCRIPT} ${outRscript} >> ${RLOG}"
  rcmdPar="/home/pkuksa/bin/R-3.2.3/bin/Rscript ${outRscriptPar} 2>&1"
  rcmdPar="${RSCRIPT} ${outRscriptPar} 2>&1"
  rcmdPar="${RSCRIPT} ${outRscriptPar} >> ${RLOG}"
  cat "${SPARPATH}/scripts/debug_header" > ${OUTDIR}/make_plots.sh
  echo "${rcmd}" >> ${OUTDIR}/make_plots.sh
  echo "${rcmdPar}" >> ${OUTDIR}/make_plots_parallel.sh
  makehtmlcmd="make_submodule_html.sh ${submoduleDir} ${OUTDIR}/figures > ${OUTDIR}/${module}_${step}.html"
  echo -e "bash ${SPARPATH}/scripts/${makehtmlcmd}" >> ${OUTDIR}/make_plots.sh
  echo -e "bash ${SPARPATH}/scripts/${makehtmlcmd}" >> ${OUTDIR}/make_plots_parallel.sh
  echo -e "exit 0" >> ${OUTDIR}/make_plots_parallel.sh

  echo "bash ${OUTDIR}/make_plots_parallel.sh" >> ${plot_script}


  cp ${finalAnnot}  ${OUTDIR}/results/byClass/peaks.byClass.Annotated.xls
  cp ${finalUnannot}  ${OUTDIR}/results/byClass/peaks.byClass.Unannotated.xls

LC_ALL=en_US.UTF-8 awk 'BEGIN{OFS="<td>";FS="\t"}{if ($0~/^#/) {if (NR>1) print "</table>\n<p>&nbsp;"; print "<table>"; inHeader=1; sep="<th>"; gsub(/^#/,"",$0)}else{inHeader=0; sep="<td>"};printf "<tr>"; if (inHeader==0) {rnaClass=$1; $1=("<a href=\"" "'${OUTDIR}'/results/byClass/peaks.byClass." rnaClass ".xls" "\" download>" rnaClass "</a>")}; printf "%s%s", sep, $1; if (inHeader==0) sep="<td style=\"text-align:right\">"; for (i=2;i<NF;i++) if (inHeader==0) printf ("%s%\047d", sep, $i+0); else {printf "%s%s", sep, $i}; printf "%s%s", sep, $NF; printf "\n";}END{print "</table>"}' ${annotSummary} > ${annotSummary}.html

annotTable=${finalAnnot/.xls/.with_conservation.xls}
unannotTable=${finalUnannot/.xls/.with_conservation.xls}
HTMLPEAKTABLE=${OUTBAM}.peak_table.html
HTMLPEAKTABLE=${OUTDIR}/results/peak_browser.html


# USE FULL TABLE
awk '{if (NR==1) gsub(/^X./,"#"); print}' "${OUTDIR}/results/Genomewide_distribution_patterns_of_small_RNA_loci_all_partition_summary_per_unannot_locus.txt" > "${OUTDIR}/results/peaks_unannot.xls"
awk '{if (NR==1) gsub(/^X./,"#"); print}' "${OUTDIR}/results/Genomewide_distribution_patterns_of_small_RNA_loci_all_partition_summary_per_annot_locus.txt" > "${OUTDIR}/results/peaks_annot.xls"

sort -k1,1 -k2,2n -k3,3n -k6,6 ${finalAnnot} -o ${finalAnnot}
sort -k1,1 -k2,2n -k3,3n -k6,6 ${finalUnannot} -o ${finalUnannot}

annotTable=${finalAnnot}
unannotTable=${finalUnannot}

bash ${SPARPATH}/scripts/make_html_peak_table.sh ${annotTable} ${unannotTable} ${GENOMEBUILD}  > ${HTMLPEAKTABLE}



runSummaryHtml=${OUTDIR}/run_summary.html
>${runSummaryHtml}
LC_ALL=en_US.UTF-8
totalAnnotReads=$(grep -e "Annotated" ${annotSummary} | awk '{print $4}' )
totalUnannotReads=$(grep -e "Unannotated" ${annotSummary} | awk '{print $4}' )
totalReads=$((totalAnnotReads+totalUnannotReads))

read avgSegmLen avgSegmExpr nSegm <<< $(awk 'BEGIN{FS="\t"}{if (NR==1) next; len=$3-$2; exprRPM=$(NF-1); l+=len; e+=exprRPM; n+=1;}END{if (n>0) {print l/n, e/n, n} else {print "0", "0", "0";}}' ${segmAllFile})



#numAnnotGenes=$(grep -e "Total" ${annotSummary} | awk '{print $2}')
numAnnotGenes=$(grep -e "^Annotated" ${annotSummary} | awk '{print $3}')
echo "<table>" >> ${runSummaryHtml}
printf "<tr><td>Reads<td class=\"alnright\">%'d\n" ${totalReads} >> ${runSummaryHtml}
printf "<tr><td>Expressed small RNA loci<td class=\"alnright\">%'d\n" ${nSegm} >> ${runSummaryHtml}
printf "<tr><td>Reads (annotated)<td class=\"alnright\">%'d\n" ${totalAnnotReads} >> ${runSummaryHtml}
printf "<tr><td>Reads (unannotated)<td class=\"alnright\">%'d\n" ${totalUnannotReads} >> ${runSummaryHtml}
printf "<tr><td>Genes (annotated)<td class=\"alnright\">%'d\n" ${numAnnotGenes} >> ${runSummaryHtml}
printf "<tr><td>Called peaks (annotated)<td class=\"alnright\">%'d\n" ${numAnnot} >> ${runSummaryHtml}
printf "<tr><td>Called peaks (unannotated)<td class=\"alnright\">%'d\n" ${numUnannot} >> ${runSummaryHtml}
printf "<tr><td>Expressed loci length (average)<td class=\"alnright\">%'.2f\n" ${avgSegmLen}  >> ${runSummaryHtml}
printf "<tr><td>Genome coverage (nucleotides)<td class=\"alnright\">%'d" ${baseCoverage} >> ${runSummaryHtml}
if [ "${avgMatLen}" != "-1" ]; then
  printf "<tr><td>Mature product size (average)<td class=\"alnright\">%'.2f\n" ${avgMatLen} >> ${runSummaryHtml}
fi
printf "<tr><td>RPM (average)<td class=\"alnright\">%'.2f" ${avgSegmExpr} >> ${runSummaryHtml}
echo "</table>" >> ${runSummaryHtml}

ENDTIME=$(date +%s)
totalTime=$((ENDTIME-STARTTIME))
processingSpeed=$(LC_ALL=en_US.UTF-8 awk 'BEGIN{printf "Processed %\047d reads in %\047d seconds (%\047d reads / second)\n", '${totalReads}', '${totalTime}', '${totalReads}'/('${totalTime}'+0.0); exit}')
printT "${processingSpeed}"

# get file sizes

for v in ${!PEAK*} ${!RAW*} ${!GENE*} ${!MATSEQ*} ${!DASHRCOMP*} ${!ENCODECOMP*} ${!NOTIN*}; do
  if [ -s ${!v} ]; then
   bytesize=`du -b ${!v} | cut -f1`;
   size=$(echo "${bytesize}" | \
   awk '{ size = $1;
         sizepower = int( (length(size)-1)/3 );
         if (sizepower < 0 ) sizepower=0;
         dispSize = size / 1024**sizepower;
         sizesString="B,KB,MB,GB,TB,PB";
         split(sizesString,sizes,",");
         #printf "%s %.2f %s", size, dispSize, sizes[sizepower+1];
         printf "%.2f %s", dispSize, sizes[sizepower+1];
       }')
    declare "size$v=$size";
    varname=size_$v;
    echo "${!v} ${!varname}";
  fi
done

# create html table / index for main results
cat ${resultsTemplate} | \
while read -r line; do
    while [[ "$line" =~ ([{][a-zA-Z]*[}]) ]] ; do
        LHS=${BASH_REMATCH[1]}
        RHS="$(eval echo "\"\$${LHS}\"")"
        line=${line//\$${LHS}/${RHS}}
    done
    echo "$line"
done > ${HTMLOUTPUT} 

awk 'BEGIN{d="'${OUTDIR}'/";}{gsub(d,"",$0); print $0;}' ${HTMLOUTPUT} > ${HTMLOUTPUT/.html/.local.html}


# clean-up
rm ${OUTDIR}/*.segm.* || true
rm ${OUTDIR}/mapping/*.segm.* || true
rm ${OUTDIR}/results/*.xls.*
rm ${OUTDIR}/*.tmp || true
rm ${OUTDIR}/mapping/*.tmp || true
rm "${OUTDIR}/module3_annot.html" || true

# clean-up
if [ "$doCleanup" = "1" ]; then
  rm "${OUTDIR}"/raw.pos.* || true
  rm "${OUTDIR}"/raw.neg.* || true
  rm "${OUTDIR}"/*.stats  || true
  #mv "${OUTDIR}/results.local.html" "${OUTDIR}/results.html" || true
  cp $configFile ${OUTDIR}/SPAR.$genomeBuild.config
  #rm "${OUTDIR}"/*.tmp.sh "${OUTDIR}/*.tmp.r"
  #rm "${OUTDIR}"/results/*.with_conservation* || true
  #rm "${OUTDIR}"/results/Genomewide*txt || true
fi


exit 0
