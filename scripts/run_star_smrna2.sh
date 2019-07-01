# map smRNA reads using STAR
# Usage:
# bash run_star_smRNA.sh <full path to sample FASTQ file>
# Example:
# bash run_star_smrna.sh /mnt/niagads/users/yyee/tissues/hsa_adpilot_smrna/hsa_adpilot_ctl02058_smrna/1.5/hsa_adpilot_ctl02058_smrna_trimmed.fastq

set -e

#source `dirname $0`/../config.sh

OUTDIR=$5
LOGSPAR=${OUTDIR}/SPAR.log
function printT
{
  echo "`date +'%b %d %H:%M:%S'` ..... $1" | tee -a ${LOGSPAR}
}


FASTQ=$1 # input FASTQ file
         # <tissue>_<SRR>_trimmed.fastq
maxMismatchCnt=$2 #${maxMismatchCnt}
multiMapMax=$3 #${maxMapCnt}
threads=$4 #${THREADS}
starOutDir=$5

# STAR output directory (absolute path)
mkdir -p ${starOutDir}


# run STAR
printT "Mapping reads"
# 03 Aug 2017: handling compressed inputs
# --readFilesCommand
#                      default: -
#                      string(s):  command line to execute for each of the input file.  This command should generate FASTA or FASTQ text and send it to stdout
#                      For example:  zcat - to uncompress .gz files, bzcat - to uncompress .bz2 files, etc

readFilesCommand=$( echo "${FASTQ}" |
                      awk '{
                            cmd="-";
                            if ($0 ~ /\.gz$/) cmd="zcat";
                            if ($0 ~ /\.bz2$/) cmd="bzcat";
                            print cmd;
                            exit;
                           }' )
if [ "${readFilesCommand}" == "-" ]; then
  ${STAR} --genomeDir ${genomeDir} --genomeLoad LoadAndKeep  --readFilesIn ${FASTQ} --runThreadN $threads --alignIntronMax 1 --outSAMattributes NH HI NM MD --outFilterMultimapNmax ${multiMapMax} --outReadsUnmapped Fastx --outFilterMismatchNmax ${maxMismatchCnt} --outFilterMatchNmin ${minMappedLength} --outFileNamePrefix ${starOutDir}/
else
  ${STAR} --genomeDir ${genomeDir} --genomeLoad LoadAndKeep  --readFilesIn ${FASTQ} --readFilesCommand "${readFilesCommand}" --runThreadN $threads --alignIntronMax 1 --outSAMattributes NH HI NM MD --outFilterMultimapNmax ${multiMapMax} --outReadsUnmapped Fastx --outFilterMismatchNmax ${maxMismatchCnt} --outFilterMatchNmin ${minMappedLength} --outFileNamePrefix ${starOutDir}/
fi

#${STAR} --genomeDir ${genomeDir} --genomeLoad LoadAndKeep  --readFilesIn ${FASTQ} --runThreadN $threads --alignIntronMax 1 --outSAMattributes NH HI NM MD --outFilterMultimapNmax ${multiMapMax} --outReadsUnmapped Fastx --outFilterMismatchNmax ${maxMismatchCnt} --outFilterMatchNmin ${minMappedLength} --outFileNamePrefix ${starOutDir}/

nFastQReads=$(grep -e "Number of input reads" ${starOutDir}/Log.final.out | cut -f 2)

# filter out reads trimmed at 5p end by more than one base
printT "Removing reads trimmed at 5p"
awk 'BEGIN{OFS="\t"; filtPos=0; filtNeg=0; alnCnt=0; totalAln=0; nFastQReads='${nFastQReads}'; readId_prev="";}
     {
       if (/^@/) {print; next}; # SAM header
       readId=$1;
       samFlag=$2; cigar=$6; strandFlag=and(samFlag,16);
       totalAln+=1
       if (strandFlag>0) # negative strand
       {
         if (cigar~/[0-9]S$/)
         {
           n=split(cigar,L,/[NMSID]/);
           clip=L[n-1]+0;
           if (clip>1) {filtNeg+=1; next}
         }
       }
       else # positive strand
       {
         if (cigar~/^[0-9]+S/)
         {
           split(cigar,L,/[NMSID]/);
           clip=L[1];
           if (clip>1) {filtPos+=1; next}
         }
       }
       #a[readId]++;
       if (cigar !~ /S/)
       {
          #if (a[readId]==1) readNoScnt++;
          if (readId != readId_prev) readNoScnt++;
          noSoftClipCnt++;
       }
       else
       {
          #if (a[readId]==1) readScnt++;
          if (readId != readId_prev) readScnt++;
          softClipCnt++;
       }
       alnCnt+=1;
       readId_prev = readId;
       print
     }
     END {
      filter5pStatsFile=(FILENAME ".filter_5p_stats");
      softClipStatsFile=(FILENAME ".softClipStats");

      filt5pPosPct = 0.00;
      filt5pNegPct = 0.00;
      alnCntPct = 0.00;
      if (totalAln>0)
      {
         filt5pPosPct = filtPos / totalAln * 100.00;
         filt5pNegPct = filtNeg / totalAln * 100.00;
         alnCntPct = alnCnt / totalAln * 100.00;
      }
      print "Alignments before removing 5p clipping", totalAln, "100.00" > filter5pStatsFile
      print "Removed 5p clipped [positive strand]:",filtPos,filt5pPosPct > filter5pStatsFile;
      print "Removed 5p clipped [negative strand]:", filtNeg, filt5pNegPct > filter5pStatsFile;
      print "Total after removing 5p clipped:",alnCnt,alnCntPct > filter5pStatsFile;
      noSpct = 0.00;
      withSpct = 0.00;
      if (alnCnt > 0)
      {
        noSpct = noSoftClipCnt / alnCnt * 100.00
        withSpct = softClipCnt / alnCnt * 100.00
      }
      print "Alignments [no soft clipping]", noSoftClipCnt, noSpct > softClipStatsFile;
      print "Alignments [with soft clipping]", softClipCnt, withSpct > softClipStatsFile;
      readNoSpct = 0.00;
      readWithSpct = 0.00;
      if (nFastQReads>0)
      {
         readNoSpct = readNoScnt / nFastQReads * 100.00
         readWithSpct = readScnt / nFastQReads * 100.00
      }
      print "Reads [no soft clipping]", readNoScnt, readNoSpct > softClipStatsFile;
      print "Reads [with soft clipping]", readScnt, readWithSpct > softClipStatsFile;
         }' ${starOutDir}/Aligned.out.sam > ${starOutDir}/Aligned.out.filtered.sam

# convert filtered SAM to BAM file format
printT "Sorting filtered SAM"
${SAMTOOLS} sort -@ $threads -O bam ${starOutDir}/Aligned.out.filtered.sam -T ${starOutDir}/Aligned.out.filtered.sort > ${starOutDir}/Aligned.out.filtered.sorted.bam

# remove soft-clipped reads
#printT "Removing soft-clipped reads"
#${SAMTOOLS} view -h ${starOutDir}/Aligned.out.filtered.sorted.bam | awk 'BEGIN{OFS="\t"}{ if ($0~/^@/) {print; next} if ($6~/S/) {softClipCnt+=1; next;} noSoftClipCnt+=1; print }END{print "Alignments [no soft clipping]", noSoftClipCnt > "'${starOutDir}'/softClipStats.txt"; print "Alignments [with soft clipping]", softClipCnt > "'${starOutDir}'/softClipStats.txt"}' | ${SAMTOOLS} view -bS - > ${starOutDir}/Aligned.out.filtered.noS.sorted.bam


# hard-clip filtered sorted BAM file
printT "Hard-clipping filtered sorted BAM"
#${SAMTOOLS} view -h -@ $threads ${starOutDir}/Aligned.out.filtered.sorted.bam | awk 'BEGIN {OFS="\t"}{if ($0~/^@/) {print; next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' | ${SAMTOOLS} view -@ 4 -bS - > ${starOutDir}/Aligned.out.filtered.hardClipped.sorted.bam

STATCIGAR=${starOutDir}/cigar.stat
>${STATCIGAR}
${SAMTOOLS} view -h -@ $threads ${starOutDir}/Aligned.out.filtered.sorted.bam | awk 'BEGIN {OFS="\t"}{if ($0~/^@/) {print; next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print; cigar=$6; cigarCnt[cigar]++;}END{for (c in cigarCnt) print c, cigarCnt[c] > "'${STATCIGAR}'"}' | ${SAMTOOLS} view -@ 4 -bS - > ${starOutDir}/Aligned.out.filtered.hardClipped.sorted.bam
sort -k2,2nr ${STATCIGAR} -o ${STATCIGAR}

# index BAM
printT "Indexing BAM"
${SAMTOOLS} index ${starOutDir}/Aligned.out.filtered.hardClipped.sorted.bam

#${SAMTOOLS} index ${starOutDir}/Aligned.out.filtered.noS.sorted.bam

${SAMTOOLS} index ${starOutDir}/Aligned.out.filtered.sorted.bam

# compute stats
printT "Computing mapping stats"
statFile=${starOutDir}/MAPSTAT.txt
>${statFile}
nFastQReads=$(grep -e "Number of input reads" ${starOutDir}/Log.final.out | cut -f 2)
nUniqMapReads=$(grep -e "Uniquely mapped reads number" ${starOutDir}/Log.final.out | cut -f 2)
nMultiMapReads=$(grep -e "Number of reads mapped to multiple loci" ${starOutDir}/Log.final.out | cut -f 2)


echo "${nFastQReads} ${nUniqMapReads} ${nMultiMapReads}" | \
   awk 'BEGIN{OFS="\t"}
        {
         print "FASTQ reads:", $1, $1/$1*100
         print "Uniquely-mapped reads:", $2, $2/$1*100
         print "Multi-mapped reads:", $3, $3/$1*100
         print "Total mapped reads:", $2+$3, ($2+$3)/$1*100
        }' >> ${statFile}
cat ${starOutDir}/*.filter_5p_stats ${starOutDir}/*.softClipStats >> ${statFile}

#nNoSReads=$(${SAMTOOLS} view ${starOutDir}/Aligned.out.filtered.noS.sorted.bam | cut -f 1 | sort -u | wc -l)
#nSandNoSReads=$(${SAMTOOLS} view ${starOutDir}/Aligned.out.filtered.sorted.bam | cut -f 1 | sort -u | wc -l)
nNoSReads=$(grep -e "Reads \[no soft" ${starOutDir}/*.softClipStats | cut -f 2)
nWithSReads=$(grep -e "Reads \[with soft" ${starOutDir}/*.softClipStats | cut -f 2)
echo -e "${nNoSReads} ${nWithSReads}" | \
awk 'BEGIN{OFS="\t"}
     {
       #print "Reads [no soft clipping]:", $1, ($1+0.0)/'${nFastQReads}'*100
       #print "Reads [with soft clipping]:", $2, ($2+0.0)/'${nFastQReads}'*100
       print "Reads [all]:", ($1+$2), ($1+$2+0.0)/'${nFastQReads}'*100
     }' >> ${statFile}

# remove SAM files; keep BAMs only
rm ${starOutDir}/*.sam

# ENCODE param
#minMappedLength=16
#${STAR} --genomeDir ${genomeDir}  --readFilesIn ${FASTQ} --runThreadN $threads --alignIntronMax 1 --outSAMattributes All --outFilterMultimapNmax ${multiMapMax} --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin ${minMappedLength} --outFileNamePrefix ${starOutDir}
