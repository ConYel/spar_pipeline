#awk 'NR%4==2' SPAR-master/testdata/adipose1_SRR772426_trimmed.fastq | sort -S2G | uniq -c | awk '{id=("read_" NR "_x" $1); read=$2; print ">"id; print read;}' > tmp.collapsed

# collapse FASTQ file to FASTA with unique reads and their counts
INFASTQ=$1
# this seems to be inefficient
<<COMMENT
awk 'NR%4==2' ${INFASTQ} | sort -S${MAXMEM} | uniq -c | \
  awk '{
         id=("read_" NR "_x" $1);
         read=$2;
         print ">"id"\n"read;
       }'
COMMENT

# this is faster, but can use a lot of memory for hashing
awk '{
       if (NR%4==2) ++read[$1];
     }
     END{  for (r in read)
           { 
              ++readNo;
              id=("read_" readNo "_x" read[r]);
              print ">"id"\n"r
           }
        }' ${INFASTQ} 

#time awk '{if (NR%4==2) ++read[$1]}END{for (r in read) { ++cnt; id=("read_" cnt "_x" read[r]); print ">"id"\n"r}}' /data/users/pkuksa/datasets/DASHR/DASHR_1.0_trimmed_fastq/merged/adipose_merged_trimmed.fastq  > /data/users/pkuksa/datasets/tmp.collapsed2
