
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

