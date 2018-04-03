INTABLE=$1 # peak table across experiments coordinate-sorted

maxGap=${2:-3} # maximum allowed gap between peak end points

skipHeader=${4:-1} # set to 0 if no header in INTABLE

hasAnnot=${3:-1}
if [ "${hasAnnot}" -eq 1 ]; then
  nSegmFields=20
else
  nSegmFields=20
fi

if [ $# -lt 1 ]; then
  echo "USAGE: $0 <coord_sorted_peak_table> <maxGap:3> <hasAnnot:0|1> <hasHeader:0|1>"
  exit 1
fi



# split table by strand
# and make sure table is sorted by start coordinate
INTABLEpos="${INTABLE}.$$.pos"
INTABLEneg="${INTABLE}.$$.neg"
cmd="cat ${INTABLE}" 
if [ ${skipHeader} -eq 1 ]; then
  cmd="tail -n +2 ${INTABLE}" 
fi

# 1. Split by strand
${cmd} | sort -k1,1 -k2,2n -k3,3n | \
  awk 'BEGIN{ posTable="'${INTABLEpos}'";
              negTable="'${INTABLEneg}'";
              skipHeader='${skipHeader}'+0; }
     {
        #if (NR==1 && skipHeader) next; # skip header
        if ($6=="+")
           print > posTable
        else 
           print > negTable
     }'

# 2. For each strand, first cluster 5', then 3' endpoints
for f in "${INTABLEpos}" "${INTABLEneg}"; do
  nFields=$(awk '{print NF; exit}' ${f})
  newStartFieldIdx=$((${nFields}+1))
  newEndFieldIdx=$((${nFields}+2))
  clustLfile=${f}.LeftClusters
  clustRfile=${f}.RightClusters
  awk 'BEGIN{OFS="\t"}
       {
        newStart=$2+$15-1; # adjust start to max 5p pos (0-based)
        newEnd=$2+$16-1 + 1; # ajust end to max 3p pos (1-based)
        if (newStart >= newEnd)
        {
           # skip "bad" locus
           print $0, newStart, newEnd > "'${f}'.bad_loci"
           next;
        }
        print $0, newStart, newEnd
       }' ${f} | sort -k1,1 -k${newStartFieldIdx},${newStartFieldIdx}n | \
  awk 'function process_cluster(clust, clust_pos, clust_expr, n)
       {

          split(clust[1],a,"\t")
          leftmost=a[2] # cluster left-most coordinate 
          split(clust[n],a,"\t")
          rightmost=a[2] # cluster right-most coordinate
          # aggregate expression per positions in the cluster
          split("",expr,":");
          for (i=1; i<=n; ++i)
             expr[clust_pos[i]]+=clust_expr[i];

          # find maximum position
          # max_pos = position corresponding to the max no reads
          max_expr = 0;
          max_pos = rightmost;
          for (p in expr)
          {
              #print p, expr[p]
              if (expr[p] > max_expr)
              {
                 max_expr = expr[p]
                 max_pos = p;
              }
          }
          #printf "Lcluster%d\t%d
       }
     BEGIN{FS="\t"; OFS="\t"; clusterNo=1; maxGap='${maxGap}'+0; startIdx='${newStartFieldIdx}'+0; endIdx='${newEndFieldIdx}'+0; clustLfile="'${clustLfile}'"}
     {
       chr = $1;
       max5pPos = $15;
       chrStart = $startIdx #+ max5pPos - 1;
       chrEnd = $endIdx;
       strand = $6;

       if ( chr_prev && (((chrStart-chrStart_prev) > maxGap) || (chr!=chr_prev)) )
       {
           n = i; # number points in the cluster 
           leftmost=clust_pos[1]
           rightmost=clust_pos[n]
           # aggregate expression per positions in the cluster
           split("",expr);
           for (i=1; i<=n; ++i)
             expr[clust_pos[i]]+=clust_expr[i];
           max_expr = 0;
           max_pos = rightmost;
           for (p in expr)
           {
              if (expr[p] > max_expr)
              {
                 max_expr = expr[p]
                 max_pos = p;
              }
           }
           if (leftmost > rightmost || max_pos > rightmost || max_pos < leftmost)
           {
             print "ERROR: clusterNo="clusterNo;
             print n, leftmost, rightmost, max_pos;
             for (i=1; i<=n; ++i)
               print clust[i];
             #exit
           }
 
           for (i=1; i<=n; ++i)
             #print clust[i], "LC_"clusterNo, "left_"leftmost, "right_"chrStart_prev, "n_"n;
             print clust[i], "LC_"clusterNo, "left_"leftmost, "right_"rightmost, "max_pos_"max_pos, "max_expr_"max_expr, "n_"n;

           print "LC_"clusterNo, leftmost, rightmost, max_pos, max_expr, n > clustLfile

           # reset cluster
           leftmost = chrStart;
           split("",clust);
           split("",clust_pos)
           split("",clust_expr)
           i = 0;
           ++clusterNo;
       }
       ++i; # number of points in the cluster
       clust[i] = $0; # store points in the cluster
       clust_pos[i] = chrStart
       max5pPct = $13
       #clust_expr[i] = $5 * max5pPct; # raw read count
       clust_expr[i] = $27 * max5pPct; # RPM
       strand_prev = strand;
       chr_prev = chr;
       chrStart_prev = chrStart;
       chrEnd_prev = chrEnd;
     }
     END{
           # report last cluster
           n = i; # number points in the cluster 
           leftmost=clust_pos[1]
           rightmost=clust_pos[n]
           # aggregate expression per positions in the cluster
           split("",expr);
           for (i=1; i<=n; ++i)
             expr[clust_pos[i]]+=clust_expr[i];
           max_expr = 0;
           max_pos = rightmost;
           for (p in expr)
           {
              if (expr[p] > max_expr)
              {
                 max_expr = expr[p]
                 max_pos = p;
              }
           }
           if (leftmost > rightmost || max_pos > rightmost || max_pos < leftmost)
           {
             print "ERROR: clusterNo="clusterNo;
             print n, leftmost, rightmost, max_pos;
             for (i=1; i<=n; ++i)
               print clust[i];
             #exit
           }
 
           for (i=1; i<=n; ++i)
             #print clust[i], "LC_"clusterNo, "left_"leftmost, "right_"chrStart_prev, "n_"n;
             print clust[i], "LC_"clusterNo, "left_"leftmost, "right_"rightmost, "max_pos_"max_pos, "max_expr_"max_expr, "n_"n;

           print "LC_"clusterNo, leftmost, rightmost, max_pos, max_expr, n > clustLfile
     }' | sort -k1,1 -k${newEndFieldIdx},${newEndFieldIdx}n | \
  awk 'BEGIN{FS="\t"; OFS="\t"; clusterNo=1; maxGap='${maxGap}'+0;startIdx='${newEndFieldIdx}'+0; clustRfile="'${clustRfile}'"}
     {
       # input is now re-sorted by the end-point
       chr = $1;
       chrStart = $startIdx; # end
       chrEnd = $2;
       strand = $6;

       if ( chr_prev && (((chrStart-chrStart_prev) > maxGap) || (chr!=chr_prev)) )
       {
           # report cluster if maxGap is exceeded or chr changed
           n = i; 
           leftmost=clust_pos[1]
           rightmost=clust_pos[n]
           if (leftmost > rightmost)
           {
             print "ERROR"
             print clust[1]
             print clust[n]
             exit
           }
           # aggregate expression per positions in the cluster
           split("",expr,":");
           for (i=1; i<=n; ++i)
             expr[clust_pos[i]]+=clust_expr[i];
           max_expr = 0;
           max_pos = rightmost;
           for (p in expr)
           {
              #print p, expr[p]
              if (expr[p] > max_expr)
              {
                 max_expr = expr[p]
                 max_pos = p;
              }
           }
           for (i=1; i<=n; ++i)
             #print clust[i], "RC_"clusterNo, "left_"leftmost, "right_"chrStart_prev, "n_"n;
             print clust[i], "RC_"clusterNo, "left_"leftmost, "right_"rightmost, "max_pos_"max_pos, "max_expr_"max_expr, "n_"n;
           print "RC_"clusterNo, leftmost, rightmost, max_pos, max_expr, n > clustRfile
           # reset cluster
           leftmost = chrStart;
           split("",clust,":");
           split("",clust_pos,":")
           split("",clust_expr,":")
           i = 0;
           ++clusterNo;
       }
       ++i; # number of points in the cluster
       clust[i] = $0; # points in the cluster
       max3pPos = $16;
       clust_pos[i] = chrStart #+ max3pPos - 1;
       #clust_expr[i] = $5; # raw read count
       same3pPct = $14;
       clust_expr[i] = $27 * same3pPct; # RPM
       strand_prev = strand;
       chr_prev = chr;
       chrStart_prev = chrStart;
       chrEnd_prev = chrEnd;
     }
     END{
           # report last cluster
           n = i; # number points in the cluster 
           leftmost=clust_pos[1]
           rightmost=clust_pos[n]
           # aggregate expression per positions in the cluster
           split("",expr);
           for (i=1; i<=n; ++i)
             expr[clust_pos[i]]+=clust_expr[i];
           max_expr = 0;
           max_pos = rightmost;
           for (p in expr)
           {
              if (expr[p] > max_expr)
              {
                 max_expr = expr[p]
                 max_pos = p;
              }
           }
           if (leftmost > rightmost || max_pos > rightmost || max_pos < leftmost)
           {
             print "ERROR: clusterNo="clusterNo;
             print n, leftmost, rightmost, max_pos;
             for (i=1; i<=n; ++i)
               print clust[i];
             #exit
           }
 
           for (i=1; i<=n; ++i)
             print clust[i], "RC_"clusterNo, "left_"leftmost, "right_"rightmost, "max_pos_"max_pos, "max_expr_"max_expr, "n_"n;

           print "RC_"clusterNo, leftmost, rightmost, max_pos, max_expr, n > clustRfile
     }' | sort -k1,1 -k2,2n -k3,3n > ${f}.withclust #${INTABLEpos} 

     awk 'BEGIN{startIdx='${newStartFieldIdx}'+5; endIdx=startIdx+6;OFS="\t"}{print $1, $startIdx, $endIdx, NR, 1, $6}' ${f}.withclust | sed -e 's/max_pos_//g' | sort -k1,1 -k2,2n -k3,3n -u | awk 'BEGIN{OFS="\t"}{$4=("C" NR); print}'> ${f}.consensus_peaks

done

cat ${INTABLEpos}.withclust ${INTABLEneg}.withclust | \
   awk 'BEGIN{OFS="\t"; consPeakFileBothStrands="'${INTABLE}'.consensus_peaks"; hasAnnot='${hasAnnot}'+0; nSegmFields='${nSegmFields}'}
        { l=$(NF-8);
          r=$(NF-2);
          gsub(/^.+_/,"",l); gsub(/^.+_/,"",r);
          chr=$1; strand=$6;
          nl = $(NF-6); nr = $NF;
          h=(chr "\t" l "\t" r "\t" strand);
          rnaClass="OTHER";
          rnaID="NA";
          if (hasAnnot==1)
          {
            rnaClass=$(nSegmFields + 5); 
            rnaID=$(nSegmFields + 4); 
          }  
          if (!clust[h])
          {
            ccnt+=1; clust[h]=("C" ccnt);
            #if (nl!=nr) {print h, nl, nr > "/dev/stderr"; exit};
            gsub(/n_/,"",nl);
            gsub(/n_/,"",nr);
            #clustScore[h]=nl;
            clustScore[h]=(nl ":" nr);
            clustAnnot[h]=(rnaClass "\t" rnaID);
          };
          print $0, clust[h]
        }
        END{ # print consensus peaks 
             for (c in clust)
             {
                split(c,a,"\t");
                cChr = a[1]; cStart=a[2]; cEnd=a[3]; cStrand=a[4];
                cID = clust[c];
                cScore = clustScore[c];
                print cChr, cStart, cEnd, cID, cScore, cStrand, clustAnnot[c] > consPeakFileBothStrands
             }
        }' | sort -k1,1 -k2,2n -k3,3n -k6,6 > ${INTABLE}.with_cluster_annot

        sort -k1,1 -k2,2n -k3,3n -k6,6 ${INTABLE}.consensus_peaks -o ${INTABLE}.consensus_peaks

echo "Consensus peaks: ${INTABLE}.consensus_peaks"
echo "Table with consensus peak annotations: ${INTABLE}.with_cluster_annot"
#cat ${INTABLEpos}.consensus_peaks ${INTABLEneg}.consensus_peaks | sort -k1,1 -k2,2n -k3,3n -k6,6 | awk 'BEGIN{OFS="\t"}{$4=("C" NR); print}' > ${INTABLE}.consensus_peaks

#rm ${INTABLEpos}*
#rm ${INTABLEneg}*
