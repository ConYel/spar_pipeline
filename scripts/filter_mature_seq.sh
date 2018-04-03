INMAT=$1 # output from ovelap_mat_seq.sh
exprThres=${2:-0.8}

${GAWK} 'BEGIN{OFS="\t"; peakIDidx=11; peakExprIdx=peakIDidx+1;
            expr_thres = '${exprThres}'+0.0; cntPeakMatProd=0;}
    {
      if ($0~/^#/) { print; next; } # print header
      peakID = $peakIDidx;
      strand = $6;
      if (peakID == peakID_prev && strand == strand_prev)
      {
        cntPeakMatProd++;
        peakMatProd[$5]=cntPeakMatProd; # store expression of the mat product
        matProd[cntPeakMatProd]=$0;
      }
      else
      {
        # process peak
        # output mature product with cumulative expression at p% of the peak
        if (peakID_prev)
        {
          n = asorti(peakMatProd, sortedExpr,"@ind_num_desc");
          matTotalExpr = 0;
          nOut = 0;
          #for (i in sortedExpr)
          #   print i, sortedExpr[i]
          #if (sortedExpr[1]/peakExpr <= 0.2) next;
          for (i = 1; i <= n; ++i)
          { 
             nOut+=1;
             matTotalExpr+=sortedExpr[i];
             #printf "matTotalExpr=%.4f\n", matTotalExpr;
             print matProd[peakMatProd[sortedExpr[i]]];
             if (matTotalExpr/peakExpr>=expr_thres)
               break;
          }
          #print peakID_prev,n,nOut
         }
         # reset mature product array
         cntPeakMatProd = 0;
         split("",peakMatProd,"");
         split("",matProd,"");
         
         peakExpr = $peakExprIdx;
         cntPeakMatProd++;
         peakMatProd[$5]=cntPeakMatProd; # store expression of the mat product
         matProd[cntPeakMatProd]=$0;
      }
      peakID_prev = peakID;
      strand_prev = strand; 
    
    }END{ # output last peak 
      if (peakID_prev)
      {
          n = asorti(peakMatProd, sortedExpr,"@ind_num_desc");
          matTotalExpr = 0;
          nOut = 0;
          #for (i in sortedExpr)
          #   print i, sortedExpr[i]
          for (i = 1; i <= n; ++i)
          { 
             nOut+=1;
             matTotalExpr+=sortedExpr[i];
             #printf "matTotalExpr=%.4f\n", matTotalExpr;
             print matProd[peakMatProd[sortedExpr[i]]];
             if (matTotalExpr/peakExpr>=expr_thres)
               break;
          }
          #print peakID_prev,n,nOut
      }
    }' ${INMAT}  | sort -k1,1 -k2,2n -k3,3n -k6,6
