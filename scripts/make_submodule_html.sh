submoduleDir=$1
WDIR=$2
ls ${submoduleDir}/*M*.r | \
  awk 'BEGIN{imgDIR="'${WDIR}'";
             #printf "<html>\n<body>\n";
             print "<div id=\"plot-names\">";
            }
     {
       plotname=$0;
       #gsub(/^.+ggplot2_/,"",plotname);
       gsub(/^.+\//,"",plotname); # remove directory
       gsub(/^M[^_]+_/,"",plotname); # remove module number
       gsub(/.r$/,"",plotname); # remove extension
       plotCategory = plotname; gsub(/_.+$/,"",plotCategory);
       gsub(/^[^_]+_/,"",plotname); # remove category name from plot name
       if (plotCategory != plotCategory_prev)
         print "<h3>"plotCategory"</h3>";

       #imgURL = ("http://tesla.pcbi.upenn.edu/~pkuksa/SPAR/" imgDIR "/" plotname ".png")
       imgURL = ("https://www.lisanwanglab.org/SPAR/" imgDIR "/" plotname ".png")

       imgTitle = imgURL;
       gsub(/^.+\//,"",imgTitle);
       gsub(/\.[^.]+$/,"",imgTitle);
       gsub(/_/," ",imgTitle);
       imgTitle=(toupper(substr(imgTitle,1,1)) tolower(substr(imgTitle,2)))
       
       img=("<a href=\"" imgURL "\" target=\"_blank\">" "<figure> <img src=\"" imgURL "\" alt=\"\"><figcaption>" imgTitle "</figcaption></figure></a>");
       plotNo++;

       plot_div[plotNo]=img; 

       plotid=("plot-" plotNo);


       print "<div class=\"plot-link\">";
       if (plotNo>1)
         print ("<a href=\"#" plotid "\" class=\"swap\">" imgTitle "</a>");
       else
         print ("<a href=\"#" plotid "\" class=\"swap active\">" imgTitle "</a>");
       print "</div>";

       plotCategory_prev = plotCategory;
     }END{#printf "</body>\n</html>"
          print  "</div>"; ## plot-names div end
        
         print "<div class=\"plot-result\">";

         for ( i=1; i<=plotNo; ++i)
         {
           plotid=("plot-" i)
           print ("<div class=\"plot-image swap\" id=\"" plotid  "\">");
           print plot_div[i];
           print "</div>";
         }
         print "</div>"; ## plot-result div end

         }' 
#> "${OUTDIR}/${module}_${step}.html"

