HTMLPEAKTABLE=$1
BASEURL=$2
awk 'BEGIN{FS="<td>"; OFS="<td>"; baseUrl="'${BASEURL}'"}
     {
       if ($0~/^<tr>/)
       {
         if ($5~/<a href/)
         {
             split($5,a,"\"");
             a[2]=(a[2] "&hgt.customText=" baseUrl);
             #$5=(a[1] a[2] a[3]);
             $5=(a[1] "\"" a[2] "\"" a[3])
 
             # add hidden div to preview / pre-render UCSC Genome Browser
         
         #    ucscUrl = a[2];
         #    gsub(/hgTracks/,"hgRenderTracks",ucscUrl);
         #    divHtml=("<div><img src=\"" ucscUrl "\"/></div>");
         #    gsub(/<\/a>/, (divHtml "</a>"), $5);
         } 
       }
       print;
     
     }' ${HTMLPEAKTABLE}
