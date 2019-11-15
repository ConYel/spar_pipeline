#!/bin/bash
# File: cp_figures_data.sh
## Konstantinos Con_Yel
## Copy all figures from resulted sample's folders of SPAR workflow to "output_fig_SPAR/"

INPUT_F="$1"
OUTPUT_F="$2"

if [ $# -lt 1 ]; then
  echo "USAGE: <INPUT_FOLDER> folder with all results <OUTPUT_FOLDER>"
  exit 1
fi


  if [ -z "$OUTPUT_F" ]; then
 echo "$OUTPUT_F"
    # if no output directory is provided, use working directory to make a folder
    OUTPUT_F="output_figures_SPAR"  #${SPARDIR}/${EXPNAME}
        mkdir "$PWD/$OUTPUT_F"
  else
    # if output directory is specified, use it
        if [ ! -d "$OUTPUT_F" ]; then
  echo   "\n i will copy them in  $OUTPUT_F"
        mkdir -pv "$OUTPUT_F"
        fi
  fi

for dir in $INPUT_F/*/figures; do  
[ -d "$dir" ] || continue
	str="${dir##${INPUT_F}/}" #remove guiding path
	filename="${str%%.trimmed*}" #remove trailing path
	cp -R -uv ${dir} ${OUTPUT_F}/
	mv ${OUTPUT_F}/figures ${OUTPUT_F}/${filename}
 done

for file in $INPUT_F/*/results/smRNA_gene_expression.xls; do  
[ -f "$file" ] || continue
  str="${file##$INPUT_F/}"
  filename="${str%%.trimmed*}"
  table="${str##*/}"
  cp -uv  $file $OUTPUT_F/${filename}/${filename}_${table}
 done
## Konstantinos Con_Yel 
