#!/bin/bash
# File: copy_paste_figures.sh
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
    OUTPUT_F="output_fig_SPAR"  #${SPARDIR}/${EXPNAME}
        mkdir "$PWD/$OUTPUT_F"
  else
    # if output directory is specified, use it
        if [ ! -d "$OUTPUT_F" ]; then
 echo   "\n i am in se $OUTPUT_F"
        mkdir -pv "$OUTPUT_F"
        fi
  fi

for png in $INPUT_F/*/figures/*.png; do
[ -f "$png" ] || continue
  str="${png##$INPUT_F/}"      #remove guiding path
  filename="${str%%.trimmed*}" #remove trailing path
  figur="${str##*/}"           #remove guiding
  figure="${figur%%.png}"      #get figure name
  cp -uv  $png $OUTPUT_F/${figure}_${filename}.png
  done
 
for html in $INPUT_F/*/figures/*.html;do
[ -f "$html" ] || continue
  str="${html##$INPUT_F/}"
  filename="${str%%.trimmed*}"
  figur="${str##*/}"
  figure="${figur%%.html}"
  cp -uv  $html $OUTPUT_F/${figure}_${filename}.html
  done

css_f=$(find $INPUT_F/ -type d  -name "css" | head -1)
js_f=$(find $INPUT_F/ -type d  -name "js" | head -1)

cp -R {$css_f,$js_f} $OUTPUT_F/

## Konstantinos Con_Yel
