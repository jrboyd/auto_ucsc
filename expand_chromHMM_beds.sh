#!/bin/bash
for f in *ments.bed; do
echo $f
#f=MCF10Abza_30_999_segments.bed
new=${f/".bed"/".expanded.bed"}
cat $f | awk '
  BEGIN {OFS="\t"; p_chr=""; p_s=-1; p_e=-1; p_state=""} 
    {chr=$1; s=$2; e=$3; state=$4;
    while (e - s >= 200){ 
      print chr,s,s+200,state; s = s + 200
    }
  }
' > $new
done
#    {if (e - s == 200) print $0}
#'
