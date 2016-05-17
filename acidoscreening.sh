#!/bin/bash
# trim centralia sequences so we can test pipeline and analysis
# use pandaseq to merge reads - requires name list (file <list.txt> in same folder as this script) of forward and reverse reads to be merged using the panda-seq program

for file in $(<Acidolist.txt)
do
    grep -w ${file} OTU_map1.uc >> ./sequences/Qiime_acido_map.txt 

done
