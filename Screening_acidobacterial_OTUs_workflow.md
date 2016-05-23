---
layout: page
title: "Screening acidobacterial sequences from Centralia sequence dataset"
comments: true
date: 2016-05-19
---

#This Workflow for Oligotyping analyses for acidobacterial community ecology in Centralia.
* Intro to computing workflow to screening OTUs belonging to the phylum Acidobacteria


#Overarching Goal
* This tutorial will contribute towards an understanding of **Oligotyping analysis for specific bacterial lineage**

#Preparing files
* Load QIIME : source /mnt/research/ShadeLab/software/loadanaconda2.sh
* OTU_map.uc
* MASTER_OTU_hdf5_filteredfailedalignments_rdp_collapse_even321798.biom

***
Picking acidobacterial sequences from full sequence dataset
***
*convert OTU_table.biom file to .txt file
```
biom convert -i MASTER_OTU_hdf5_filteredfailedalignments_rdp.biom -o Classic_OTU_collapse.txt --to-tsv --header-key="taxonomy"
```

#Pick OTUs belonging to the Acidobacteria from Classic_OTU_collapse.txt file (only includes OTU ID belonging to the acidobacterial group)
```
grep "Acidobacteria" Classic_OTU_collapse.txt > AcidoOTUs.txt
awk '{print$1}' AcidoOTUs.txt >> Acidolist.txt
```
```
Number of acidobacterial sequences = 2,221,021 sequences
Number of OTUs = 2,409 otus
3 OTUs (13 sequences) only observed in the Mock community.
```

#Change uclust OTU map to Qiime map.
```
python Convert_UC_2_QiimeJS.py OTU_map.uc > OTU_map_qiime.txt
```

#Grep sequences belong to the acidobacterial OTUs from mapping file.
```
#!/bin/bash
# trim centralia sequences so we can test pipeline and analysis
# use pandaseq to merge reads - requires name list (file <list.txt> in same folder as this script) of forward and reverse reads to be merged using the panda-seq program

for file in $(<Acidolist.txt)
do
    grep -w ^${file} OTU_map_qiime.txt >> ./sequences/${file}_map.txt 

done
```
or
```
#!/bin/bash
# trim centralia sequences so we can test pipeline and analysis
# use pandaseq to merge reads - requires name list (file <list.txt> in same folder as this script) of forward and reverse reads to be merged using the panda-seq program

for file in $(<Acidolist.txt)
do
    grep -w ^${file} OTU_map_qiime.txt >> ./sequences/total_qiime_acido_map.txt 

done
```

#Pick acidobacterial sequences from total sequence set
```
filter_fasta.py -f combined_merged.fna -m sequences_JSmethod/total_qiime_acido_map.txt -o Acido_seqs_from_combined_merg.fasta
```
***
from OTU table = 2,217,760 sequences
In the Acido_seqs_from_merg_Refs.fasta = 2,172,325 sequences
Acido_seqs_from_combined_merg.fasta = 2,172,325 sequences
balance = 45,435 sequences???
***


***
#Next step
***
* Following the Meren Lab's protocol? (link below)
```
http://merenlab.org/projects/oligotyping/
```

***
Other concept
***
#OTU picking based on the RDP database
```
Reference-based OTU picking using usearch_global: Cluster sequences at 97% identity to the RDP database, version (aligned_trainset14_032015.rdp.fasta)
```
```
usearch -usearch_global nocrap_denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db aligned_trainset14_032015.rdp.fasta -notmatchedfq RefNoMatch_nocrap_denoised_nosigs_uniques_combined_merged.fastq -strand plus -uc RefMatchOTUMap_nocrap_denoised_nosigs_uniques_combined_merged.uc -dbmatched RDP_97_rep_set_matched.fa
```


# tips
*if you want to change name or specific charactor, use "sed" command.
```
sed -i 's/string to be replaced/string to replace it with/g' file to edit
```
```
(example)
sed -i 's/OTU_dn_/OTUdn/g' OTU_map.uc
```

#Prepare fasta file from combined_merged.fastq file
```
convert_fastaqual_fastq.py -c fastq_to_fastaqual -f combined_merged.fastq -o fastaqual
grep -w ">" 97_otus.fasta >> combined_merged.fna
cut -d ";" -f 1 fastaqual/combined_merged.fna > fastaqual/CleanedHeaders_combined_green_merged.fna
```
