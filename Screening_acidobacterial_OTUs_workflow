---
layout: page
title: "Screening acidobacterial sequences from Centralia sequence dataset"
comments: true
date: 2016-05-19
---

##This Workflow for Oligotyping analyses for acidobacterial community ecology in Centralia.

#Intro to computing workflow to screening OTUs belonging to the phylum Acidobacteria


##Overarching Goal
* This tutorial will contribute towards an understanding of **Oligotyping analysis for specific bacterial lineage**

##preparing files
* source /mnt/research/ShadeLab/software/loadanaconda2.sh
*	OTU_map.uc
*	MASTER_OTU_hdf5_filteredfailedalignments_rdp_collapse_even321798.biom
*
*

***
Picking acidobacterial sequences from full sequence dataset
***
1. convert OTU_table.biom file to .txt file
biom convert -i MASTER_OTU_hdf5_filteredfailedalignments_rdp_collapse_even321798.biom -o Classic_OTU_even321798.txt --to-tsv --header-key="taxonomy"

2. pick OTU ID belonging to the Acidobacteria from classic OTU.txt file (using excel)

3. make “Acidolist.txt file” (only includes OTU ID belonging to the acidobacterial group). (using excel)

4. screening acidobacterial OTU map from “OTU_map.uc file” using shell script “acidoscreening.sh file”
#!/bin/bash
# trim centralia sequences so we can test pipeline and analysis
# use pandaseq to merge reads - requires name list (file <list.txt> in same folder as $
for file in $(<Acidolist.txt)
do
    grep -w ${file} OTU_map1.uc >> ./sequences/Qiime_acido_map.txt
done

5. convert “OTU_dn_###” to “OTUdn” in the Qiime_acido_map.txt file
# because, to convert mapping file to Qiime format using Convert_UC_2_Qiime.py, but when I try to do it, “OTU_dn_” makes problem to running it.
Find and replace it using text editor.

6. convert uclust file to qiime format file.
python Convert_UC_2_Qiime2.py sequences/Qiime_acido_map.txt > Qiime_map.txt

7. pick acidobacterial tax_assignment information from “MASTER_RepSeqs_filteredfailedalignments_tax_assignments_acido.txt” file, using acidotaxscreening.sh file
#!/bin/bash
# trim centralia sequences so we can test pipeline and analysis
# use pandaseq to merge reads - requires name list (file <list.txt> in same folder as this script) of forward and reverse reads to be merged using the panda-seq program
for file in $(<Acidolist.txt)
do
    grep -w ${file} MASTER_RepSeqs_filteredfailedalignments_tax_assignments_acido.txt >> ./sequences/tax_assignments_acido.txt 
done

8. prepare fasta file from uniques_combined_merged.fastq file
convert_fastaqual_fastq.py -c fastq_to_fastaqual -f uniques_combined_merged.fastq -o fastaqual
cut -d ";" -f 1 uniques_combined_merged.fna > CleanedHeaders_uniques_combined_merged.fna

9. pick acidobacterial sequences from total sequence set
filter_fasta.py -f fastaqual/CleanedHeaders.fna -m ../QIIME_ACIDO_MAP2.txt -o Acido_from_RefNo.fasta


OTU picking based on the RDP database
# Reference-based OTU picking using usearch_global: Cluster sequences at 97% identity to the RDP database, version (aligned_trainset14_032015.rdp.fasta)
usearch -usearch_global nocrap_denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db aligned_trainset14_032015.rdp.fasta -notmatchedfq RefNoMatch_nocrap_denoised_nosigs_uniques_combined_merged.fastq -strand plus -uc RefMatchOTUMap_nocrap_denoised_nosigs_uniques_combined_merged.uc -dbmatched RDP_97_rep_set_matched.fa
