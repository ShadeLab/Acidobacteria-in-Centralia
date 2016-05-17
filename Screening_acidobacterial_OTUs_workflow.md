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

#Preparing files
* Load QIIME : source /mnt/research/ShadeLab/software/loadanaconda2.sh
* OTU_map.uc
* MASTER_OTU_hdf5_filteredfailedalignments_rdp_collapse_even321798.biom



***
Picking acidobacterial sequences from full sequence dataset
***
1. convert OTU_table.biom file to .txt file
'''
biom convert -i MASTER_OTU_hdf5_filteredfailedalignments_rdp_collapse_even321798.biom -o Classic_OTU_even321798.txt --to-tsv --header-key="taxonomy"
'''

2. take OTUs belonging to the Acidobacteria from classic OTU.txt file (using excel, based on the taxonomic classification)
   Number of acidobacterial sequences = 1,597,270 sequences

3. make "Acidolist.txt" file (only includes OTU ID belonging to the acidobacterial group). (using excel)

4. screening acidobacterial OTU map from “OTU_map.uc file” using shell script "acidoscreening.sh" file
#!/bin/bash
# trim centralia sequences so we can test pipeline and analysis
# use pandaseq to merge reads - requires name list (file <list.txt> in same folder as $
for file in $(<Acidolist.txt)
do
    grep -w ${file} OTU_map1.uc >> ./sequences/Qiime_acido_map.txt
done

## problem = grep OTU informations that are not involved in Acidobacteria
***
---
Example
1 : H	4526	253	99.6	+	0	0	512I253M614I	C11D01.1718	210306
2 : H	6544	253	98.4	+	0	0	516I253M695I	C08D02.1718	1111657
:H	5804	253	98.8	+	0	0	508I253M760I	C08D02.29889	1718
:H	1718	253	98.0	+	0	0	517I253M729I	C08D02.52093	581409
:H	1718	253	97.6	+	0	0	517I253M729I	C08D02.89689	581409
:H	1718	253	98.0	+	0	0	517I253M729I	C08D02.101541	581409
:H	5804	253	99.6	+	0	0	508I253M760I	C08D02.138214	1718
:H	5804	253	100.0	+	0	0	508I253M760I	C08D02.138217	1718
:H	8239	253	97.6	+	0	0	476I253M725I	C11D03.1718	4352522
:H	5613	253	97.6	+	0	0	439I253M697I	C12D03.1718	61819
:H	9118	253	99.6	+	0	0	472I253M627I	C10D02.1718	111933
:H	1718	253	97.2	+	0	0	517I253M729I	C10D02.15445	581409
:H	1718	253	97.2	+	0	0	517I253M729I	C10D02.73053	581409
:H	1312	253	98.4	+	0	0	485I253M691I	C18D05.1718	671400
:H	1718	253	97.2	+	0	0	517I253M729I	C18D05.27623	581409
---
***


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
