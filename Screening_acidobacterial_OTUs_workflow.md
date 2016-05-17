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
biom convert -i MASTER_OTU_hdf5_filteredfailedalignments_rdp_collapse_even321798.biom -o Classic_OTU_even321798.txt --to-tsv --header-key="taxonomy"
```

#Pick OTUs belonging to the Acidobacteria from classic OTU_even321798.txt file (only includes OTU ID belonging to the acidobacterial group)
```
grep "Acidobacteria" Classic_OTU_even321798.txt > AcidoOTUs.txt
awk '{print$1}' AcidoOTUs.txt >> Acidolist.txt
```
```
Number of acidobacterial sequences = 1,597,270 sequences
Number of OTUs = 2391 otus
```

#screening acidobacterial OTU map from “OTU_map.uc file” using shell script "acidoscreening.sh" file
```
#!/bin/bash
# trim centralia sequences so we can test pipeline and analysis
# use pandaseq to merge reads - requires name list (file <list.txt> in same folder as $
for file in $(<Acidolist.txt)
do
   grep -w ${file} OTU_map1.uc >> ./sequences/Qiime_acido_map.txt
done
```

***
problem = grep OTU informations that are not involved in Acidobacteria
***
```
Example
H	4526	253	99.6	+	0	0	512I253M614I	C11D01.1718	210306
H	6544	253	98.4	+	0	0	516I253M695I	C08D02.1718	1111657
H	5804	253	98.8	+	0	0	508I253M760I	C08D02.29889	1718
H	1718	253	98.0	+	0	0	517I253M729I	C08D02.52093	581409
H	1718	253	97.6	+	0	0	517I253M729I	C08D02.89689	581409
H	1718	253	98.0	+	0	0	517I253M729I	C08D02.101541	581409
H	5804	253	99.6	+	0	0	508I253M760I	C08D02.138214	1718
H	5804	253	100.0	+	0	0	508I253M760I	C08D02.138217	1718
H	8239	253	97.6	+	0	0	476I253M725I	C11D03.1718	4352522
H	5613	253	97.6	+	0	0	439I253M697I	C12D03.1718	61819
H	9118	253	99.6	+	0	0	472I253M627I	C10D02.1718	111933
H	1718	253	97.2	+	0	0	517I253M729I	C10D02.15445	581409
H	1718	253	97.2	+	0	0	517I253M729I	C10D02.73053	581409
H	1312	253	98.4	+	0	0	485I253M691I	C18D05.1718	671400
H	1718	253	97.2	+	0	0	517I253M729I	C18D05.27623	581409
```
***

*New command for grep OTUs
```
awk '$10 == 1718' OTU_map1.uc
```
```
H       5804    253     98.8    +       0       0       508I253M760I    C08D02.29889    1718
H       5804    253     99.6    +       0       0       508I253M760I    C08D02.138214   1718
H       5804    253     100.0   +       0       0       508I253M760I    C08D02.138217   1718
H       5804    253     97.2    +       0       0       508I253M760I    C18D05.159942   1718
H       5804    253     97.2    +       0       0       508I253M760I    C16D03.31522    1718
H       5804    253     97.2    +       0       0       508I253M760I    C16D03.31524    1718
H       5804    253     97.2    +       0       0       508I253M760I    C09D04.84959    1718
H       5804    253     97.2    +       0       0       508I253M760I    C02D03.137417   1718
H       5804    253     97.2    +       0       0       508I253M760I    C16D01.51912    1718
H       5804    253     97.2    +       0       0       508I253M760I    C07D03.97862    1718
H       5804    253     97.2    +       0       0       508I253M760I    C06D01.145156   1718
H       5804    253     97.2    +       0       0       508I253M760I    Mock.14422      1718
H       5804    253     97.2    +       0       0       508I253M760I    Mock.52525      1718
H       5804    253     97.2    +       0       0       508I253M760I    C01D03.151278   1718
H       5804    253     97.2    +       0       0       508I253M760I    C01D03.151281   1718
H       5804    253     97.2    +       0       0       508I253M760I    C01D03.183601   1718
H       5804    253     97.2    +       0       0       508I253M760I    C01D03.183603   1718
H       5804    253     99.2    +       0       0       508I253M760I    C17D01.38150    1718
H       5804    253     99.2    +       0       0       508I253M760I    C17D01.38151    1718
H       5804    253     97.2    +       0       0       508I253M760I    C18D06.122883   1718
H       5804    253     97.2    +       0       0       508I253M760I    C14D03.5889     1718
H       5804    253     97.2    +       0       0       508I253M760I    C08D03.186861   1718
H       5804    253     97.2    +       0       0       508I253M760I    C08D03.186867   1718
H       5804    253     97.2    +       0       0       508I253M760I    C09D06.74692    1718
H       5804    253     97.2    +       0       0       508I253M760I    C10D01.55291    1718
H       5804    253     97.2    +       0       0       508I253M760I    C10D01.190466   1718
```
***
***
Q1. Can we make it automate, too?
***

#Change “OTU_dn_###” to “OTU” in the Qiime_acido_map.txt file
```
sed -i 's/OTU_dn_/OTUdn/g' OTU_map_NC.uc
```


#Convert uclust mapping file to qiime format file.
```
python Convert_UC_2_Qiime2.py sequences/Qiime_acido_map.txt > QIIME_ACIDO_MAP2.txt
```
***
Converted file has only 2066 OTUs with sequence name. Something is problem.
***

#Prepare fasta file from combined_merged.fastq file
```
convert_fastaqual_fastq.py -c fastq_to_fastaqual -f combined_merged.fastq -o fastaqual
grep -w ">" 97_otus.fasta >> combined_merged.fna
cut -d ";" -f 1 fastaqual/combined_merged.fna > fastaqual/CleanedHeaders_combined_green_merged.fna
```
***
But, it only changes Denovo OTUs, not the greengene references's one.
***

#Pick acidobacterial sequences from total sequence set
```
filter_fasta.py -f fastaqual/CleanedHeaders_combined_green_merged.fna -m ../QIIME_ACIDO_MAP2.txt -o Acido_from_RefNo.fasta
```
* Result: "Acido_from_RefNo.fasta" has about 95,000 sequences were picked. It probably only from the Denovo OTUs.
***
Solution?
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
