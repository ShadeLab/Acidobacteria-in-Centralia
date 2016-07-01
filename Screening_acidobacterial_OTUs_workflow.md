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
* Load QIIME
```
source /mnt/research/ShadeLab/software/loadanaconda2.sh
```
* OTU_map.uc
* MASTER_OTU_hdf5_filteredfailedalignments_rdp.biom

***
Picking acidobacterial sequences from full sequence dataset
***
*convert OTU_table.biom file to .txt file
```
biom convert -i MASTER_OTU_hdf5_filteredfailedalignments_rdp.biom -o Classic_OTU_rdp.txt --to-tsv --header-key="taxonomy"
```

#Pick OTUs belonging to the Acidobacteria from Classic_OTU_collapse.txt file (only includes OTU ID belonging to the acidobacterial group)
```
grep "Acidobacteria" Classic_OTU_rdp.txt > AcidoOTUs.txt
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

for file in $(<Acidolist.txt)
do
    grep -w ^${file} OTU_map_qiime.txt >> ./sequences/${file}_map.txt 

done
```
or
```
#!/bin/bash

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

#grep sequences belonging to each OTUs individually (using qsub script) 
```
#! /bin/bash --login
# Time job will take to execute (HH:MM:SS format)
#PBS -l walltime=010:00:00
# Memory needed by the job
#PBS -l mem=264Gb
# Number of shared memory nodes required and the number of processors per node
#PBS -l nodes=1:ppn=8
# Make output and error files the same file
#PBS -j oe
# Send an email when a job is aborted, begins or ends
#PBS -m abe
# Give the job a name
#PBS -N batch_usearch_centralia_03dec2015
# _______________________________________________________________________#

cd ${PBS_O_WORKDIR}

### i.  make fasttree
     source /mnt/research/ShadeLab/software/loadanaconda2.sh
     python Convert_UC_2_QiimeJS.py screening_otus/AcidoOTUs/${file}.txt > screening_otus/AcidoOTUs/${file}_qiime_map.txt

     filter_fasta.py -f OTUfromRDP/fastaqual/combined_merged.fna -m screening_otus/AcidoOTUs/${file}_qiime_map.txt -o screening_otus/Acidosequences/${file}.fasta

 _______________________________________________________________________#
# PBS stats
cat ${PBS_NODEFILE}
env | grep PBS
qstat -f ${PBS_JOBID}
```

***
#Next step
***
* Following the Meren Lab's protocol? (link below)
```
http://merenlab.org/projects/oligotyping/
```
***
Commands
***
```
o-trim-uninformative-columns-from-alignment Acido_seqs_from_combined_merg_aligned.fasta

o-smart-trim Acido_seqs_from_combined_merg_aligned.fasta-TRIMMED -E -o Acido_seqs_from_combined_merg_aligned_sub_trimmed_from_begining.fasta

o-smart-trim Acido_seqs_from_combined_merg_aligned_sub_trimmed_from_begining.fasta -S -o Acido_align_fixed.fasta
```

***
Test with 1000 subsampled aligned sequences
***
```
seqtk sample -s 100 Acido_seqs_from_combined_merg_aligned.fasta 1000 > Acido_seqs_from_combined_merg_aligned_sub.fasta

o-trim-uninformative-columns-from-alignment Acido_seqs_from_combined_merg_aligned_sub.fasta

o-smart-trim Acido_seqs_from_combined_merg_aligned_sub.fasta-TRIMMED -E -o Acido_seqs_from_combined_merg_aligned_sub_trimmed_from_begining.fasta

o-smart-trim Acido_seqs_from_combined_merg_aligned_sub_trimmed_from_begining.fasta -S -o Acido_align_fixed.fasta

entropy-analysis Acido_align_fixed.fasta

oligotype Acido_align_fixed.fasta Acido_align_fixed.fasta-ENTROPY -c 2 -M 10 --quick
```


# test 2
***
Preparing "RDP reference sequences" fasta file. (Informations from Mothur homepage)
***
Down load reference file from Mothur website (http://www.mothur.org/wiki/RDP_reference_files)
```
wget -N http://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/RDPClassifier_16S_trainsetNo14_rawtrainingdata.zip
unzip -o RDPClassifier_16S_trainsetNo14_rawtrainingdata.zip
mv RDPClassifier_16S_trainsetNo14_rawtrainingdata/* ./
```
RDP duplicated the sequence of GL982576_U010303693 in these files and so we need to remove the second copy. (with Mothur)
```
sed 's/[\|]/_/' trainset14_032015.fasta > trainset14_032015.temp.fasta

echo "GL982576_U010303693" > duplicate.accnos

* using Mothur from here.
mothur "#get.seqs(accnos=duplicate.accnos, fasta=trainset14_032015.temp.fasta); unique.seqs()"

mothur "#remove.seqs(accnos=duplicate.accnos, fasta=trainset14_032015.temp.fasta)"

cat trainset14_032015.temp.pick.fasta trainset14_032015.temp.pick.unique.fasta > trainset14_032015.rdp.fasta
```
Taxonomy file and the fasta file that will be our reference.
```
grep ">" trainset14_032015.rdp.fasta | cut -c 2- > trainset14_032015_rmdup.tax
```
Get Taxonomy file properly formatted (Using R code)
```
R code

tax_file <- scan(file="trainset14_032015_rmdup.tax", what="", sep="\n", quiet=TRUE)

accession <- gsub("^(\\S*).*", "\\1", tax_file) #some are separated by tabs or spaces or both

taxonomy <- gsub(".*(Root.*)", "\\1", tax_file)
taxonomy <- gsub(" ", "_", taxonomy)    #remove spaces and replace with '_'
taxonomy <- gsub("\t", "", taxonomy)    #remove extra tab characters
taxonomy <- gsub("[^;]*_incertae_sedis$", "", taxonomy)
taxonomy <- gsub('\"', '', taxonomy) #remove quote marks

levels <- read.table(file="trainset14_db_taxid.txt", sep="*", stringsAsFactors=FALSE)
subs <- levels[grep("sub", levels$V5),]
sub.names <- subs$V2

tax.split <- strsplit(taxonomy, split=";")

remove.subs <- function(tax.vector){
    return(tax.vector[which(!tax.vector %in% sub.names)])
}

no.subs <- lapply(tax.split, remove.subs)
no.subs.str <- unlist(lapply(no.subs, paste, collapse=";"))
no.subs.str <- gsub("^Root;(.*)$", "\\1;", no.subs.str)

write.table(cbind(as.character(accession), no.subs.str), "trainset14_032015.rdp.tax", row.names=F, col.names=F, quote=F, sep="\t")
```
Add sequences involved in the mitochondria or sequences from eukaryotes into sequence dataset
```
wget -N http://mothur.org/w/images/2/24/Trainset10_082014.pds.tgz
tar xvzf Trainset10_082014.pds.tgz
mv trainset10_082014.pds/trainset10_082014* ./
rm -rf trainset10_082014.pds Trainset10_082014.pds.tgz
```
Pull out the extra sequences that are in the pds files
```
mothur "#get.lineage(fasta=trainset10_082014.pds.fasta, taxonomy=trainset10_082014.pds.tax, taxon=Eukaryota-Mitochondria)"
```
This last command gets us the extra “pds” sequences that we can now use to paste on to the end of the normal RDP training set
```
cat trainset14_032015.rdp.tax trainset10_082014.pds.pick.tax > trainset14_032015.pds.tax
cat trainset14_032015.rdp.fasta trainset10_082014.pds.pick.fasta > trainset14_032015.pds.fasta
```
Confirm sequence dataset
```
wc -l *.pds.tax

*results
##    10773 trainset10_082014.pds.tax
##    10801 trainset14_032015.pds.tax
##    21574 total
```
Alignment of this reference file using RDP pipeline
```
https://pyro.cme.msu.edu/index.jsp
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
