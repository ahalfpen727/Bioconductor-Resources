#!bin/bash

# use samtools to index fasta genome

samtools faidx human_hg19.fa

#then you can generate the whole genome bed file by entering, using awk for instance, a "0" column between the first 2 columns of the .fai file generated previously

awk '{print $1 "\t0\t" $2}' human_hg19.fa.fai > human_hg19.bed