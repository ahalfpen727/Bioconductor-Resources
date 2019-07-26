	
 gunzip -c reference/human/hg19.fa.gz > reference/human/hg19.fa
$ time bowtie2-build reference/human/hg19.fa reference/human/hg19
   mkdir -p /var/data/bi/reference/human
$ rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ /var/data/bi/reference/human
$ cat reference/human/chr*.fa.gz > /var/data/bi/reference/human/hg19.fa.gz
    $ mkdir -p /var/data/bi/reference/human
    $ rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ /var/data/bi/reference/human
    $ cat reference/human/chr*.fa.gz > /var/data/bi/reference/human/hg19.fa.gz
    1
    2
    3
    	
    $ mkdir -p /var/data/bi/reference/human
    $ rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ /var/data/bi/reference/human
    $ cat reference/human/chr*.fa.gz > /var/data/bi/reference/human/hg19.fa.gz
    Use Table Browser to download UCSC gene annotations for hg19 in GTF format.

Alternatively, you can download a pre-build packaging of raw sequences and various annotation information:
Shell
$ mkdir -p /var/data/bi/reference/prebuild
$ wget -bqc -P /var/data/bi/reference/prebuild/ ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh38.tar.gz
$ tar -zxvf /var/data/bi/reference/prebuild/Homo_sapiens_Ensembl_GRCh37.tar.gz
1
2
3
	
$ mkdir -p /var/data/bi/reference/prebuild
$ wget -bqc -P /var/data/bi/reference/prebuild/ ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz
$ tar -zxvf /var/data/bi/reference/prebuild/Homo_sapiens_Ensembl_GRCh37.tar.gz

 gunzip -c reference/human/hg19.fa.gz > reference/human/hg19.fa
$ time bowtie2-build reference/human/hg19.fa reference/human/hg19
 
time bowtie2 -p 8 -x reference/human/hg19 -1 input_1.fastq -2 input_2.fastq > output.sam
	
When runnning TopHat with the -G option, it will build a Bowtie index from the provided reference annotation each time. TopHat will output the following lines when doing so:

[..] Building transcriptome data files tophat_out/tmp/genes
[..] Building Bowtie index from genes.fa

Building the index can take quite some time. You can pre-build this index once and use it afterwards each time you run TopHat. To do so, invoke TopHat just with the -G and the --transcriptome-index parameters. TopHat will build an index in the provided location.

$ tophat -G /var/data/bi/reference/prebuild/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf --transcriptome-index /var/data/bi/reference/prebuild/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/transciptome_index/genes /var/data/bi/reference/prebuild/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome

Note that you still have to provide the <bowtie_index> argument. Read the TopHat manual as well. Output:

Building transcriptome files with TopHat v2.1.0
-----------------------------------------------
Checking for Bowtie
Bowtie version:	 2.2.6.0
Checking for Bowtie index files (genome)..
Checking for reference FASTA file
Building transcriptome data files /var/data/bi/reference/prebuild/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/transciptome_index/genes
Building Bowtie index from genes.fa
-----------------------------------------------
Transcriptome files prepared. This was the only task requested.

Files created:

$ ls -lh /var/data/bi/reference/prebuild/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/transciptome_index/
total 1,7G
-rw-rw-r-- 1 root root 125M Jan  3 17:08 genes.1.bt2
-rw-rw-r-- 1 root root  69M Jan  3 17:08 genes.2.bt2
-rw-rw-r-- 1 root root 1,7M Jan  3 16:56 genes.3.bt2
-rw-rw-r-- 1 root root  69M Jan  3 16:56 genes.4.bt2
-rw-rw-r-- 1 root root 305M Jan  3 16:56 genes.fa
-rw-rw-r-- 1 root root  26M Jan  3 16:56 genes.fa.tlst
-rwxrwxr-x 1 root root 874M Jan  3 16:55 genes.gff
-rw-rw-r-- 1 root root 125M Jan  3 17:19 genes.rev.1.bt2
-rw-rw-r-- 1 root root  69M Jan  3 17:19 genes.rev.2.bt2
-rw-rw-r-- 1 root root   24 Jan  3 16:56 genes.ver


 
