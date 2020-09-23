#################################################################################################
#                                                                                               #
#                               Fix gtfs and build RSEM reference                               #       
#                                                                                               #
#   Download from ucsc tables the BED file for ucsc genes and the corresponding kgXref file     #
#                                                                                               #
#################################################################################################
source ~xz18w/ve
~xz18w/bin/group_bed_v2.py -i hg19_ucsc.bed -o hg19_ucsc_clusters_new.bed
grep "_TR.1" hg19_ucsc_clusters_new.bed | awk -F "	chr" '{print "chr"$2}' > hg19_ucsc_clusters_new_tr1.bed
awk -F "\t" '{print $1,$2,$3,$4,int($5),$6,$7,$8,$9,$10,$11,$12}' hg19_ucsc_clusters_new_tr1.bed | sed "s/ /\t/g" > hg19_ucsc_clusters_new_tr1_mod.bed
bedToGenePred hg19_ucsc_clusters_new_tr1_mod.bed stdout | genePredToGtf file stdin stdout -honorCdsStat | sed "s/stdin/Gene/" > hg19_ucsc_clusters_new_tr1_mod.gtf
perl name_gtf.pl hg19_kgXref hg19_ucsc_clusters_new_tr1_mod.gtf > hg19_ucsc_final.gtf
awk -F "gene_id" '{print $2}' hg19_ucsc_final.gtf | sed "s/\"//g;s/ //g;s/;//g;s/transcript_id/\t/" | sort | uniq > hg19_gene_transcript_ref.rsem
egrep -v "random|chrUn|hap" hg19_ucsc_final.gtf > truechr_hg19_ucsc_final.gtf
awk -F "gene_id" '{print $2}' truechr_hg19_ucsc_final.gtf | sed "s/\"//g;s/ //g;s/;//g;s/transcript_id/\t/" | sort | uniq > truechr_hg19_gene_transcript_ref.rsem
module load RSEM/1.2.11
rsem-prepare-reference --gtf truechr_hg19_ucsc_final.gtf --transcript-to-gene-map truechr_hg19_gene_transcript_ref.rsem --bowtie2 --bowtie2-path /project/umw_garberlab/bin/bowtie2-2.1.0 ../human/hg19/hg19.fa ./hg19/hg19_ref_rsem
