#!/usr/bin/perl

my $db=$ARGV[0];
my $version=$ARGV[1];

my $dir="/share/data/umw_biocore/genome_data";
`mkdir -p $dir/$db/$version/make_$version`;
my $fine=1;
while ($fine)
{
 `./prepRSEM.bash $db $version>$dir/$db/$version/make_$version/res.txt 2>&1`;
 $a=`grep N[MR] $dir/$db/$version/make_$version/res.txt`;
 if ($a=~/Transcript\s(N[MR]_\d+)\sis/ || $a=~/cannot find (N[MR]_.+)\'s gene_id/)
 {
   print $1."\n";
   `grep -v $1 $dir/$db/$version/ucsc.gtf > $dir/$db/$version/ucs1c.gtf`;
   `mv $dir/$db/$version/ucs1c.gtf $dir/$db/$version/ucsc.gtf`;
 }
 else
 {
   $fine=0;
 }
}
