This directory contains the Dec. 2013 (GRCh38/hg38) assembly of the
    human genome (hg38, GRCh38 Genome Reference Consortium Human Reference 38 (GCA_000001405.2)),
    as well as repeat annotations and GenBank sequences.

For more information about this assembly, please note the NCBI resources:
    http://www.ncbi.nlm.nih.gov/genome/51
    http://www.ncbi.nlm.nih.gov/genome/assembly/883148
    http://www.ncbi.nlm.nih.gov/bioproject/31257

Please note the extra 'analysis set' sequences as
described in: analysisSet/README.txt

Files included in this directory:

hg38.2bit - contains the complete human/hg38 genome sequence
    in the 2bit file format.  Repeats from RepeatMasker and Tandem Repeats
    Finder (with period of 12 or less) are shown in lower case; non-repeating
    sequence is shown in upper case.  The utility program, twoBitToFa (available
    from the kent src tree), can be used to extract .fa file(s) from
    this file.  A pre-compiled version of the command line tool can be
    found at:
        http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
    See also:
        http://genome.ucsc.edu/admin/git.html
	http://genome.ucsc.edu/admin/jk-install.html

hg38.agp.gz - Description of how the assembly was generated from
    fragments.

hg38.chromFa.tar.gz - The assembly sequence in one file per chromosome.
    Repeats from RepeatMasker and Tandem Repeats Finder (with period
    of 12 or less) are shown in lower case; non-repeating sequence is
    shown in upper case.

hg38.chromFaMasked.tar.gz - The assembly sequence in one file per chromosome.
    Repeats are masked by capital Ns; non-repeating sequence is shown in
    upper case.

hg38.fa.gz - "Soft-masked" assembly sequence in one file.
    Repeats from RepeatMasker and Tandem Repeats Finder (with period
    of 12 or less) are shown in lower case; non-repeating sequence is
    shown in upper case.

hg38.fa.masked.gz - "Hard-masked" assembly sequence in one file.
    Repeats are masked by capital Ns; non-repeating sequence is shown in
    upper case.

hg38.fa.out.gz - RepeatMasker .out file.  RepeatMasker was run with the
    -s (sensitive) setting.
    June 20 2013 (open-4-0-3) version of RepeatMasker
    RepBase library: RELEASE 20130422

hg38.fa.align.gz - RepeatMasker .align file.  RepeatMasker was run with the
    -s (sensitive) setting.
    June 20 2013 (open-4-0-3) version of RepeatMasker
    RepBase library: RELEASE 20130422

hg38.trf.bed.gz - Tandem Repeats Finder locations, filtered to keep repeats
    with period less than or equal to 12, and translated into UCSC's BED
    format.

md5sum.txt - checksums of files in this directory

mrna.fa.gz - Human mRNA from GenBank. This sequence data is updated
    once a week via automatic GenBank updates.

refMrna.fa.gz - RefSeq mRNA from the same species as the genome.
    This sequence data is updated once a week via automatic GenBank
    updates.

upstream1000.fa.gz - Sequences 1000 bases upstream of annotated
    transcription starts of RefSeq genes with annotated 5' UTRs.
    This file is updated weekly so it might be slightly out of sync with
    the RefSeq data which is updated daily for most assemblies.

upstream2000.fa.gz - Same as upstream1000, but 2000 bases.

upstream5000.fa.gz - Same as upstream1000, but 5000 bases.

xenoMrna.fa.gz - GenBank mRNAs from species other than that of
    the genome. This sequence data is updated once a week via automatic
    GenBank updates.


hg38.chrom.sizes - Two-column tab-separated text file containing assembly
    sequence names and sizes.

------------------------------------------------------------------
If you plan to download a large file or multiple files from this
directory, we recommend that you use ftp rather than downloading the
files via our website. To do so, ftp to hgdownload.cse.ucsc.edu
[username: anonymous, password: your email address], then cd to the
directory goldenPath/hg38/bigZips. To download multiple files, use
the "mget" command:

    mget <filename1> <filename2> ...
    - or -
    mget -a (to download all the files in the directory)

Alternate methods to ftp access.

Using an rsync command to download the entire directory:
    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/ .
For a single file, e.g. chromFa.tar.gz
    rsync -avzP 
        rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/chromFa.tar.gz .

Or with wget, all files:
    wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/*'
With wget, a single file:
    wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/chromFa.tar.gz' 
        -O chromFa.tar.gz

To unpack the *.tar.gz files:
    tar xvzf <file>.tar.gz
To uncompress the fa.gz files:
    gunzip <file>.fa.gz

-----------------------------------------------------------------------------
GenBank Data Usage

The GenBank database is designed to provide and encourage access within
the scientific community to the most up to date and comprehensive DNA
sequence information. Therefore, NCBI places no restrictions on the use
or distribution of the GenBank data. However, some submitters may claim
patent, copyright, or other intellectual property rights in all or a
portion of the data they have submitted. NCBI is not in a position to
assess the validity of such claims, and therefore cannot provide comment
or unrestricted permission concerning the use, copying, or distribution
of the information contained in GenBank.
-----------------------------------------------------------------------------
