## Configuration file for rnaseq_pipeline.sh
##
## Place this script in a working directory and edit it accordingly.
##
## The default configuration assumes that the user unpacked the 
## chrX_data.tar.gz file in the current directory, so all the input
## files can be found in a ./chrX_data sub-directory

#how many CPUs to use on the current machine?
NUMCPUS=8

#### Program paths ####

## optional BINDIR, using it here because these programs are installed in a common directory
#BINDIR=/usr/local/bin
#HISAT2=$BINDIR/hisat2
#STRINGTIE=$BINDIR/stringtie
#SAMTOOLS=$BINDIR/samtools

#if these programs are not in any PATH directories, please edit accordingly:
HISAT2=$(which hisat2)
STRINGTIE=$(which stringtie)
SAMTOOLS=$(which samtools)

#### File paths for input data
### Full absolute paths are strongly recommended here.
## Warning: if using relatives paths here, these will be interpreted 
## relative to the  chosen output directory (which is generally the 
## working directory where this script is, unless the optional <output_dir>
## parameter is provided to the main pipeline script)

## Optional base directory, if most of the input files have a common path
# BASEDIR="/home/johnq/RNAseq_protocol/chrX_data"
BASEDIR=$(pwd -P)/chrX_data

FASTQLOC="$BASEDIR/samples"
GENOMEIDX="$BASEDIR/indexes/chrX_tran"
GTFFILE="$BASEDIR/genes/chrX.gtf"
PHENODATA="$BASEDIR/geuvadis_phenodata.csv"

TEMPLOC="./tmp" #this will be relative to the output directory

## list of samples 
## (only paired reads, must follow _1.*/_2.* file naming convention)
reads1=(${FASTQLOC}/*_1.*)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1./_2.}")
