Full genome sequences for Homo sapiens (UCSC version hg38)
Description

Full genome sequences for Homo sapiens (Human) as provided by UCSC (hg38, Dec. 2013) and stored in Biostrings objects.
Note

This BSgenome data package was made from the following source data files:

hg38.2bit from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/
  

See ?BSgenomeForge and the BSgenomeForge vignette (vignette("BSgenomeForge")) in the BSgenome software package for how to make a BSgenome data package.
Author(s)

The Bioconductor Dev Team
See Also

    BSgenome objects and the available.genomes function in the BSgenome software package.

    DNAString objects in the Biostrings package.

    The BSgenomeForge vignette (vignette("BSgenomeForge")) in the BSgenome software package for how to make a BSgenome data package.

Examples

BSgenome.Hsapiens.UCSC.hg38
genome <- BSgenome.Hsapiens.UCSC.hg38
seqlengths(genome)
genome$chr1  # same as genome[["chr1"]]

## ---------------------------------------------------------------------
## Extract the upstream sequences
## ---------------------------------------------------------------------
## The upstream sequences located in
##   http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/
## are based on RefSeq genes (RefSeq Genes track in the Genome Browser).
## Upstream sequences based on UCSC genes (UCSC Genes track in the
## Genome Browser) can easily be extracted from the full genome
## sequences with:

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
knownGene_txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
knownGene_up1000seqs <- extractUpstreamSeqs(genome, knownGene_txdb)

## Or, to get upstream sequences based on RefSeq genes:

refGene_txdb <- makeTxDbFromUCSC("hg38", "refGene")
refGene_up1000seqs <- extractUpstreamSeqs(genome, refGene_txdb)

## Note that you can make a TxDb object from various annotation
## resources. See the makeTxDbFromUCSC(), makeTxDbFromBiomart(), and
## makeTxDbFromGFF() functions in the GenomicFeatures package for more
## information.
## IMPORTANT: Make sure you use a TxDb package (or TxDb object) that
## contains a gene model based on hg38 or on a compatible genome (i.e.
## a genome with sequences identical to the sequences in hg38). See
## ?extractUpstreamSeqs in the GenomicFeatures package for more
## information.

## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------
## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")

Results


R version 3.3.1 (2016-06-21) -- "Bug in Your Hair"
Copyright (C) 2016 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> library(BSgenome.Hsapiens.UCSC.hg38)
Loading required package: BSgenome
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: 'BiocGenerics'

The following objects are masked from 'package:parallel':

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from 'package:stats':

    IQR, mad, xtabs

The following objects are masked from 'package:base':

    Filter, Find, Map, Position, Reduce, anyDuplicated, append,
    as.data.frame, cbind, colnames, do.call, duplicated, eval, evalq,
    get, grep, grepl, intersect, is.unsorted, lapply, lengths, mapply,
    match, mget, order, paste, pmax, pmax.int, pmin, pmin.int, rank,
    rbind, rownames, sapply, setdiff, sort, table, tapply, union,
    unique, unsplit

Loading required package: S4Vectors
Loading required package: stats4

Attaching package: 'S4Vectors'

The following objects are masked from 'package:base':

    colMeans, colSums, expand.grid, rowMeans, rowSums

Loading required package: IRanges
Loading required package: GenomeInfoDb
Loading required package: GenomicRanges
Loading required package: Biostrings
Loading required package: XVector
Loading required package: rtracklayer
> png(filename="/home/ddbj/snapshot/RGM3/R_BC/result/BSgenome.Hsapiens.UCSC.hg38/package.Rd_%03d_medium.png", width=480, height=480)
> ### Name: BSgenome.Hsapiens.UCSC.hg38
> ### Title: Full genome sequences for Homo sapiens (UCSC version hg38)
> ### Aliases: BSgenome.Hsapiens.UCSC.hg38-package
> ###   BSgenome.Hsapiens.UCSC.hg38 Hsapiens
> ### Keywords: package data
> 
> ### ** Examples
> 
> BSgenome.Hsapiens.UCSC.hg38
Human genome:
# organism: Homo sapiens (Human)
# provider: UCSC
# provider version: hg38
# release date: Dec. 2013
# release name: Genome Reference Consortium GRCh38
# 455 sequences:
#   chr1                    chr2                    chr3                   
#   chr4                    chr5                    chr6                   
#   chr7                    chr8                    chr9                   
#   chr10                   chr11                   chr12                  
#   chr13                   chr14                   chr15                  
#   ...                     ...                     ...                    
#   chrUn_KI270744v1        chrUn_KI270745v1        chrUn_KI270746v1       
#   chrUn_KI270747v1        chrUn_KI270748v1        chrUn_KI270749v1       
#   chrUn_KI270750v1        chrUn_KI270751v1        chrUn_KI270752v1       
#   chrUn_KI270753v1        chrUn_KI270754v1        chrUn_KI270755v1       
#   chrUn_KI270756v1        chrUn_KI270757v1                               
# (use 'seqnames()' to see all the sequence names, use the '$' or '[[' operator
# to access a given sequence)
> genome <- BSgenome.Hsapiens.UCSC.hg38
> seqlengths(genome)
                   chr1                    chr2                    chr3 
              248956422               242193529               198295559 
                   chr4                    chr5                    chr6 
              190214555               181538259               170805979 
                   chr7                    chr8                    chr9 
              159345973               145138636               138394717 
                  chr10                   chr11                   chr12 
              133797422               135086622               133275309 
                  chr13                   chr14                   chr15 
              114364328               107043718               101991189 
                  chr16                   chr17                   chr18 
               90338345                83257441                80373285 
                  chr19                   chr20                   chr21 
               58617616                64444167                46709983 
                  chr22                    chrX                    chrY 
               50818468               156040895                57227415 
                   chrM     chr1_GL383518v1_alt     chr1_GL383519v1_alt 
                  16569                  182439                  110268 
    chr1_GL383520v2_alt     chr1_KI270759v1_alt     chr1_KI270760v1_alt 
                 366580                  425601                  109528 
    chr1_KI270761v1_alt     chr1_KI270762v1_alt     chr1_KI270763v1_alt 
                 165834                  354444                  911658 
    chr1_KI270764v1_alt     chr1_KI270765v1_alt     chr1_KI270766v1_alt 
                  50258                  185285                  256271 
    chr1_KI270892v1_alt     chr2_GL383521v1_alt     chr2_GL383522v1_alt 
                 162212                  143390                  123821 
    chr2_GL582966v2_alt     chr2_KI270767v1_alt     chr2_KI270768v1_alt 
                  96131                  161578                  110099 
    chr2_KI270769v1_alt     chr2_KI270770v1_alt     chr2_KI270771v1_alt 
                 120616                  136240                  110395 
    chr2_KI270772v1_alt     chr2_KI270773v1_alt     chr2_KI270774v1_alt 
                 133041                   70887                  223625 
    chr2_KI270775v1_alt     chr2_KI270776v1_alt     chr2_KI270893v1_alt 
                 138019                  174166                  161218 
    chr2_KI270894v1_alt     chr3_GL383526v1_alt     chr3_JH636055v2_alt 
                 214158                  180671                  173151 
    chr3_KI270777v1_alt     chr3_KI270778v1_alt     chr3_KI270779v1_alt 
                 173649                  248252                  205312 
    chr3_KI270780v1_alt     chr3_KI270781v1_alt     chr3_KI270782v1_alt 
                 224108                  113034                  162429 
    chr3_KI270783v1_alt     chr3_KI270784v1_alt     chr3_KI270895v1_alt 
                 109187                  184404                  162896 
    chr3_KI270924v1_alt     chr3_KI270934v1_alt     chr3_KI270935v1_alt 
                 166540                  163458                  197351 
    chr3_KI270936v1_alt     chr3_KI270937v1_alt     chr4_GL000257v2_alt 
                 164170                  165607                  586476 
    chr4_GL383527v1_alt     chr4_GL383528v1_alt     chr4_KI270785v1_alt 
                 164536                  376187                  119912 
    chr4_KI270786v1_alt     chr4_KI270787v1_alt     chr4_KI270788v1_alt 
                 244096                  111943                  158965 
    chr4_KI270789v1_alt     chr4_KI270790v1_alt     chr4_KI270896v1_alt 
                 205944                  220246                  378547 
    chr4_KI270925v1_alt     chr5_GL339449v2_alt     chr5_GL383530v1_alt 
                 555799                 1612928                  101241 
    chr5_GL383531v1_alt     chr5_GL383532v1_alt     chr5_GL949742v1_alt 
                 173459                   82728                  226852 
    chr5_KI270791v1_alt     chr5_KI270792v1_alt     chr5_KI270793v1_alt 
                 195710                  179043                  126136 
    chr5_KI270794v1_alt     chr5_KI270795v1_alt     chr5_KI270796v1_alt 
                 164558                  131892                  172708 
    chr5_KI270897v1_alt     chr5_KI270898v1_alt     chr6_GL000250v2_alt 
                1144418                  130957                 4672374 
    chr6_GL000251v2_alt     chr6_GL000252v2_alt     chr6_GL000253v2_alt 
                4795265                 4604811                 4677643 
    chr6_GL000254v2_alt     chr6_GL000255v2_alt     chr6_GL000256v2_alt 
                4827813                 4606388                 4929269 
    chr6_GL383533v1_alt     chr6_KB021644v2_alt     chr6_KI270758v1_alt 
                 124736                  185823                   76752 
    chr6_KI270797v1_alt     chr6_KI270798v1_alt     chr6_KI270799v1_alt 
                 197536                  271782                  152148 
    chr6_KI270800v1_alt     chr6_KI270801v1_alt     chr6_KI270802v1_alt 
                 175808                  870480                   75005 
    chr7_GL383534v2_alt     chr7_KI270803v1_alt     chr7_KI270804v1_alt 
                 119183                 1111570                  157952 
    chr7_KI270805v1_alt     chr7_KI270806v1_alt     chr7_KI270807v1_alt 
                 209988                  158166                  126434 
    chr7_KI270808v1_alt     chr7_KI270809v1_alt     chr7_KI270899v1_alt 
                 271455                  209586                  190869 
    chr8_KI270810v1_alt     chr8_KI270811v1_alt     chr8_KI270812v1_alt 
                 374415                  292436                  282736 
    chr8_KI270813v1_alt     chr8_KI270814v1_alt     chr8_KI270815v1_alt 
                 300230                  141812                  132244 
    chr8_KI270816v1_alt     chr8_KI270817v1_alt     chr8_KI270818v1_alt 
                 305841                  158983                  145606 
    chr8_KI270819v1_alt     chr8_KI270820v1_alt     chr8_KI270821v1_alt 
                 133535                   36640                  985506 
    chr8_KI270822v1_alt     chr8_KI270900v1_alt     chr8_KI270901v1_alt 
                 624492                  318687                  136959 
    chr8_KI270926v1_alt     chr9_GL383539v1_alt     chr9_GL383540v1_alt 
                 229282                  162988                   71551 
    chr9_GL383541v1_alt     chr9_GL383542v1_alt     chr9_KI270823v1_alt 
                 171286                   60032                  439082 
   chr10_GL383545v1_alt    chr10_GL383546v1_alt    chr10_KI270824v1_alt 
                 179254                  309802                  181496 
   chr10_KI270825v1_alt    chr11_GL383547v1_alt    chr11_JH159136v1_alt 
                 188315                  154407                  200998 
   chr11_JH159137v1_alt    chr11_KI270826v1_alt    chr11_KI270827v1_alt 
                 191409                  186169                   67707 
   chr11_KI270829v1_alt    chr11_KI270830v1_alt    chr11_KI270831v1_alt 
                 204059                  177092                  296895 
   chr11_KI270832v1_alt    chr11_KI270902v1_alt    chr11_KI270903v1_alt 
                 210133                  106711                  214625 
   chr11_KI270927v1_alt    chr12_GL383549v1_alt    chr12_GL383550v2_alt 
                 218612                  120804                  169178 
   chr12_GL383551v1_alt    chr12_GL383552v1_alt    chr12_GL383553v2_alt 
                 184319                  138655                  152874 
   chr12_GL877875v1_alt    chr12_GL877876v1_alt    chr12_KI270833v1_alt 
                 167313                  408271                   76061 
   chr12_KI270834v1_alt    chr12_KI270835v1_alt    chr12_KI270836v1_alt 
                 119498                  238139                   56134 
   chr12_KI270837v1_alt    chr12_KI270904v1_alt    chr13_KI270838v1_alt 
                  40090                  572349                  306913 
   chr13_KI270839v1_alt    chr13_KI270840v1_alt    chr13_KI270841v1_alt 
                 180306                  191684                  169134 
   chr13_KI270842v1_alt    chr13_KI270843v1_alt    chr14_KI270844v1_alt 
                  37287                  103832                  322166 
   chr14_KI270845v1_alt    chr14_KI270846v1_alt    chr14_KI270847v1_alt 
                 180703                 1351393                 1511111 
   chr15_GL383554v1_alt    chr15_GL383555v2_alt    chr15_KI270848v1_alt 
                 296527                  388773                  327382 
   chr15_KI270849v1_alt    chr15_KI270850v1_alt    chr15_KI270851v1_alt 
                 244917                  430880                  263054 
   chr15_KI270852v1_alt    chr15_KI270905v1_alt    chr15_KI270906v1_alt 
                 478999                 5161414                  196384 
   chr16_GL383556v1_alt    chr16_GL383557v1_alt    chr16_KI270853v1_alt 
                 192462                   89672                 2659700 
   chr16_KI270854v1_alt    chr16_KI270855v1_alt    chr16_KI270856v1_alt 
                 134193                  232857                   63982 
   chr17_GL000258v2_alt    chr17_GL383563v3_alt    chr17_GL383564v2_alt 
                1821992                  375691                  133151 
   chr17_GL383565v1_alt    chr17_GL383566v1_alt    chr17_JH159146v1_alt 
                 223995                   90219                  278131 
   chr17_JH159147v1_alt    chr17_JH159148v1_alt    chr17_KI270857v1_alt 
                  70345                   88070                 2877074 
   chr17_KI270858v1_alt    chr17_KI270859v1_alt    chr17_KI270860v1_alt 
                 235827                  108763                  178921 
   chr17_KI270861v1_alt    chr17_KI270862v1_alt    chr17_KI270907v1_alt 
                 196688                  391357                  137721 
   chr17_KI270908v1_alt    chr17_KI270909v1_alt    chr17_KI270910v1_alt 
                1423190                  325800                  157099 
   chr18_GL383567v1_alt    chr18_GL383568v1_alt    chr18_GL383569v1_alt 
                 289831                  104552                  167950 
   chr18_GL383570v1_alt    chr18_GL383571v1_alt    chr18_GL383572v1_alt 
                 164789                  198278                  159547 
   chr18_KI270863v1_alt    chr18_KI270864v1_alt    chr18_KI270911v1_alt 
                 167999                  111737                  157710 
   chr18_KI270912v1_alt    chr19_GL000209v2_alt    chr19_GL383573v1_alt 
                 174061                  177381                  385657 
   chr19_GL383574v1_alt    chr19_GL383575v2_alt    chr19_GL383576v1_alt 
                 155864                  170222                  188024 
   chr19_GL949746v1_alt    chr19_GL949747v2_alt    chr19_GL949748v2_alt 
                 987716                  729520                 1064304 
   chr19_GL949749v2_alt    chr19_GL949750v2_alt    chr19_GL949751v2_alt 
                1091841                 1066390                 1002683 
   chr19_GL949752v1_alt    chr19_GL949753v2_alt    chr19_KI270865v1_alt 
                 987100                  796479                   52969 
   chr19_KI270866v1_alt    chr19_KI270867v1_alt    chr19_KI270868v1_alt 
                  43156                  233762                   61734 
   chr19_KI270882v1_alt    chr19_KI270883v1_alt    chr19_KI270884v1_alt 
                 248807                  170399                  157053 
   chr19_KI270885v1_alt    chr19_KI270886v1_alt    chr19_KI270887v1_alt 
                 171027                  204239                  209512 
   chr19_KI270888v1_alt    chr19_KI270889v1_alt    chr19_KI270890v1_alt 
                 155532                  170698                  184499 
   chr19_KI270891v1_alt    chr19_KI270914v1_alt    chr19_KI270915v1_alt 
                 170680                  205194                  170665 
   chr19_KI270916v1_alt    chr19_KI270917v1_alt    chr19_KI270918v1_alt 
                 184516                  190932                  123111 
   chr19_KI270919v1_alt    chr19_KI270920v1_alt    chr19_KI270921v1_alt 
                 170701                  198005                  282224 
   chr19_KI270922v1_alt    chr19_KI270923v1_alt    chr19_KI270929v1_alt 
                 187935                  189352                  186203 
   chr19_KI270930v1_alt    chr19_KI270931v1_alt    chr19_KI270932v1_alt 
                 200773                  170148                  215732 
   chr19_KI270933v1_alt    chr19_KI270938v1_alt    chr20_GL383577v2_alt 
                 170537                 1066800                  128386 
   chr20_KI270869v1_alt    chr20_KI270870v1_alt    chr20_KI270871v1_alt 
                 118774                  183433                   58661 
   chr21_GL383578v2_alt    chr21_GL383579v2_alt    chr21_GL383580v2_alt 
                  63917                  201197                   74653 
   chr21_GL383581v2_alt    chr21_KI270872v1_alt    chr21_KI270873v1_alt 
                 116689                   82692                  143900 
   chr21_KI270874v1_alt    chr22_GL383582v2_alt    chr22_GL383583v2_alt 
                 166743                  162811                   96924 
   chr22_KB663609v1_alt    chr22_KI270875v1_alt    chr22_KI270876v1_alt 
                  74013                  259914                  263666 
   chr22_KI270877v1_alt    chr22_KI270878v1_alt    chr22_KI270879v1_alt 
                 101331                  186262                  304135 
   chr22_KI270928v1_alt     chrX_KI270880v1_alt     chrX_KI270881v1_alt 
                 176103                  284869                  144206 
    chrX_KI270913v1_alt  chr1_KI270706v1_random  chr1_KI270707v1_random 
                 274009                  175055                   32032 
 chr1_KI270708v1_random  chr1_KI270709v1_random  chr1_KI270710v1_random 
                 127682                   66860                   40176 
 chr1_KI270711v1_random  chr1_KI270712v1_random  chr1_KI270713v1_random 
                  42210                  176043                   40745 
 chr1_KI270714v1_random  chr2_KI270715v1_random  chr2_KI270716v1_random 
                  41717                  161471                  153799 
 chr3_GL000221v1_random  chr4_GL000008v2_random  chr5_GL000208v1_random 
                 155397                  209709                   92689 
 chr9_KI270717v1_random  chr9_KI270718v1_random  chr9_KI270719v1_random 
                  40062                   38054                  176845 
 chr9_KI270720v1_random chr11_KI270721v1_random chr14_GL000009v2_random 
                  39050                  100316                  201709 
chr14_GL000194v1_random chr14_GL000225v1_random chr14_KI270722v1_random 
                 191469                  211173                  194050 
chr14_KI270723v1_random chr14_KI270724v1_random chr14_KI270725v1_random 
                  38115                   39555                  172810 
chr14_KI270726v1_random chr15_KI270727v1_random chr16_KI270728v1_random 
                  43739                  448248                 1872759 
chr17_GL000205v2_random chr17_KI270729v1_random chr17_KI270730v1_random 
                 185591                  280839                  112551 
chr22_KI270731v1_random chr22_KI270732v1_random chr22_KI270733v1_random 
                 150754                   41543                  179772 
chr22_KI270734v1_random chr22_KI270735v1_random chr22_KI270736v1_random 
                 165050                   42811                  181920 
chr22_KI270737v1_random chr22_KI270738v1_random chr22_KI270739v1_random 
                 103838                   99375                   73985 
 chrY_KI270740v1_random        chrUn_GL000195v1        chrUn_GL000213v1 
                  37240                  182896                  164239 
       chrUn_GL000214v1        chrUn_GL000216v2        chrUn_GL000218v1 
                 137718                  176608                  161147 
       chrUn_GL000219v1        chrUn_GL000220v1        chrUn_GL000224v1 
                 179198                  161802                  179693 
       chrUn_GL000226v1        chrUn_KI270302v1        chrUn_KI270303v1 
                  15008                    2274                    1942 
       chrUn_KI270304v1        chrUn_KI270305v1        chrUn_KI270310v1 
                   2165                    1472                    1201 
       chrUn_KI270311v1        chrUn_KI270312v1        chrUn_KI270315v1 
                  12399                     998                    2276 
       chrUn_KI270316v1        chrUn_KI270317v1        chrUn_KI270320v1 
                   1444                   37690                    4416 
       chrUn_KI270322v1        chrUn_KI270329v1        chrUn_KI270330v1 
                  21476                    1040                    1652 
       chrUn_KI270333v1        chrUn_KI270334v1        chrUn_KI270335v1 
                   2699                    1368                    1048 
       chrUn_KI270336v1        chrUn_KI270337v1        chrUn_KI270338v1 
                   1026                    1121                    1428 
       chrUn_KI270340v1        chrUn_KI270362v1        chrUn_KI270363v1 
                   1428                    3530                    1803 
       chrUn_KI270364v1        chrUn_KI270366v1        chrUn_KI270371v1 
                   2855                    8320                    2805 
       chrUn_KI270372v1        chrUn_KI270373v1        chrUn_KI270374v1 
                   1650                    1451                    2656 
       chrUn_KI270375v1        chrUn_KI270376v1        chrUn_KI270378v1 
                   2378                    1136                    1048 
       chrUn_KI270379v1        chrUn_KI270381v1        chrUn_KI270382v1 
                   1045                    1930                    4215 
       chrUn_KI270383v1        chrUn_KI270384v1        chrUn_KI270385v1 
                   1750                    1658                     990 
       chrUn_KI270386v1        chrUn_KI270387v1        chrUn_KI270388v1 
                   1788                    1537                    1216 
       chrUn_KI270389v1        chrUn_KI270390v1        chrUn_KI270391v1 
                   1298                    2387                    1484 
       chrUn_KI270392v1        chrUn_KI270393v1        chrUn_KI270394v1 
                    971                    1308                     970 
       chrUn_KI270395v1        chrUn_KI270396v1        chrUn_KI270411v1 
                   1143                    1880                    2646 
       chrUn_KI270412v1        chrUn_KI270414v1        chrUn_KI270417v1 
                   1179                    2489                    2043 
       chrUn_KI270418v1        chrUn_KI270419v1        chrUn_KI270420v1 
                   2145                    1029                    2321 
       chrUn_KI270422v1        chrUn_KI270423v1        chrUn_KI270424v1 
                   1445                     981                    2140 
       chrUn_KI270425v1        chrUn_KI270429v1        chrUn_KI270435v1 
                   1884                    1361                   92983 
       chrUn_KI270438v1        chrUn_KI270442v1        chrUn_KI270448v1 
                 112505                  392061                    7992 
       chrUn_KI270465v1        chrUn_KI270466v1        chrUn_KI270467v1 
                   1774                    1233                    3920 
       chrUn_KI270468v1        chrUn_KI270507v1        chrUn_KI270508v1 
                   4055                    5353                    1951 
       chrUn_KI270509v1        chrUn_KI270510v1        chrUn_KI270511v1 
                   2318                    2415                    8127 
       chrUn_KI270512v1        chrUn_KI270515v1        chrUn_KI270516v1 
                  22689                    6361                    1300 
       chrUn_KI270517v1        chrUn_KI270518v1        chrUn_KI270519v1 
                   3253                    2186                  138126 
       chrUn_KI270521v1        chrUn_KI270522v1        chrUn_KI270528v1 
                   7642                    5674                    2983 
       chrUn_KI270529v1        chrUn_KI270530v1        chrUn_KI270538v1 
                   1899                    2168                   91309 
       chrUn_KI270539v1        chrUn_KI270544v1        chrUn_KI270548v1 
                    993                    1202                    1599 
       chrUn_KI270579v1        chrUn_KI270580v1        chrUn_KI270581v1 
                  31033                    1553                    7046 
       chrUn_KI270582v1        chrUn_KI270583v1        chrUn_KI270584v1 
                   6504                    1400                    4513 
       chrUn_KI270587v1        chrUn_KI270588v1        chrUn_KI270589v1 
                   2969                    6158                   44474 
       chrUn_KI270590v1        chrUn_KI270591v1        chrUn_KI270593v1 
                   4685                    5796                    3041 
       chrUn_KI270741v1        chrUn_KI270742v1        chrUn_KI270743v1 
                 157432                  186739                  210658 
       chrUn_KI270744v1        chrUn_KI270745v1        chrUn_KI270746v1 
                 168472                   41891                   66486 
       chrUn_KI270747v1        chrUn_KI270748v1        chrUn_KI270749v1 
                 198735                   93321                  158759 
       chrUn_KI270750v1        chrUn_KI270751v1        chrUn_KI270752v1 
                 148850                  150742                   27745 
       chrUn_KI270753v1        chrUn_KI270754v1        chrUn_KI270755v1 
                  62944                   40191                   36723 
       chrUn_KI270756v1        chrUn_KI270757v1 
                  79590                   71251 
> genome$chr1  # same as genome[["chr1"]]
  248956422-letter "DNAString" instance
seq: NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
> 
> ## ---------------------------------------------------------------------
> ## Extract the upstream sequences
> ## ---------------------------------------------------------------------
> ## The upstream sequences located in
> ##   http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/
> ## are based on RefSeq genes (RefSeq Genes track in the Genome Browser).
> ## Upstream sequences based on UCSC genes (UCSC Genes track in the
> ## Genome Browser) can easily be extracted from the full genome
> ## sequences with:
> 
> library(TxDb.Hsapiens.UCSC.hg38.knownGene)
Loading required package: GenomicFeatures
Loading required package: AnnotationDbi
Loading required package: Biobase
Welcome to Bioconductor

    Vignettes contain introductory material; view with
    'browseVignettes()'. To cite Bioconductor, see
    'citation("Biobase")', and for packages 'citation("pkgname")'.

> knownGene_txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
> knownGene_up1000seqs <- extractUpstreamSeqs(genome, knownGene_txdb)
> 
> ## Or, to get upstream sequences based on RefSeq genes:
> 
> refGene_txdb <- makeTxDbFromUCSC("hg38", "refGene")
Download the refGene table ... OK
Download the refLink table ... OK
Extract the 'transcripts' data frame ... OK
Extract the 'splicings' data frame ... OK
Download and preprocess the 'chrominfo' data frame ... OK
Prepare the 'metadata' data frame ... OK
Make the TxDb object ... OK
Warning message:
In .extractCdsLocsFromUCSCTxTable(ucsc_txtable, exon_locs) :
  UCSC data anomaly in 545 transcript(s): the cds cumulative length is
  not a multiple of 3 for transcripts 'NM_001305275' 'NM_017940'
  'NM_001289974' 'NM_001291281' 'NM_001134939' 'NM_001301371'
  'NM_016178' 'NM_001145051' 'NM_001128929' 'NM_001075' 'NM_001144767'
  'NM_001322371' 'NM_032470' 'NM_004197' 'NM_032454' 'NM_016098'
  'NM_001788' 'NM_001172437' 'NM_001184961' 'NM_015068' 'NM_001159995'
  'NM_001159999' 'NM_001160001' 'NM_001005336' 'NM_001288737'
  'NM_001288738' 'NM_001288739' 'NM_004408' 'NM_020469' 'NM_001001676'
  'NM_033380' 'NM_053005' 'NM_001013356' 'NM_173600' 'NM_006400'
  'NM_001130048' 'NM_001318849' 'NM_015296' 'NM_006220' 'NM_001282494'
  'NM_001282490' 'NM_001301302' 'NM_002537' 'NM_001278425' 'NM_052892'
  'NM_130464' 'NM_001277332' 'NM_182705' 'NM_001291471' 'NM_001291472'
  'NM_001291473' 'NM_001291474' 'NM_001291475' 'NM_001123392'
  'NM_001291462' 'NM_001291463' 'NM_001291465' 'NM_000068'
  'NM_001174080' 'NM_023035' 'NM_001736' 'NM_001301020' 'NM_00415 [... truncated]
> refGene_up1000seqs <- extractUpstreamSeqs(genome, refGene_txdb)
> 
> ## Note that you can make a TxDb object from various annotation
> ## resources. See the makeTxDbFromUCSC(), makeTxDbFromBiomart(), and
> ## makeTxDbFromGFF() functions in the GenomicFeatures package for more
> ## information.
> ## IMPORTANT: Make sure you use a TxDb package (or TxDb object) that
> ## contains a gene model based on hg38 or on a compatible genome (i.e.
> ## a genome with sequences identical to the sequences in hg38). See
> ## ?extractUpstreamSeqs in the GenomicFeatures package for more
> ## information.
> 
> ## ---------------------------------------------------------------------
> ## Genome-wide motif searching
> ## ---------------------------------------------------------------------
> ## See the GenomeSearching vignette in the BSgenome software
> ## package for some examples of genome-wide motif searching using
> ## Biostrings and the BSgenome data packages:
> #if (interactive())
>     vignette("GenomeSearching", package="BSgenome")
