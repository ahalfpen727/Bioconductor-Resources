## updating derfinder R package:
## read things directly from bam file
## rather than creating the text files and such.
## FIRST ATTEMPT
## 11/20/12

## updated 12/11/12

library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)

# chromosome lengths:
#http://www.ncrna.org/glocal/cgi-bin/hgTracks?chromInfoPage=
# fix this later
# I think you can actually also get this from the sam header

## make a GRanges object representing each bp
## chr should be the chromosome as represented in your bam file.
## for now, works on our test bam file.
makeTabSkeleton = function(chr){
  chrlengths = c(249250621, 243199373, 198022430, 191154276, 180915260,
      171115067, 159138663, 146364022, 141213431, 135534747, 135006516,
      133851895, 115169878, 107349540, 102531392, 90354753, 81195210,
      78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
  if(is.character(chr)){
    if(chr!="X" & chr!="Y") stop("chromosome must be 1-22, X, or Y")
    if(chr=="X") chrnum=23
    if(chr=="Y") chrnum=24
  }
  if(is.numeric(chr)) chrnum=chr
  maxpos = chrlengths[chrnum]
  chr.by.bp = GRanges(seqnames = Rle(chr), ranges = IRanges(start = c(1:(maxpos-1)), end = c(2:maxpos)) )
  return(chr.by.bp)
}

#source("~/Google Drive/hopkins/research/memorycheck.R")
#ch22.by.bp = makeTabSkeleton(22)

#getCoverage = function(bamfile, tabskeleton){
#  aln = readGappedAlignments(bamfile)
#  bpcounts = countOverlaps(tabskeleton, aln)
#  return(bpcounts)
#}

#toybam = "~/Google Drive/hopkins/research/_toRNAdo-project/_Rpackage/tests/toyfile.ba
## ok this works.  I think this is going to be a horrible memory-hog though.  could we loop through, do all the samples, and also do the t-test??


getCoverage <- function(sampbams, sampnames, countfname, chr, cutoff = 10,
                           nzmeds = FALSE, chunk = FALSE, chunksize = NULL, group, adjustvars = NULL, scalefac = 32){
  numsamps = length(sampbams)
  srles = list()
  alns = list()
  for(i in 1:numsamps){
    alns[[i]] = readGappedAlignments(sampbams[i])
  }
  maxpos = as.numeric(seqlengths(alns[[1]])[which(names(seqlengths(alns[[1]])) == chr)])
  #passfilt = rep(0,length(posrange))
  countfile = file(countfname,"w")
  cat(paste("pos",paste(sampnames,sep="\t"), sep="\t"),file=countfile)
  cat("\n",file=countfile)
  pos = 1
  while(pos <= maxpos){
    posgranges = GRanges(seqnames = Rle(chr), ranges = IRanges(start = pos, end = pos+1))
    sampcounts = NULL
    for(i in 1:length(sampbams)){
      sampcounts[i] = countOverlaps(posgranges, alns[[i]])
    }

    ## only write this row out if it passes the filter.
    if(median(sampcounts)>=cutoff){
      cat(paste(pos,paste(sampcounts,collapse="\t"),sep="\t"),file=countfile)
      cat("\n",file=countfile)
      #passfilt[pos] = 1
    }

    ## update each sample's RLE object:
    for(s in 1:numsamps){
      if(pos==1){
        srles[[s]] = rle(sampcounts[s])
      }
      if(pos>1){
        if(sampcounts[s] == srles[[s]]$values[length(srles[[s]]$values)]){
          srles[[s]]$lengths[length(srles[[s]]$lengths)] <- srles[[s]]$lengths[length(srles[[s]]$lengths)] + 1
        }
        if(sampcounts[s] != srles[[s]]$values[length(srles[[s]]$values)]){
          srles[[s]]$lengths = c(srles[[s]]$lengths, 1)
          srles[[s]]$values = c(srles[[s]]$values, sampcounts[s])
        }
      } #end pos>1 case

    } #end updating RLEs
    pos = pos+1
  } #end while loop (over pos)

  close(countfile)

  sampmeds = NULL
  for(s in 1:numsamps){
    rle.ord = list(lengths = srles[[s]]$lengths[order(srles[[s]]$values)], values = sort(srles[[s]]$values))
    medind = which(cumsum(rle.ord$lengths) >= sum(rle.ord$lengths)/2)[1]
    sampmeds[s] = rle.ord$values[medind]
  }


  if(!chunk){
    counts = read.table(countfile, sep="\t", header=T)
    if(!is.null(adjustvars)){
      string1 = ""
      for(i in 1:dim(adjustvars)[2]){
        eval(parse(text=paste("av",i," <- adjustvars[,",i,"]",sep="")))
        string1 = paste(string1, paste("av",i,sep=""),sep="+")
      }
      eval(parse(text=paste("x = model.matrix(~group+sampmeds",string1,")",sep="")))
    } #end adjustvars case
    if(is.null(adjustvars)) x = model.matrix(~group+sampmeds)
    mymat = counts[,-1] #remove pos
    mymat = log2(mymat+scalefac)
    Amean = rowMeans(mymat)
    fit = lmFit(mymat, x)
    return(list(ebobject=list(coefficients=fit$coefficients, stdev.unscaled = fit$stdev.unscaled, sigma = fit$sigma, df.residual = fit$df.residual, Amean = Amean), pos=counts$pos))
  } #end function if not chunking.


} #end function




### test ###
setwd("/amber2/scratch/jleek/orbFrontal/results/testbams/")
sampbams = c("orbFrontalF11-small.bam", "orbFrontalF1-small.bam", "orbFrontalF23-small.bam",
  "orbFrontalF2-small.bam", "orbFrontalF32-small.bam", "orbFrontalF33-small.bam",
  "orbFrontalF3-small.bam", "orbFrontalF40-small.bam", "orbFrontalF42-small.bam",
  "orbFrontalF43-small.bam", "orbFrontalF47-small.bam", "orbFrontalF53-small.bam",
  "orbFrontalF55-small.bam", "orbFrontalF56-small.bam", "orbFrontalF58-small.bam")
sampnames = unlist(lapply(sampbams, function(x) substr(x,1,nchar(x)-10)))
maxpos = 249250621
countfname = "chr1-count-table.txt"
chr = 1
group = c(1,1,0,0,1,0,1,0,1,1,1,1,0,0,1) #sex labels

k = getCoverage(sampbams, sampnames, countfname, chr=1, group=group)
save(k,file="k.rda")





