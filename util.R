library(ShortRead)

#setwd("C:/devel/data/small_RNA")
#setwd("/Users/radek/reactor/phd/small_RNA")

data.dir <- "/Users/radek/reactor/phd/small_RNA"
require(ShortRead)

removeAdapters <- function(fq) {
    # remoe poly A, if any
    polyA   <- "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    polyADist   <- srdistance(sread(fq), polyA)[[1]]
    fqClean <- fq[polyADist >= 25]

    # Filter N's, these bases are not relevant when mapping the sequence anyway
    filterN <- nFilter(threshold=3)
    fqClean <- fqClean[filterN(sread(fqClean))]
    
    # Removing the adapter
    adapter <- "ATGGAATTCTCGGGTGCCA"
    adapterPad <- paste(substr(adapter,1,12), paste(rep("N", 25),collapse=""), sep="")
    maxMm <- round(seq_len(12)*0.2)
    maxMmPad <-  c(maxMm,rep(maxMm[length(maxMm)], 25))
    coordinates <- trimLRPatterns(Rpattern=adapterPad, subject=sequences, max.Rmismatch=maxMmPad, 
    		ranges=TRUE, with.Rindels=FALSE, Rfixed=FALSE)
    seqTrim <- DNAStringSet(sequences, start=start(coordinates), end=end(coordinates))
    qualTrim <- BStringSet(quality, start=start(coordinates), end=end(coordinates))

    qualTrim <- FastqQuality(qualTrim)
    fqTrim <- ShortReadQ(sread=seqTrim, quality=qualTrim, id=id(fqClean))
    fqTrim
}