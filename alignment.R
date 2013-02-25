
# another in series of RNA tutorials
# based on: http://www.bioconductor.org/help/course-materials/2010/BioC2010/Exercises-SimpleRNAseqUseCase.pdf
# http://www.bioconductor.org/help/course-materials/2012/useR2012/Bioconductor-tutorial.pdf
# very informative:
# http://www.bioconductor.org/help/course-materials/2012/Bressanone2012/2012-07-03-Anders-Sequencing-alignment.pdf

source("util.R")

library(Biostrings)

#fq <- readFastq(".", pattern = "B470001.1.fq")
#detail(fq)

library(Biostrings)

#setwd("/Users/radek/reactor/phd/small_RNA/clean")
file <- "B470001.1.fq"
fqClean <- readFastq(".", pattern = file)

experimentSeqs <- RNAStringSet(sread(fqClean))
# filter the records down only to cel (C. Elegans)
# how to loop through RNAStringSet?
for (seq in experimentSeqs) {
    print(seq)
    break
}

# get hairpin data
hairpinSeqs <- read.RNAStringSet(file=path.seqs, format="fasta")
hairpinNames <- names(hairpinSeqs)

# filter down to only c. elegans


# seq matching using vmatchPattern
for (seq in as.character(experimentSeqs)) {
    # get only > 18 long sequences
    mindex <- vmatchPattern(seq, hairpinSeqs, max.mismatch=2)
    # use vcountPattern instead - since we are interested inly in counts anyway
    len <- length(unlist(mindex))
    if (len > 0) {
        print(seq)
        print(len)
        #print(sum(countIndex(mindex)))
        if (len > 5) len <- 5
        print(unlist(mindex)[1:len])
        # check the counters?
    }
    #break;
}

# 1. ensure we count the records correctly
# 2. 

# store the counts


