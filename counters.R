# find the excel file with the list of experiments
# use vmatch pattern
# count the patterns per experiment

# get hairpin data
library(Biostrings)
require(ShortRead)

hairpin.path = "D:/devel/data/hairpin.fa"
hairpinSeqs <- readRNAStringSet(file=hairpin.path, format="fasta")

cElegHairpings <- hairpinSeqs[grep("^cel-", names(hairpinSeqs))] 

#setwd("D:/devel/data/small_RNA")
#file <- "B470002.1.fq"

getCounts <- function(file, occurences) {
    fq <- readFastq("./rna_no_adapters/", pattern = file)
    top = tables(sread(fq), n = length(unique(sread(fq))))
    # get only the sequences with at least this number of occurences
    experimentSeqs = top$top[top$top >= occurences]
    vcount <- rep(0, length(experimentSeqs))
    names(vcount) <- names(experimentSeqs)
    for (key in 1:length(experimentSeqs)) {
        seq <- names(vcount)[key]
        vcount[key] <- sum(vcountPattern(RNAString(DNAString(seq)), cElegHairpings, max.mismatch=2))
        if (key %% 1000 == 0) {
            print(key)
            flush.console()
        }
    }
    vcount
}

files <- list.files()
for(file in files) {
    vcount <- getCounts(file, occurences = 5)
    # store the counters
    dput(vcount[vcount > 0], paste("counters/", file, sep=""))
}


# normalize table counts using PreprocessCore library (see normalize.quantiles)


