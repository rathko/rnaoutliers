# find the excel file with the list of experiments
# use vmatch pattern
# count the patterns per experiment

# get hairpin data
hairpin.path = "D:/devel/data/hairpin.fa"
hairpinSeqs <- readRNAStringSet(file=hairpin.path, format="fasta")

cElegHairpings <- hairpinSeqs[grep("^cel-", names(hairpinSeqs))] 

#setwd("C:/devel/data/small_RNA")
file <- "B470001.1.fq"
fq <- readFastq("./output/", pattern = file)
experimentSeqs <- RNAStringSet(sread(fq))


#pdict <- PDict(DNAStringSet(experimentSeqs))
vcount <- rep(0, length(experimentSeqs))
names(vcount) <- as.character(experimentSeqs)
for (key in 1:length(experimentSeqs)) {
#for (key in 1:100) {
    #print(key)
    seq <- names(vcount)[key]
    #print(seq)
    vcount[key] <- sum(vcountPattern(seq, cElegHairpings, max.mismatch=2))
    print(vcount[key])
}

# how to speed up the above?