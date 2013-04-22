# find the excel file with the list of experiments
# count the patterns per experiment

folder.project <- "/Users/radek/reactor/rnaoutliers"
source(paste(folder.project, "config.R", sep="/"))

# get hairpin data
hairpinSeqs <- readRNAStringSet(file=hairpin.path, format="fasta")

# get only C. Elegans data
cElegHairpings <- hairpinSeqs[grep("^cel-", names(hairpinSeqs))] 

#setwd("D:/devel/data/small_RNA")
#file <- "B470002.1.fq"

# get counts for all files within a given condition
getCounts <- function(files, occurences, path) {
    condCounts <- list()
    top_m <- NULL
    for (file in files) {
        fq <- readFastq(path, pattern = file)
        top = tables(sread(fq), n = length(unique(sread(fq))))
        condCounts[[file]] <- top$top
        row_m <- matrix(top$top, nrow=length(top$top), ncol=1)
        colnames(row_m) <- file
        rownames(row_m) <- names(top$top)
        if (is.null(top_m)) {
            top_m <- row_m
        } else {
            top_m <- cbindX(top_m, row_m)
        }
    }
    # filter out all sequences that has less than mean number of occurences
    top_m[rowSums(top_m) > occurences,])

    for (file in files) {
        fq <- readFastq(path, pattern = file)
        top = tables(sread(fq), n = length(unique(sread(fq))))
        experimentSeqs = top$top[top$top >= occurences]
        vcount <- rep(0, length(experimentSeqs))
        names(vcount) <- names(experimentSeqs)
        for (key in 1:length(experimentSeqs)) {
            seq <- names(vcount)[key]
            vcount[key] <- sum(vcountPattern(RNAString(DNAString(seq)), 
                cElegHairpings, max.mismatch=2))
            if (key %% 1000 == 0) {
                print(key)
                flush.console()
            }
        }
        vcount
    }

}

# for every condition
conditions <- unique(samples$condition)
# for testing use just one condition
condition <- "N2C"
for (condition in conditions) {
    files <- as.character(samples[samples$condition == condition,4])
    condCounts <- list()
    for (file in files) {
        #files <- list.files(folder.clean)
        #for(file in files) {
        vcount <- getCounts(file, occurences = 5, folder.clean)
        condCounts[[file]] <- vcount
        # TODO: only discard counters if they have less than 5 occurences 
        # across the conditions
        #vcount.positive <- vcount[vcount > 0]
        #dput(vcount.positive, paste(folder.counters, file, sep=""))
    }
}

files <- list.files(folder.counters)
count.m <- NULL
for (file in files) {
    vcount.positive <- dget(paste(folder.counters, file, sep=""))

    row.m <- matrix(vcount.positive, nrow=length(vcount.positive), ncol=1)
    colnames(row.m) <- file
    rownames(row.m) <- names(vcount.positive)

    if (is.null(count.m)) {
        count.m <- row.m
    } else {
        count.m <- cbindX(count.m, row.m)
    }
}
count.m <- na.omit(count.m)

# normalize table counts using PreprocessCore library (see normalize.quantiles)
# is this correct?
count.q <- normalize.quantiles(count.m)
rownames(count.q) <- rownames(count.m)
colnames(count.q) <- colnames(count.m)

# 

