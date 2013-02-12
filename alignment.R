
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

####

# assess quality of files

#fastqDir <- file.path(data.dir, "fq")
fastqDir <- data.dir
fastqFiles <- dir(fastqDir, full=TRUE)

qas0 <- Map(function(fl, nm) {
        fq <- FastqSampler(fl)
        qa(yield(fq), nm)
}, fastqFiles, sub("_subset.fastq", "", basename(fastqFiles)))
qas <- do.call(rbind, qas0)
rpt <- report(qas, dest=tempfile())
browseURL(rpt)

#A report from a larger subset of the experiment is available
#rpt <- system.file("GSM461176_81_qa_report", "index.html", package="useR2012")
#browseURL(rpt)

# alignemnt with positions
# http://manuals.bioinformatics.ucr.edu/home/ht-seq

mydict <- DNAStringSet(sapply(1:10, function(x) paste(sample(c("A","T","G","C"), 8, replace=T), collapse=""))) # Creates random sample sequences.
names(mydict)<- paste("d", 1:10, sep="") # Names the sequences.

# see adapters.R for details on how to remove the adapters
experimentSeqs <- RNAStringSet(sread(fqTrim))
# filter the records down only to cel (C. Elegans)
# how to loop through RNAStringSet
for (seq in experimentSeqs) {
    print(seq)
    break
}

# get hairpin data
hairpinSeqs <- read.RNAStringSet(file=path.seqs, format="fasta")
# filter down to only c. elegans


# seq matching using vmatchPattern
for (seq in as.character(experimentSeqs)) {
    # get only > 18 long sequences
    mindex <- vmatchPattern(seq, hairpinSeqs, max.mismatch=2)
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

# store the c


