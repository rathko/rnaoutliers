
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

# instead of random use the actual sequences
# in addition make sure we remove the adapters before the alignment (see adapters.R)
fq <- readFastq(".", pattern = "B470001.1.fq")
sread(fq)

mypdict <- PDict(mydict) # Creates a PDict dictionary. Allows only sequences of same length. 

# get hairpin data
seqs <- read.RNAStringSet(file=path.seqs, format="fasta")
# match initial dictionary against every sequence from the hairpin (vmatchPDict is not implemented yet)
for (seq in as.character(seqs)) {
    mysearch <- matchPDict(mypdict, DNAString(RNAString(seq)), max.mismatch=0) # Searches all dictionary entries against chromosome. 
    #print(length(unlist(mysearch)))
    print(seq)
    len <- length(unlist(mysearch))
    if (len > 0) {
        #if (len > 5) len <- 5
        print(unlist(mysearch)[1:len])
    }
    #break;
}

## (2) Accessing the Results
unlist(mysearch)[1:10] # Returns first 10 matches.


# alignment counters

#path.mir <- "/Users/radek/reactor/phd/small_RNA/"
path.mir <- "/Users/radek/download/"
# cd ~/download; wget http://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz; gzip -d ./hairpin.fa.gz
path.seqs <- paste(path.mir,"hairpin.fa",sep="")
#seqs <- read.DNAStringSet(file=path.seqs, format="fastq")
seqs <- read.RNAStringSet(file=path.seqs, format="fasta")

# what should be seed.match?
seed.match <- "ACACUCC"
seq.ids <- names(seqs)
seq.lengths <- nchar(seqs)










