
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



# notes on adapters removal
fq <- readFastq(".", pattern = "B470001.1.fq")
detail(fq)

# checking the quality
# estimate the frequency of bases per cycle
abc <- alphabetByCycle(sread(fq), alphabet=c("A","T","G","C","N"))
abcFreq <- abc / length(fq)

cols <- c("green","blue","yellow","red","grey")
plot(NULL, xlim=c(1,37), ylim=c(0,1), xlab="Cycle", ylab="Frequency", main="Base frequency by cycle")
matlines(t(abcFreq), col=cols, lty=1)
legend("topleft",legend=rownames(abcFreq), lty=1, col=cols, bty="n", horiz=TRUE)

# look at the Phred quality scores
numQuals <- as(quality(fq), "matrix")
head(numQuals)
# see this as a figure
meanQ <- colMeans(numQuals)
minQ  <- apply(numQuals,2,min)
maxQ  <- apply(numQuals,2,max)
plot(meanQ, xlim=c(1,37), ylim=c(1,40), xlab="Cycle", 
  	ylab="Phred Quality", main="Max/Mean/Min quality by cycle", col="blue")
points(minQ, col="red")
points(maxQ, col="green")

# tables
tables(fq)

# remoe poly A, if any
polyA   <- "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
polyADist   <- srdistance(sread(fq), polyA)[[1]]
hist(polyADist, xlim=c(1,37), main="Histogram of distances to polyA", col="skyblue")

fqClean <- fq[polyADist >= 25]
length(fqClean)

# Filter N's, these bases are not relevant when mapping the sequence anyway
filterN <- nFilter(threshold=3)
sread(fqClean[!filterN(sread(fqClean))])

fqClean <- fqClean[filterN(sread(fqClean))]

# Removing the adapter
adapter <- "ATGGAATTCTCGGGTGCCA"
#adapter <- "TCGTATGCCGTCTTCTGCTTGTT"
sequences <- sread(fqClean)
quality <- quality(quality(fqClean))

adapterPad <- paste(substr(adapter,1,12), paste(rep("N", 25),collapse=""), sep="")
adapterPad

maxMm <- round(seq_len(12)*0.2)
maxMmPad <-  c(maxMm,rep(maxMm[length(maxMm)], 25))
maxMmPad

coordinates <- trimLRPatterns(Rpattern=adapterPad, subject=sequences, max.Rmismatch=maxMmPad, 
		ranges=TRUE, with.Rindels=FALSE, Rfixed=FALSE)
seqTrim <- DNAStringSet(sequences, start=start(coordinates), end=end(coordinates))
qualTrim <- BStringSet(quality, start=start(coordinates), end=end(coordinates))

qualTrim <- FastqQuality(qualTrim)
fqTrim <- ShortReadQ(sread=seqTrim, quality=qualTrim, id=id(fqClean))

# check trimmed sequences
tables(fqTrim)$top

# TODO: compare this to the sRNA workbench results

# in the next step run these analysis on all files
####### Actual analysis
