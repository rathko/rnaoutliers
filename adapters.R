# Fasta format: http://www.biostars.org/p/3398/
# For a good tutorial see: http://rcourse.iop.kcl.ac.uk/2011/4thu/session4/Course_notes.html

#setwd("C:/devel/data/small_RNA")
#setwd("/Users/radek/reactor/phd/small_RNA")
folder.project <- "/Users/radek/reactor/rnaoutliers/"
folder.raw <- "/Users/radek/reactor/phd/small_RNA/"

require(ShortRead)
source(paste(folder.project, "util.R", sep=""))


# just for testing
files <- c("B470004.1.fq", "B470005.1.fq", "B470006.1.fq")

files <- list.files(folder.raw)
for(file in files) {
    print(paste("Reading file: ", file))
	fq <- readFastq(folder.raw, pattern = file)
	fqClean <- removeAdapters(fq)
	# store the cleaned files
	fileClean <- paste(folder.raw, "clean/", file, sep="")
	writeFastq(fqClean, fileClean)
}

