# Fasta format: http://www.biostars.org/p/3398/
# For a good tutorial see: http://rcourse.iop.kcl.ac.uk/2011/4thu/session4/Course_notes.html

#setwd("C:/devel/data/small_RNA")
#setwd("/Users/radek/reactor/phd/small_RNA")
require(ShortRead)

files <- list.files()
for(file in files) {
    print(paste("Reading file: ", file))
	fq <- readFastq(".", pattern = file)
	fqClean <- removeAdapters(fq)
	# store the cleaned files
	fileClean <- paste("clean/", file, sep="")
	writeFastq(fqClean, fileClean)
}



