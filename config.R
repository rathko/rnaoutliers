
folder.project <- "/Users/radek/reactor/rnaoutliers"

library(Biostrings)
library(ShortRead)
library(preprocessCore)
library(gdata)

hairpin.path <- "/Users/radek/download/hairpin.fa"
folder.clean <- "/Users/radek/reactor/phd/small_RNA/clean/"
folder.counters <- "/Users/radek/reactor/phd/small_RNA/counters/"

# all common configuration
samples <- read.csv(paste(folder.project, "conditions.csv", sep="/"), sep=",", 
	strip.white = TRUE)
