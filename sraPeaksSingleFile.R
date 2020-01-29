#!/usr/bin/env Rscript
library(matrixStats)
library('getopt')


spec = matrix(c(
  'plotfile', 'f', 1, "character",
  'gcf', 'g', 1, "character",
  'help' , 'h', 0, "logical",
  'smooth' , 's', 0, "character",
  'file_path', 'p', 2, "character",
  'out_name', 'o', 2, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

if ( !is.null(opt$help) ) {
#print("NOT WORKING YET", quote = F)
print("sraPeaksSingleFile.R version 1.0", quote = F)
print(" ", quote = F)
print("Use sraPeaksSingleFile.R <options> -f <file list> -g <gca accession>", quote = F)
print(" ", quote = F)
print("Options:", quote = F)
print("  -f <file list> The list of the plot files to work with for calling regions", quote = F)
print("  -g <gca accession> The accession of the gff file that contains info to use in the new gff file", quote = F)
print("  -p <file path> The location of the other files and the output file", quote = F)
print("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input", quote = F)
print("  ", quote = F)
q(status=1)
}

# opt$plotfile <- "SRR1185095.plot"
# opt$gcf <- "GCA_000017745.1"
# opt$file_path <- "~/phd/RNASeq/escherichia/"
# opt$output <- "SRR1185095_test"

if ( is.null(opt$plotfile) ) {
  print("Error: -f <file list> is required.", quote = F)
  q(status=1)
}
#>> change this so that it is not needed or works better
if ( is.null(opt$gcf) ) {
  print("Error: -g <gca accession> is required.", quote = F)
  q(status=1)
}  
if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$out_name ) ) { opt$out_name = opt$plotfile }

library(tidyverse)
library(tjnFunctions)
  


gffMain <- readLines(paste(opt$file_path, "/", opt$gcf, ".gff", sep = ""))
gffMain <- data.frame(text = gffMain)
genomeInfo <- as.character(gffMain[8,1])
genomeBuild <- as.character(gffMain[4,1])
genomeSpecies <- as.character(gffMain[9,1])
accession <- strsplit(genomeInfo, " ")[[1]][2]

dat <- read.table(paste(opt$file_path, opt$plotfile, sep = "/"))

##include the removal of the cds regions to reduce a step
#>> removeCDSregions()

##Fwd ####
dfFwd <- sraSmoothandsort(dat = dat, col.num = 2, binwidth = 25)

dfFwd <- dfFwd%>%arrange(as.numeric(start))%>%filter(max.value > 5, length >=50)
dfShortFwd <- dfFwd%>%filter(length < 500, max.value < 1000)

##Rev ####
dfRev <- sraSmoothandsort(dat = dat, col.num = 1, binwidth = 25)

dfRev <- dfRev%>%arrange(as.numeric(start))%>%filter(max.value > 5, length >=50)
dfShortRev <- dfRev%>%filter(length < 500, max.value < 1000)

##Make gff file ####

gffFwd <- dfFwd%>%filter(max.value > 5)%>%mutate(strand = "+", 
                           source = "sraAlignedncRNAExpression",
                           seqname = accession,
                           strand = "+",
                           feature = "ncRNA",
                           frame = ".",
                           attribute = paste("ID=rna_fwd_", row_number(), sep = ""))%>%
  select(seqname, source, feature, start, end, median.val, strand, frame, attribute)

gffRev <- dfRev%>%filter(max.value > 5)%>%mutate(strand = "-", 
                           source = "sraAlignedncRNAExpression",
                           seqname = accession,
                           strand = "+",
                           feature = "ncRNA",
                           frame = ".",
                           attribute = paste("ID=rna_rev_", row_number(), sep = ""))%>%
  select(seqname, source, feature, start, end, median.val, strand, frame, attribute)%>%
  arrange(as.numeric(start))
gffRev <- gffRev%>%mutate(strand = rep("-", nrow(gffRev)))


gff <- gffFwd%>%bind_rows(gffRev)%>%arrange(as.numeric(start))


fileConn<-file(paste(opt$file_path, "/",opt$out_name, "_sra_calls_peaks.gff", sep = ""))
writeLines(c("##gff-version 3",
             "#!gff-spec-version 1.21",
             "#!processor R script (local) with manual add of top section",
             genomeBuild,
             paste("#!genome-build-accession NCBI_Assembly:", opt$gcf, sep = ""),
             paste("#!annotation-date ", Sys.Date(), sep = ""), 
             "#!annotation-source sraPlotSummary.R (local version)",
             genomeInfo,
             genomeSpecies), fileConn)
close(fileConn)



write.table(x = gff, file = paste(opt$file_path, "/",opt$out_name, "_sra_calls_peaks.gff", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)


