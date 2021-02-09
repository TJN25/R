#!/usr/bin/env Rscript
library('getopt')


spec = matrix(c(
  'cmscanOutput', 'f', 1, "character",
  'gcf', 'g', 1, "character",
  'help' , 'h', 0, "logical",
  'file_path', 'p', 2, "character",
  'out_name', 'o', 2, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)
#
# opt$cmscanOutput <- "GCA_000017745.1.tblout"
# opt$gcf <- "GCA_000017745.1"
# opt$file_path <- "~/phd/RNASeq/escherichia/"
# opt$output <- "escherichia_test"

if ( !is.null(opt$help) ) {
  cat("cmscanToGffWrapper.R version 1.0\n\n")
  cat("Use cmscanToGffWrapper.R <options> -f <cmscan ouptut file> -g <gff file>\n\n")
  cat("Options:\n")
  cat("  -f <cmscan ouptut file> The file that contains the cmscan output\n")
  cat("  -g <gff file> The file that contains the gff data. Do not inclue the gff file extension\n")
  cat("  -f <file path> The location of the other files and the output file\n")
  cat("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the gca input\n")
   q(status=1)
}

if ( is.null(opt$cmscanOutput) ) {
  cat("Error: -f <cmscan ouptut file> is required.\n")
  q(status=1)
}
if ( is.null(opt$gcf) ) {
  cat("Error: -g <gff file> is required.\n")
  q(status=1)
}

library(tidyverse)
library(comparativeSRA)

if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$output ) ) { opt$output = opt$gcf }

rfamRes <- read.table(paste(opt$file_path, opt$cmscanOutput, sep = "/"), header = F, comment.char = "#",quote = "", fill = T)


gff <- cmscanToGff(rfamRes = rfamRes)


gffMain <- readLines(paste(opt$file_path, "/", opt$gcf, ".gff", sep = ""))
gffMain <- data.frame(text = gffMain)
genomeInfo <- as.character(gffMain[8,1])
genomeBuild <- as.character(gffMain[4,1])
genomeSpecies <- as.character(gffMain[9,1])
accession <- strsplit(genomeInfo, " ")[[1]][2]

fileConn<-file(paste(opt$file_path, "/",opt$output, "_ncRNA.gff", sep = ""))
writeLines(c("##gff-version 3",
             "#!gff-spec-version 1.21",
             "#!processor R script (local)",
             genomeBuild,
             paste("#!genome-build-accession NCBI_Assembly:", opt$gcf, sep = ""),
             paste("#!annotation-date ", Sys.Date(), sep = ""),
             "#!annotation-source cmscan (rFam) (local version)",
             genomeInfo,
             genomeSpecies), fileConn)
close(fileConn)


write.table(x = gff, file = paste(opt$file_path, "/",opt$output, "_ncRNA.gff", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)
