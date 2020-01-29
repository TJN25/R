#!/usr/bin/env Rscript
suppressMessages(library('getopt'))


spec = matrix(c(
  'sra', 's', 1, "character",
  'help' , 'h', 0, "logical",
  'reference' , 'r', 1, "character",
  'alignment' , 'a', 1, "character",
  'file_path', 'p', 2, "character",
  'out_name', 'o', 2, "character",
  'seq1', 'x', 2, "character",
  'seq2', 'y', 2, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat("sraPlotGenomeAlignment.R version 1.0\n")
  cat(" \n")
  cat("Use sraPlotGenomeAlignment.R <options> -s <sra plot file> -r <reference plot file> -a <alignment file>\n")
  cat(" \n")
  cat("Options:\n")
  cat("  -f <file path> The location of the other files and the output file\n")
  cat("  -s <sra plot file> The sra plot file to be reordered\n")
  cat("  -r <reference plot file> The sra file that the rearranged file needs to match\n")
  cat("  -a <alignment file> The file that will be used to rearrange the sra file (This needs to include the extension)\n")
  cat("  -x <seq1 column> The column number for seq 1\n")
  cat("  -y <seq2 column> The column number for seq 2\n")
  cat("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input\n")
  cat("  \n")
  q(status=1)
}

if ( is.null(opt$sra) ) {
  cat("Error: -s <sra plot file> is required.\n")
  q(status=1)
}
if ( is.null(opt$reference) ) {
  cat("Error: -r <reference file> is required.\n")
  
  q(status=1)
}

if ( is.null(opt$alignment) ) {
  cat("Error: -a <alignment file> is required.\n")
  
  q(status=1)
}

suppressMessages(library(tidyverse))
suppressMessages(library(tjnFunctions))


if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$out_name ) ) { opt$out_name = opt$sra }
if ( is.null(opt$x ) ) {  opt$x = "1" }
if ( is.null(opt$y ) ) { opt$y = "2" }


opt$x <- as.numeric(opt$x)
opt$y <- as.numeric(opt$y)


sraName <- opt$sra
gffName <- opt$gff
filePath <- opt$file_path



reference <- read.table(opt$reference)
sra <- read.table(opt$sra)
alignment <- read.table(opt$alignment, header = T)

sraReordered <- reorderSraPlot(sra = sra, alignment = alignment, reference = reference, seqA = opt$x, seqB = opt$y)


cat(paste("Writing the plot output to ", opt$file_path, "/", opt$out_name, "_reformat.plot\n", sep = ""))
write.table(x = sraReordered%>%select(a, b), file = paste(opt$file_path, "/", opt$out_name, "_reformat.plot", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
