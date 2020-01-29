#!/usr/bin/env Rscript
suppressMessages(library('getopt'))

spec = matrix(c(
  'sra', 'f', 1, "character",
  'help' , 'h', 0, "logical",
  'stranded' , 's', 0, "logical",
  'quiet' , 'q', 0, "logical",
  'gff' , 'g', 1, "character",
  'file_path', 'p', 2, "character",
  'range', 'r', 2, "integer",
  'out_name', 'o', 2, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat("run_rnaPeakCalling.R version 1.0\n")
  cat(" \n")
  cat("Use run_rnaPeakCalling.R <options> -f <sra plot file> -g <gff file>\n")
  cat(" \n")
  cat("Options:\n")
  cat("  -f <sra plot file> The file that contains the plot data. Do not inclue the .plot file extension\n")
  cat("  -g <gff file> The file that contains the gff data. Do not inclue the gff file extension\n")
  cat("  -s <stranded data> The data is stranded\n")
  cat("  -q <quiet> Do not print any updates\n")
  cat("  -p <file path> The location of the other files and the output file\n")
  cat("  -r <protein coding range> The number of nucleotides either side of a CDS region that should also be set to zero\n")
  cat("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input\n")
  q(status=1)
}

if ( is.null(opt$sra) ) {
  cat("Error: -f <sra plot file> is required.\n")
  q(status=1)
}




suppressMessages(library(tidyverse))
suppressMessages(library(tjnFunctions))
###--- column 1 is reverse and column 2 is forward ---###

if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$range ) ) { opt$range = 50 }
if ( is.null(opt$out_name ) ) { opt$out_name = opt$sra }
if(is.null(opt$stranded)){
  stranded <- F
}else{
  stranded <- T
}

if(is.null(opt$quiet)){
  quiet <- F
}else{
  quiet <- T
}

sraName <- opt$sra
gffName <- opt$gff
filePath <- opt$file_path

ptm <- proc.time()




plotDat  <- tryCatch({
  suppressWarnings(plotDat <- read.table(paste(filePath, "/", sraName, ".plot", sep = "")))
  plotDat
}, error =  function(e) {
  cat(paste("Error: ", opt$file_path, "/", opt$sra, ".plot not found.\n", sep = ""))
  q(status=1)
})

total <- (sum(plotDat$V1) + sum(plotDat$V2))/1000000
print(15/total)
print(total)
