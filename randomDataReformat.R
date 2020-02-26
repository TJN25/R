#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library('getopt'))

spec = matrix(c(
  'help' , 'h', 0, "logical",
  'quiet' , 'q', 0, "logical",
  'gff' , 'g', 1, "character",
  'input' , 'i', 1, "character",
  'file_path', 'p', 2, "character",
  'out_name', 'o', 2, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(" version 1.0\n")
  cat(" \n")
  cat("Use  <options> -f <sra plot file> -g <gff file>\n")
  cat(" \n")
  cat("Options:\n")
  cat("  -g <gff file> The file that contains the gff data. Do not inclue the gff file extension\n")
  cat("  -p <file path> The location of the other files and the output file\n")
  cat("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input\n")
  q(status=1)
}


if ( is.null(opt$gff) ) {
  cat("Error: -g <gff file> is required.\n")
  q(status=1)
}



suppressMessages(library(tidyverse))
suppressMessages(library(tjnFunctions))
###--- column 1 is reverse and column 2 is forward ---###

if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$out_name ) ) { opt$out_name = opt$input }

if ( is.null(opt$input ) ) { opt$input = opt$gff }



gffName <- opt$gff
filePath <- opt$file_path








##change this
sraDat  <- tryCatch({
suppressWarnings(sraDat <- read.table(paste(filePath, "/", opt$input, "_random_data.txt", sep = ""), sep = "\t", header = T))
  sraDat
}, error =  function(e) {
  cat(paste("Error: ", opt$file_path, "/", opt$input, sep = ""))
  q(status=1)
})


gffMain <- readLines(paste(filePath, "/gff_files/", gffName, ".gff", sep = ""))
gffMain <- data.frame(text = gffMain)
genomeInfo <- as.character(gffMain[8,1])
genomeBuild <- as.character(gffMain[4,1])
genomeSpecies <- as.character(gffMain[9,1])
accession <- strsplit(genomeInfo, " ")[[1]][2]

#Change this
gff <- sraDat%>%mutate(strand = strandR,
                       start = startR,
                       stop = endR,
                       source = "random_sequence",
                       seqname = accession,
                       median.val = 1,
                       feature = type,
                       frame = ".",
                       attribute = paste("ID=random_rna_", row_number(), sep = ""))%>%
  select(seqname, source, feature, start, stop, median.val, strand, frame, attribute) %>%
  arrange(as.numeric(start))



gff <- gff%>%filter(start != 0)




fileConn<-file(paste(filePath, "/", opt$out_name, "_random_sra.gff", sep = ""))
writeLines(c("##gff-version 3",
             "#!gff-spec-version 1.21",
             "#!processor R script (local) with manual add of top section",
             genomeBuild,
             paste("#!genome-build-accession NCBI_Assembly:", opt$gff, sep = ""),
             paste("#!annotation-date ", Sys.Date(), sep = ""),
             "#!annotation-source randomDataReformat.R (local version)",
             genomeInfo,
             genomeSpecies), fileConn)
close(fileConn)

cat(paste("Writing the gff output to ", filePath, "/", opt$out_name, "_random_sra.gff\n", sep = ""))
write.table(x = gff, file = paste(filePath, "/", opt$out_name, "_random_sra.gff", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)





