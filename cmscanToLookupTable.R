#!/usr/bin/env Rscript
library('getopt')

cmscanToLookup <- function(rfamRes){
  colnames(rfamRes) <- c("idx", "target.name", "accession1", "query.name", "accession2", "clan.name", "mdl", "mdl.from",  "mdl.to", "seq.from",   "seq.to",
                         "strand", "trunc", "pass",   "gc",  "bias",  "score",   "E.value", "inc", "olp", "anyidx", "afrct1", "afrct2", "winidx", "wfrct1", "wfrct2", "description.of.target")
  
  rfamRes <- rfamRes%>%
    mutate(seq.from2 = seq.from)%>%
    mutate(seq.from = ifelse(seq.from > seq.to, seq.to, seq.from))%>%
    mutate(seq.to = ifelse(seq.from2 > seq.to, seq.from2, seq.to))
  
  gff <- data.frame(seqname = rfamRes$query.name,
                    source = rep("rfam", nrow(rfamRes)),
                    feature = rep("ncRNA", nrow(rfamRes)),
                    start = rfamRes$seq.from,
                    end = rfamRes$seq.to,
                    score = rfamRes$score,
                    strand = rfamRes$strand,
                    frame = rep(".", nrow(rfamRes)),
                    attribute = rfamRes$accession1)
  return(gff)
}


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
library(tjnFunctions)

if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$output ) ) { opt$output = opt$gcf }

rfamRes <- read.table(paste(opt$file_path, opt$cmscanOutput, sep = "/"), header = F, comment.char = "#",quote = "", fill = T)


gff <- cmscanToLookup(rfamRes = rfamRes)
gff <- gff %>% arrange(start)


write.table(x = gff, file = paste(opt$file_path, "/",opt$output, "_ncRNA.lookup", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")