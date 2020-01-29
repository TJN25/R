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
if ( is.null(opt$gff) ) {
  cat("Error: -g <gff file> is required.\n")
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


gffDat  <- tryCatch({
  suppressWarnings(gffDat <- read.table(paste(filePath, "/", gffName, ".gff", sep = ""), sep = "\t", fill = T, comment.char = "#", quote = ""))
  gffDat
}, error =  function(e) {
  cat(paste("Error: ", opt$file_path, "/", opt$gff, ".gff not found.\n", sep = ""))
  q(status=1)
})


plotDat  <- tryCatch({
  suppressWarnings(plotDat <- read.table(paste(filePath, "/", sraName, ".plot", sep = "")))
  plotDat
}, error =  function(e) {
  cat(paste("Error: ", opt$file_path, "/", opt$sra, ".plot not found.\n", sep = ""))
  q(status=1)
})

total <- (sum(plotDat$V1) + sum(plotDat$V2))/1000000

colnames(gffDat) <- c("sequence", "source", "feature", "start", "end", "score", "strand", "phase", "Atrribute")

plotDatncRNA  <- tryCatch({
  ##Change this path or put the header file in the working directory
  suppressWarnings(plotDatncRNA <- read.table(paste(filePath, "/", sraName, "_ncRNA.plot", sep = "")))

  plotDatncRNA
}, error =  function(e) {
  plotDat <- read.table(paste(filePath, "/", sraName, ".plot", sep = ""))
  if(quiet == F){
    cat("Running removeCDSregions\n")

  }
  plotDatncRNA <- removeCDSregions(plotDat = plotDat, gffDat = gffDat, stranded = stranded, time.it = T)

  plotDatncRNA
} )

plotDatncRNA$V1 <- plotDatncRNA$V1/total
plotDatncRNA$V2 <- plotDatncRNA$V2/total



cat("Running rnPeakCalling\n")

cat("Calling forward\n")
callsDatFwd <- rnaPeakCalling(dat = plotDatncRNA, col.num = 2, small_peaks = F, plot_threshold = 15/total)

cat("Calling reverse\n")
callsDatRev <- rnaPeakCalling(dat = plotDatncRNA, col.num = 1, small_peaks = F, plot_threshold = 15/total)


callsDatFwd <- callsDatFwd%>%mutate(strand = "+")
callsDatRev <- callsDatRev%>%mutate(strand = "-")

runningTime <- proc.time() - ptm
  printRunningTime(runningTime = runningTime)

  if("feature.length" %in% colnames(callsDatFwd) == F){
    print(colnames(callsDatFwd))
    print(head(callsDatFwd))
    cat("Warning: feature.length column not found in callsDatFwd.\n")
    quitStatus <- T
    callsDatFwdTmp <- callsDatFwd%>%filter(start != 0)%>%
      mutate(feature.length = stop - start)%>%
      mutate(feature.score = feature.length*mean.score)%>%
      filter(feature.score > 3)
  }else{
  callsDatFwdTmp <- callsDatFwd%>%filter(start != 0)%>%
    mutate(feature.score = feature.length*mean.score)%>%
    filter(feature.score > 3)
  }
  if("feature.length" %in% colnames(callsDatRev) == F){
    print(colnames(callsDatRev))
    print(head(callsDatRev))
    cat("Warning: feature.length column not found in callsDatRevTmp.\n")
    quitStatus <- T
    callsDatRevTmp <- callsDatRev%>%filter(start != 0)%>%
      mutate(feature.length = stop - start)%>%
      mutate(feature.score = feature.length*mean.score)%>%
      filter(feature.score > 3)
  }else{
  callsDatRevTmp <- callsDatRev%>%filter(start != 0)%>%
    mutate(feature.score = feature.length*mean.score)%>%
    filter(feature.score > 3)
  }
  
  # if(quitStatus == T){
  #   q(status=1)
  # }

gffMain <- readLines(paste(filePath, "/", gffName, ".gff", sep = ""))
gffMain <- data.frame(text = gffMain)
genomeInfo <- as.character(gffMain[8,1])
genomeBuild <- as.character(gffMain[4,1])
genomeSpecies <- as.character(gffMain[9,1])
accession <- strsplit(genomeInfo, " ")[[1]][2]



gffFwd <- callsDatFwdTmp%>%mutate(strand = "+",
                                                 source = "sraAlignedncRNAExpression",
                                                 seqname = accession,
                                                 median.val = round(mean.score*100),
                                                 feature = "ncRNA",
                                                 frame = ".",
                                                 attribute = paste("ID=rna_fwd_", row_number(), sep = ""))%>%
  select(seqname, source, feature, start, stop, median.val, strand, frame, attribute)

gffRev <- callsDatRevTmp%>%mutate(strand = "-",
                                                 source = "sraAlignedncRNAExpression",
                                                 seqname = accession,
                                                 median.val = round(mean.score*100),
                                                 feature = "ncRNA",
                                                 frame = ".",
                                                 attribute = paste("ID=rna_rev_", row_number(), sep = ""))%>%
  select(seqname, source, feature, start, stop, median.val, strand, frame, attribute)%>%
  arrange(as.numeric(start))



gff <- gffFwd%>%bind_rows(gffRev)%>%arrange(as.numeric(start))
gff <- gff%>%filter(start != 0)




fileConn<-file(paste(filePath, "/", opt$out_name, "_sra_calls.gff", sep = ""))
writeLines(c("##gff-version 3",
             "#!gff-spec-version 1.21",
             "#!processor R script (local) with manual add of top section",
             genomeBuild,
             paste("#!genome-build-accession NCBI_Assembly:", opt$gff, sep = ""),
             paste("#!annotation-date ", Sys.Date(), sep = ""),
             "#!annotation-source sraPlotSummary.R (local version)",
             genomeInfo,
             genomeSpecies), fileConn)
close(fileConn)

cat(paste("Writing the gff output to ", filePath, "/", opt$out_name, "_sra_calls.gff\n", sep = ""))
write.table(x = gff, file = paste(filePath, "/", opt$out_name, "_sra_calls.gff", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)


