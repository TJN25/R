#!/usr/bin/env Rscript
library(matrixStats)
library('getopt')


spec = matrix(c(
  'list', 'l', 1, "character",
  'gcf', 'g', 1, "character",
  'help' , 'h', 0, "logical",
  'smooth' , 's', 0, "character",
  'file_path', 'p', 2, "character",
  'out_name', 'o', 2, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

if ( !is.null(opt$help) ) {
  #print("NOT WORKING YET", quote = F)
  print("sraPlotSummary.R version 1.0", quote = F)
  print(" ", quote = F)
  print("Use sraPlotSummary.R <options> -l <file list> -g <gca accession>", quote = F)
  print(" ", quote = F)
  print("Options:", quote = F)
  print("  -l <file list> The list of the plot files to work with for calling regions", quote = F)
  print("  -g <gca accession> The accession of the gff file that contains info to use in the new gff file", quote = F)
  print("  -p <file path> The location of the other files and the output file", quote = F)
  print("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input", quote = F)
  print("  ", quote = F)
  q(status=1)
}

opt$list <- "serratia_ncRNA_sra_plot_list.txt"
opt$gcf <- "GCA_000438825.1"
opt$file_path <- "~/phd/RNASeq/serratia/"
opt$output <- "serratia_test"

if ( is.null(opt$list) ) {
  print("Error: -l <file list> is required.", quote = F)
  q(status=1)
}
#>> change this so that it is not needed or works better
if ( is.null(opt$gcf) ) {
  print("Error: -g <gca accession> is required.", quote = F)
  q(status=1)
}  
if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$out_name ) ) { opt$out_name = opt$list }

library(tidyverse)
library(tjnFunctions)




fileList <- readLines(paste(opt$file_path, opt$list, sep = "/"))
fileList <- fileList[1]
#fileList <- readLines("~/phd/RNASeq/serratia/serratia_ncRNA_sra_plot_list.txt")
#fileList <- readLines("~/phd/RNASeq/xanthomonas/xanthomonas_ncRNA_sra_plot_list.txt")
#gffMain <- readLines("~/phd/RNASeq/serratia/GCA_000438825.1.gff")
gffMain <- readLines(paste(opt$file_path, "/", opt$gcf, ".gff", sep = ""))
gffMain <- data.frame(text = gffMain)
genomeInfo <- as.character(gffMain[8,1])
genomeBuild <- as.character(gffMain[4,1])
genomeSpecies <- as.character(gffMain[9,1])
accession <- strsplit(genomeInfo, " ")[[1]][2]

dat <- read.table(paste(opt$file_path, fileList[1], sep = "/"))
datMatFwd <- matrix(0, nrow = nrow(dat), ncol = length(fileList))
datMatRev <- matrix(0, nrow = nrow(dat), ncol = length(fileList))

for(i in 1:length(fileList)){
  print(i)
  tmp <- read.table(paste(opt$file_path, fileList[i], sep = "/"))
  datMatFwd[,i] <- tmp[,1]
  datMatRev[,i] <- tmp[,2]
}

datFwd <- as.data.frame(datMatFwd)
datFwd$min.val <- rowMins(as.matrix(datFwd))
datFwd$median.val <- rowMedians(as.matrix(datFwd))
datFwd$max.val <- rowMaxs(as.matrix(datFwd))

datRev <- as.data.frame(datMatRev)
datRev$min.val <- rowMins(as.matrix(datRev))
datRev$median.val <- rowMedians(as.matrix(datRev))
datRev$max.val <- rowMaxs(as.matrix(datRev))


minDat <- as.data.frame(cbind(datFwd[,length(fileList) + 1], datRev[,length(fileList) + 1]))
medianDat <- as.data.frame(cbind(datFwd[,length(fileList) + 2], datRev[,length(fileList) + 2]))
maxDat <- as.data.frame(cbind(datFwd[,length(fileList) + 3], datRev[,length(fileList) + 3]))




##Fwd ####
minDfFwd <- sraSmoothandsort(dat = minDat, binwidth = 50)
medianDfFwd <- sraSmoothandsort(dat = medianDat, binwidth = 50)
maxDfFwd <- sraSmoothandsort(dat = maxDat, binwidth = 50)


minDfFwd <- minDfFwd%>%mutate(group = "Min")
medianDfFwd <- medianDfFwd%>%mutate(group = "Median")%>%arrange(as.numeric(start))
maxDfFwd <- maxDfFwd%>%mutate(group = "Max")

dfFwd <- minDfFwd%>%bind_rows(medianDfFwd, maxDfFwd)%>%mutate(group = as.factor(group))%>%filter(max.value > 5)
dfFwd <- minDfFwd%>%bind_rows(medianDfFwd, maxDfFwd)%>%mutate(group = as.factor(group))%>%filter(max.value > 5)

dfShortFwd <- dfFwd%>%filter(length < 500, max.value < 1000)

##Rev ####
minDfRev <- sraSmoothandsort(dat = minDat, col.num = 2, binwidth = 50)
medianDfRev <- sraSmoothandsort(dat = medianDat, col.num = 2, binwidth = 50)
maxDfRev <- sraSmoothandsort(dat = maxDat, col.num = 2, binwidth = 50)


minDfRev <- minDfRev%>%mutate(group = "Min")
medianDfRev <- medianDfRev%>%mutate(group = "Median")%>%arrange(as.numeric(start))
maxDfRev <- maxDfRev%>%mutate(group = "Max")

dfRev <- minDfRev%>%bind_rows(medianDfRev, maxDfRev)%>%mutate(group = as.factor(group))%>%filter(max.value > 5)
dfShortRev <- dfRev%>%filter(length < 500, max.value < 1000)
# 
# ##lengths ####
# ggplot(data = dfFwd, aes(x=length, y = ..density.., fill = group)) +
#   geom_histogram(colour="black", binwidth = 25) +
#   facet_grid(group ~ .) +
#   theme(legend.position = "none")
# 
# ggplot(data = dfShortFwd, aes(x=length, y = ..density.., fill = group)) +
#   geom_histogram(colour="black", binwidth = 25) +
#   facet_grid(group ~ .) +
#   theme(legend.position = "none")
# 
# ##median ####
# ggplot(data = dfFwd, aes(x=median.val, y = ..density.., fill = group)) +
#   geom_histogram(colour="black", binwidth = 25) +
#   facet_grid(group ~ .) +
#   theme(legend.position = "none")
# 
# ggplot(data = dfShortFwd, aes(x=median.val, y = ..density.., fill = group)) +
#   geom_histogram(colour="black", binwidth = 25) +
#   facet_grid(group ~ .) +
#   theme(legend.position = "none")
# 
# ##max ####
# ggplot(data = dfFwd, aes(x=max.value, y = ..density.., fill = group)) +
#   geom_histogram(colour="black", binwidth = 25) +
#   facet_grid(group ~ .) +
#   theme(legend.position = "none")
# 
# ggplot(data = dfShortFwd, aes(x=max.value, y = ..density.., fill = group)) +
#   geom_histogram(colour="black", binwidth = 25) +
#   facet_grid(group ~ .) +
#   theme(legend.position = "none")
# 
# 
##Make gff file ####

gffFwd <- medianDfFwd%>%filter(max.value > 5)%>%mutate(strand = "+", 
                                                       source = "sraAlignedncRNAExpression",
                                                       seqname = accession,
                                                       strand = "+",
                                                       feature = "ncRNA",
                                                       frame = ".",
                                                       attribute = "NA")%>%
  select(seqname, source, feature, start, end, median.val, strand, frame, attribute)

gffRev <- medianDfRev%>%filter(max.value > 5)%>%mutate(strand = "-", 
                                                       source = "sraAlignedncRNAExpression",
                                                       seqname = accession,
                                                       strand = "+",
                                                       feature = "ncRNA",
                                                       frame = ".",
                                                       attribute = "NA")%>%
  select(seqname, source, feature, start, end, median.val, strand, frame, attribute)%>%
  arrange(as.numeric(start))
gffRev <- gffRev%>%mutate(strand = rep("-", nrow(gffRev)))


gff <- gffFwd%>%bind_rows(gffRev)%>%arrange(as.numeric(start))


fileConn<-file("~/phd/RNASeq/serratia/GCA_000438825.1_ncRNA_plotCalls.gff")
writeLines(c("##gff-version 3",
             "#!gff-spec-version 1.21",
             "#!processor R script (local) with manual add of top section",
             genomeBuild,
             paste("#!genome-build-accession NCBI_Assembly:", opt$gcf, sep = ""),
             paste("#!annotation-date ", Sys.Date(), sep = T), 
             "#!annotation-source sraPlotSummary.R (local version)",
             genomeInfo,
             genomeSpecies), fileConn)
close(fileConn)



write.table(x = gff, file = "~/phd/RNASeq/serratia/GCA_000438825.1_ncRNA_plotCalls.gff", row.names = F, col.names = F, quote = F, sep = "\t", append = T)


#####
#write.table(x = minDat, file = "~/phd/RNASeq/serratia/serratia_min_ncRNA.plot", quote = F, row.names = F, col.names = F, sep = "\t")
#write.table(x = medianDat, file = "~/phd/RNASeq/serratia/serratia_median_ncRNA.plot", quote = F, row.names = F, col.names = F, sep = "\t")
#write.table(x = maxDat, file = "~/phd/RNASeq/serratia/serratia_max_ncRNA.plot", quote = F, row.names = F, col.names = F, sep = "\t")

