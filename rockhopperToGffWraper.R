#!/usr/bin/env Rscript
library('getopt')

##Set up####
##labelled as cmscan however the input is from rockhopper and will require different formatting
spec = matrix(c(
  'cmscanOutput', 'f', 1, "character",
  'gcf', 'g', 1, "character",
  'help' , 'h', 0, "logical",
  'file_path', 'p', 2, "character",
  'out_name', 'o', 2, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)
# 
 opt$cmscanOutput <- "SRR1185095_NC_009801_transcripts.txt"
 opt$gcf <- "GCA_000017745.1"
 opt$file_path <- "~/phd/RNASeq/escherichia/"
 opt$output <- "SRR1185095"

##Help options####
if ( !is.null(opt$help) ) {
  print("NOT CORRECT")
  # print("removeProteinCoding.R version 1.0", quote = F)
  # print(" ", quote = F)
  # print("Use removeProteinCoding.R <options> -s <sra plot file> -g <gff file>", quote = F)
  # print(" ", quote = F)
  # print("Options:", quote = F)
  # print("  -s <sra plot file> The file that contains the plot data. Do not inclue the .plot file extension", quote = F)
  # print("  -g <gff file> The file that contains the gff data. Do not inclue the gff file extension", quote = F)
  # print("  -f <file path> The location of the other files and the output file", quote = F)
  # print("  -r <protein coding range> The number of nucleotides either side of a CDS region that should also be set to zero", quote = F)
  # print("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input", quote = F)
  q(status=1)
}
if ( is.null(opt$cmscanOutput) ) {
  print("NOT CORRECT")
  # print("Error: -s <sra plot file> is required.", quote = F)
  # print(" ", quote = F)
  # print("removeProteinCoding.R version 1.0", quote = F)
  # print("Use removeProteinCoding.R <options> -s <sra plot file> -g <gff file>", quote = F)
  # print(" ", quote = F)
  # print("Options:", quote = F)
  # print("  -s <sra plot file> The file that contains the plot data. Do not inclue the .plot file extension", quote = F)
  # print("  -g <gff file> The file that contains the gff data. Do not inclue the gff file extension", quote = F)
  # print("  -f <file path> The location of the other files and the output file", quote = F)
  # print("  -r <protein coding range> The number of nucleotides either side of a CDS region that should also be set to zero", quote = F)
  # print("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input", quote = F)
  q(status=1)
}
if ( is.null(opt$gcf) ) {
  print("NOT CORRECT")
  # print("Error: -s <sra plot file> is required.", quote = F)
  # print(" ", quote = F)
  # print("removeProteinCoding.R version 1.0", quote = F)
  # print("Use removeProteinCoding.R <options> -s <sra plot file> -g <gff file>", quote = F)
  # print(" ", quote = F)
  # print("Options:", quote = F)
  # print("  -s <sra plot file> The file that contains the plot data. Do not inclue the .plot file extension", quote = F)
  # print("  -g <gff file> The file that contains the gff data. Do not inclue the gff file extension", quote = F)
  # print("  -f <file path> The location of the other files and the output file", quote = F)
  # print("  -r <protein coding range> The number of nucleotides either side of a CDS region that should also be set to zero", quote = F)
  # print("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input", quote = F)
  q(status=1)
}

library(tidyverse)
library(tjnFunctions)

if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$output ) ) { opt$output = opt$gcf }
#####
 
 
 rockhopperToGff <- function(rockhopperRes){
   
   
   rhProteins <- rockhopperRes%>%filter(!is.na(Translation.Start))
   rhncRNA <- rockhopperRes%>%filter(is.na(Translation.Start))
   
   rhProteins <- rhProteins%>%mutate(start = ifelse(Translation.Start < Translation.Stop, Translation.Start,Translation.Stop),
                                     stop = ifelse(Translation.Start < Translation.Stop, Translation.Stop,Translation.Start))%>%select(-Transcription.Start, 
                                                                                                  -Translation.Start, 
                                                                                                  -Translation.Stop,
                                                                                                  -Transcription.Stop)%>%
     mutate(feature.type = "CDS")
   
   rhncRNA <- rhncRNA%>%mutate(start = ifelse(Transcription.Start < Transcription.Stop, Transcription.Start,Transcription.Stop),
                               stop = ifelse(Transcription.Start < Transcription.Stop, Transcription.Stop,Transcription.Start))%>%select(-Transcription.Start, 
                                                                                                -Translation.Start, 
                                                                                                -Translation.Stop,
                                                                                                -Transcription.Stop)%>%
     mutate(feature.type = "ncRNA")
   
   
   rockhopper <- rhProteins%>%bind_rows(rhncRNA)
   
   gff <- data.frame(seqname = rep(accession, nrow(rockhopper)),
                     source = rep("rockhopper", nrow(rockhopper)),
                     feature = rockhopper$feature.type,
                     start = rockhopper$start,
                     end = rockhopper$stop,
                     score = rockhopper$Expression.1,
                     strand = rockhopper$Strand,
                     frame = rep(".", nrow(rockhopper)),
                     attribute = paste("ID=",rockhopper$Name, "_", row_number(), sep = ""))
   

   return(gff)
 }
 
 ##import exsisting gff file in order to get required information
 gffMain <- readLines(paste(opt$file_path, "/", opt$gcf, ".gff", sep = ""))
 gffMain <- data.frame(text = gffMain)
 genomeInfo <- as.character(gffMain[8,1])
 genomeBuild <- as.character(gffMain[4,1])
 genomeSpecies <- as.character(gffMain[9,1])
 accession <- strsplit(genomeInfo, " ")[[1]][2]
 
 
 ##read in the rockhopper output
rockhopperRes <- read.table(paste(opt$file_path, opt$cmscanOutput, sep = "/"), header = F, comment.char = "#",quote = "", fill = T, sep = "\t")
colnames(rockhopperRes) <- c("Transcription.Start", "Translation.Start", "Translation.Stop", "Transcription.Stop", "Strand", "Name", "Synonym", "Product",	"Expression.1")





gff <- rockhopperToGff(rockhopperRes = rockhopperRes)




fileConn<-file(paste(opt$file_path, "/",opt$output, "_rockhopper.gff", sep = ""))
writeLines(c("##gff-version 3",
             "#!gff-spec-version 1.21",
             "#!processor R script (local)",
             genomeBuild,
             paste("#!genome-build-accession NCBI_Assembly:", opt$gcf, sep = ""),
             paste("#!annotation-date ", Sys.Date(), sep = ""), 
             "#!annotation-source RockHopper (local version)",
             genomeInfo,
             genomeSpecies), fileConn)
close(fileConn)


write.table(x = gff, file = paste(opt$file_path, "/",opt$output, "_rockhopper.gff", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)
