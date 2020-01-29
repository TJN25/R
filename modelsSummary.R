#!/usr/bin/env Rscript
##Initial Inputs ####
options(warn = -1)

#!/usr/bin/env Rscript
##Initial Inputs ####
options(warn = -1)

# args <- c(
#   "--bacteria","Bdellovibrio_out.txt",
#   "--phage", "phage_info_19-06-2018.txt",
#   "--genus", "Bdellovibrio",
#   "--out_name", "Bdellovibrio_out_test",
#   "--aa", "Bdellovibrio_out.amino_acids.txt",
#   "--wd", "/Volumes/2TTJN/Virus_Host_Prediction/",
#   "--trna", "Bdellovibrio_out.tRNA.txt",
#   "--cds", "Bdellovibrio_out.CDS.txt",
#   "--out_type", "plot")

#args <- c("--genus", "Lactococcus")
#args <- c("--genus", "Clostridium")

args <- commandArgs(trailingOnly = T)

gc_plotting = T
order_aa = T

if(length(args) < 1){
  print("   Required fields:", quote = F)
  print("      --input (bacterial info file),", quote = F )
  print(" ", quote = F)
  print("   Other fields:", quote = F)
  print("      --out_name (output file prefix),", quote = F)
  quit()
}
i <- 11
#print(args)
for(i in 1:length(args)){
  input_type <- args[i]
  if(input_type == "--help"){
    print("   Required fields:", quote = F)
    print("      --input (bacterial info file),", quote = F )
    print(" ", quote = F)
    print("   Other fields:", quote = F)
    print("      --out_name (output file prefix),", quote = F)
    quit()

  }else if(input_type == "--input"){
    input_name <- args[i + 1]
    print(args[i + 1])
  }else if(input_type == "--out_name"){
    input_out_name <- args[i + 1]
}

}

##Required


if(exists("input_name") == F){
  print("Error: --input is a required field.")
  quit()
}


if(exists("input_out_name") == F){
  # print(paste("Writing to ", input_genus, "_out", sep = ""))
  input_out_name <- paste(input_name, "_out", sep = "")
}




print("Starting.")
print(paste("Writing to ", input_out_name, sep = ""))
library(tidyverse)
library(ggplot2)
library(tjnFunctions)




models <- readLines(input_name)
models <- data.frame(text = models, stringsAsFactors = F)

models <- models%>%
  separate(col = text, into = c("id", "i1"), sep = " ", remove = F, extra = "merge")%>%
  separate(col = i1, into = c("protein", "i2"), sep = "\\[", remove = T, extra = "merge")%>%
  separate(col = i2, into = c("genome", "model"), sep = "\\]", remove = T, extra = "merge")

hypotheicals <- models%>%filter(grepl("hypothetical", protein) == T)
ORFs <- models%>%filter(grepl("ORF", protein) == T)
VPs <- models%>%filter(grepl("VP", protein) == T)
other <- models%>%filter(grepl("hypothetical", protein) == F)%>%filter(grepl("ORF", protein) == F)%>%filter(grepl("VP", protein) == F)


hypotheicals <- hypotheicals%>%bind_rows(ORFs)%>%bind_rows(VPs)

hypotheicals$protein <- rep(" ", nrow(hypotheicals))

models <- hypotheicals%>%bind_rows(other)


proteins <- models%>%group_by(model)%>%summarise(protein.list = paste(unique(protein), collapse = ","))
proteinCounts <- models%>%group_by(model)%>%select(model, protein)%>%unique()%>%summarise(protein.count = n())

genomes <- models%>%group_by(model)%>%summarise(genome.list = paste(unique(genome), collapse = ","))
genomeCounts <- models%>%group_by(model)%>%select(model, genome)%>%unique()%>%summarise(genome.count = n())


modelSummary <- proteins%>%left_join(proteinCounts)%>%left_join(genomes)%>%left_join(genomeCounts)%>%
  arrange(model)

write.table(modelSummary, input_out_name, quote = F, col.names = F, row.names = F, sep = "\t")
