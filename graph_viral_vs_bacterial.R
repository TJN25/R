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
  print("      --bacteria (bacterial info file),", quote = F )
  print("      --phage (phage info file),", quote = F)
  print("      --genus (genus name), ", quote = F)
  print("      --wd (working directory), ", quote = F)
  print(" ", quote = F)
  print("   Other fields:", quote = F)
  print("      --out_name (output file prefix),", quote = F)
  print("      --out_type (plot or write),", quote = F)
  print("      --aa (amino acid freq file for an order codon plot),", quote = F)
  print(" ", quote = F)
  print("   Needed for GC plots:", quote = F)
  print("      --trna (bacterial tRNA results file)", quote = F)
  print("      --cds (bacterial CDS results file)", quote = F)
  quit()
}
i <- 11
#print(args)
for(i in 1:length(args)){
  input_type <- args[i]
  if(input_type == "--help"){
    print("   Required fields:", quote = F)
    print("      --bacteria (bacterial info file),", quote = F )
    print("      --phage (phage info file),", quote = F)
    print("      --genus (genus name), ", quote = F)
    print("      --wd (working directory), ", quote = F)
    print(" ", quote = F)
    print("   Other fields:", quote = F)
    print("      --out_name (output file prefix),", quote = F)
    print("      --out_type (plot or write),", quote = F)
    print("      --aa (amino acid freq file for an order codon plot),", quote = F)
    print(" ", quote = F)
    print("   Needed for GC plots:", quote = F)
    print("      --trna (bacterial tRNA results file)", quote = F)
    print("      --cds (bacterial CDS results file)", quote = F)
    quit()

  }else if(input_type == "--bacteria"){
    input_bacteria <- args[i + 1]
    print(args[i + 1])
    print(input_bacteria)
  }else if(input_type == "--phage"){
    input_phage <- args[i + 1]
  }else if(input_type == "--genus"){
    input_genus <- args[i + 1]
  }else if(input_type == "--out_name"){
    input_out_name <- args[i + 1]
  }else if(input_type == "--aa"){
    input_aa <- args[i + 1]
  }else if(input_type == "--wd"){
    input_wd <- args[i + 1]
  }else if(input_type == "--trna"){
    input_trna <- args[i + 1]
  }else if(input_type == "--cds"){
    input_cds <- args[i + 1]
  }else if(input_type == "--out_type"){
    input_out_type <- args[i + 1]
  }

}

##Required
if(exists("input_genus") == F){
  print("Error: --genus is a required field.")
  quit()

}

##Optional
if(exists("input_bacteria") == F){
  input_bacteria <- paste(input_genus, "_all.txt", sep = "")
  #print("Warning: --bacteria not provided. Trying genus name to open file")

}


if(exists("input_phage") == F){
  input_phage <- "phage_info_19-06-2018.txt"
  #print("Warning: --phage not provided. Using default")

}

if(exists("input_wd") == F){

  input_wd <- "/Users/thomasnicholson/phd/virus_host_prediction/"

  print("Warning: Working directory set to: /Volumes/2TTJN/Virus_Host_Prediction/")
  print("Use --wd to change this.")
}else{
  print(input_wd)
}

if(exists("input_out_name") == F){
  # print(paste("Writing to ", input_genus, "_out", sep = ""))
  input_out_name <- paste(input_genus, "_out", sep = "")
}

if(exists("input_aa") == F){
  input_aa <- paste(input_genus, "_all.amino_acids.txt", sep = "")
  #print("Warning: Amino acids are not ordered by abundance as the file for ordering has not been provided")
}

if(exists("input_cds") == F){
  input_cds <- paste(input_genus, "_CDS_output.txt", sep = "")
  #  print("Warning: GC plot will not be produced as the CDS file is missing")
}

if(exists("input_trna") == F){
  input_trna <- paste(input_genus, "_tRNA_output.txt", sep = "")
  #  print("Warning: GC plot will not be produced as the tRNA file is missing")
}

if(exists("input_out_type") == F){
  #print(paste("Writing plot to ", input_out_name, sep = ""))
  input_out_type = "write"
}



print("Starting.")
print(paste("Writing to ", input_out_name, sep = ""))
library(tidyverse)
library(ggplot2)
library(tjnFunctions)




setwd(input_wd)






## Import and Setup ####
codonLookup <- read.table("/Users/thomasnicholson/phd/virus_host_prediction/codon_lookup.txt", header = F, sep = "\t", comment.char = "", quote = "", as.is = T, fill = T)

try(aminoAcidFreq <- read.table(input_aa, header = T, sep = "\t", comment.char = "", quote = "", as.is = T, fill = T), silent=TRUE)

if(exists("aminoAcidFreq") == F){
  input_aa <- paste(input_genus, "_out.amino_acids.txt", sep = "")
  #print("Warning: Amino acids are not ordered by abundance as the file for ordering has not been provided")
  try(aminoAcidFreq <- read.table(input_aa, header = T, sep = "\t", comment.char = "", quote = "", as.is = T, fill = T), silent=TRUE)
  if(exists("aminoAcidFreq") == F){
    print("Warning: Amino acids are not ordered by abundance as the file for ordering has not been provided")
    print("Cannot find an amino acid file. Specify with --aa")
    order_aa = F
  }
}

if(order_aa == T){
  aminoSummary <- createAASummary(aa = aminoAcidFreq, cl = codonLookup)
}


try(bacterial <- read.table(input_bacteria, header = F, sep = "\t", comment.char = "", quote = "", as.is = T, fill = T), silent=TRUE)


if(exists("bacterial") == F){
  input_bacteria <- paste(input_genus, "_out.txt", sep = "")
  #print("Warning: --bacteria not provided. Trying genus name to open file")
  try(bacterial <- read.table(input_bacteria, header = F, sep = "\t", comment.char = "", quote = "", as.is = T, fill = T), silent=TRUE)

  if(exists("bacterial") == F){
    print("Cannot find a bacterial file. Specify with --bacteria")
    quit()
  }
}

bacterial <- bacterialColumnNames(bac = bacterial)



if(gc_plotting == T){
  try(tRNA <- read.table(input_trna, as.is = T, comment.char = "", quote = "", header = F), silent=TRUE)

  if(exists("tRNA") == F){
    input_trna <- paste(input_genus, "_out.tRNA.txt", sep = "")
    try(tRNA <- read.table(input_trna, as.is = T, comment.char = "", quote = "", header = F), silent=TRUE)
    if(exists("tRNA") == F){
      print("Warning: GC plot will not be produced as the tRNA file is missing")
      print("Cannot find a tRNA file. Specify with --trna")
      gc_plotting = F
    }

  }

  # tRNA <- read.table("Bacillus_tRNA_output.txt")
  try(CDS <- read.table(input_cds, sep = "|", as.is = T, comment.char = "", quote = "", header = F), silent = T)
  if(exists("CDS") == F){
    input_cds <- paste(input_genus, "_out.CDS.txt", sep = "")
    #  print("Warning: GC plot will not be produced as the CDS file is missing")
    try(CDS <- read.table(input_cds, sep = "|", as.is = T, comment.char = "", quote = "", header = F), silent = T)
    if(exists("CDS") == F){
      print("Warning: GC plot will not be produced as the CDS file is missing")
      print("Cannot find a CDS file. Specify with --cds")
      gc_plotting = F
    }
  }
}

try(phage <- read.table(input_phage, header = F, sep = "\t", comment.char = "", quote = "", as.is = T, fill = T), silent = T)
if(exists("phage") == F){
  print("Cannot find a phage file. Specify with --phage")
  quit()
}



phage <- phageColumnNames(ph = phage)
phage <- phage%>%filter(genome_genera == input_genus)

bacCodons <- bacterialCodonsFormat(bac = bacterial)
phageCodons <- phageCodonsFormat(ph = phage)

bacDinucleotides <- bacterialDinucleotideFormat(bac = bacterial)
phageDinucleotides <- phageDinucleotideFormat(ph = phage)

codons <- codonSummary(ph = phageCodons, bac = bacCodons, codon_order = order_aa)

dinucleotides <- dinucleotideSummary(ph = phageDinucleotides, bac = bacDinucleotides)


plotCodons(cc = codons, output = input_out_type)

plotDinucleotides(cc = dinucleotides, output = input_out_type)

if(gc_plotting == T){
  CDS$V2 <- as.numeric(CDS$V2)
  plotGC(ph = phage, tr = tRNA, cd = CDS, output = input_out_type)
}

print("Done.")

