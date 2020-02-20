#!/usr/bin/env Rscript

# getopts -----------------------------------------------------------------


suppressMessages(library('getopt'))


spec = matrix(c(
  'gff1' , 'r', 1, "character",
  'gff2', 'g', 1, "character",
  'help' , 'h', 0, "logical",
  'initial_data' , 'i', 0, "logical",
  'intergenic' , 'j', 0, "logical",
  'alignment' , 'a', 1, "character",
  'file_path', 'p', 2, "character",
  'out_name', 'o', 2, "character",
  'id1', 'x', 2, "character",
  'id2', 'y', 2, "character",
  'seq1', 's', 2, "character",
  'seq2', 't', 2, "character"
), byrow=TRUE, ncol=4)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat("sraPlotGenomeAlignment.R version 1.0\n")
  cat(" \n")
  cat("Use sraPlotGenomeAlignment.R <options> -s <sra plot file> -r <reference plot file> -a <alignment file>\n")
  cat(" \n")
  cat("Options:\n")
  cat("  -f <file path> The location of the other files and the output file\n")
  cat("  -r <reference gff file> The gff file that the rearranged file needs to match\n")
  cat("  -g <other gff file> The gff file to be reordered\n")
  cat("  -a <alignment file> The file that will be used to rearrange the sra file (This needs to include the extension). If missing the files will not be aligned.\n")
  cat("  -x <id1> The id for gff 1\n")
  cat("  -y <id2> The id for gff 2\n")
  cat("  -s <seq1 column> The column number for seq 1\n")
  cat("  -t <seq2 column> The column number for seq 2\n")
  cat("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input\n")
  cat("  -i <initial data> the input data is the output from combine_gff_files.R\n")
  cat("  -j <intergenic data> the file locations are different\n")
  cat("  \n")
  q(status=1)
}

if ( !is.null(opt$initial_data) ) {
  initial_data <- T
}else{
  initial_data <- F
}





if ( is.null(opt$gff2) ) {
  cat("Error: -g <other gff file> is required.\n")
  q(status=1)
}
if ( is.null(opt$gff1) ) {
  cat("Error: -r <reference gff file> is required.\n")
  q(status=1)
}

if ( is.null(opt$alignment) ) {
  align <- F
  cat("Warning: -a <alignment file> not specified. Files will not be aligned.\n")


}else{
  align <- T
}

# packages ----------------------------------------------------------------


suppressMessages(library(tidyverse))
suppressMessages(library(comparativeSRA))

# Test setup --------------------------------------------------------------
test_setup <- F
if(test_setup == T){
  if(initial_data == F){
opt$gff1 <- "~/phd/RNASeq/combined_gff_files/esch_1-2-3_merged.gff"
opt$gff2 <- "~/phd/RNASeq/combined_gff_files/GCA_000017745.1-GCA_900186905.1_merged.gff"
opt$alignment <- "~/phd/RNASeq/alignments/escherichia-salmonella.backbone"
opt$id1 <- "esch_1-2-3"
opt$id2 <- "GCA_900186905.1"
opt$out_name <- "esch_1-2-3-15"
opt$file_path <- "~/phd/RNASeq/combined_gff_files/"
align <- F
initial_data <- F
}else{
  initial_data <- T
opt$gff1 <- "~/phd/RNASeq/escherichia/GCA_000017745_data/GCA_000017745.1_new_calls.txt"
opt$gff2 <- "~/phd/RNASeq/escherichia/GCA_000017765.1_data/GCA_000017765.1_new_calls.txt"
opt$alignment <- "~/phd/RNASeq/escherichia/escherichia.backbone"
opt$id1 <- "GCA_000017745.1"
opt$id2 <- "GCA_000017765.1"
opt$out_name <- "escherichia_1-2"
}
}



# Defining variables ------------------------------------------------------



if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$id1 ) ) {  opt$id1 = opt$gff1 }
if ( is.null(opt$id2 ) ) { opt$id2 = opt$gff2 }
if ( is.null(opt$seq1 ) ) {  opt$seq1 = "1" }
if ( is.null(opt$seq2 ) ) { opt$seq2= "2" }

if(align){
  if ( is.null(opt$out_name ) ) { opt$out_name = paste(opt$alignment, "_", opt$seq1, "-", opt$seq2, sep = "") }
  
}else{
  placeholer_1 <- unlist(strsplit(opt$gff1, "_"))
  placeholer_2 <- unlist(strsplit(opt$gff2, "_"))
  if(placeholer_1[1] == placeholer_2[[1]]){
    if ( is.null(opt$out_name ) ) { opt$out_name = paste(placeholer_1[1], "_", placeholer_1[2], "-", placeholer_2[2], sep = "") }
    
  }else{
    if ( is.null(opt$out_name ) ) { opt$out_name = paste(opt$gff1,  "-", opt$gff2, sep = "") }
    
  }  
  
}

filePath <- opt$file_path

if(align){
if(grepl("/", opt$alignment) == F){
  if(grepl(".backbone", opt$alignment) == F){
  opt$alignment <- paste("~/phd/RNASeq/alignments/backbones/", opt$alignment, ".backbone", sep = "")
  }else{
    opt$alignment <- paste("~/phd/RNASeq/alignments/backbones/", opt$alignment, sep = "")
    
  }
}
}

# Main section ------------------------------------------------------------
if(initial_data == T){
  cat(paste("Analysing initial calls from ", "~/phd/RNASeq/new_calls/", opt$gff1, "_new_calls.txt and ", "~/phd/RNASeq/new_calls/", opt$gff2, "_new_calls.txt\n", sep = ""))
  
  gff1 <- read.table(paste("~/phd/RNASeq/new_calls/", opt$gff1, "_new_calls.txt", sep = ""), sep = "\t", header = T, as.is = T)
  gff2 <- read.table(paste("~/phd/RNASeq/new_calls/", opt$gff2, "_new_calls.txt", sep = ""), sep = "\t", header = T, as.is = T)

  
  ncRNAgff <- alignAndCombine(reference = opt$alignment,
                                      gff1 = gff1,
                                      gff2 = gff2,
                                      filenum1 = opt$id1,
                                      filenum2 = opt$id2,
                                      seqA = opt$s,
                                      seqB = opt$t,
                                quiet = T)
ncRNAgff <- ncRNAgff %>% mutate(set_val = 1)
mergedData <- mergeSRA(ncRNAgff = ncRNAgff,
                       filenum1 = opt$id1,
                       filenum2 = opt$id2,
                       initial_data = initial_data, 
                       align = T)

mergedData <- mergedData%>%mutate(change = ifelse(start < end, F, T))%>%
  mutate(start.tmp = end)%>%
  mutate(end.tmp = start)%>%
  mutate(start = ifelse(change == T, start.tmp, start))%>%
  mutate(end = ifelse(change == T, end.tmp, end))%>%
  select(-start.tmp, -end.tmp, -change)



mergedData <- mergedData%>%filter(!is.na(sequence))
mergedData[is.na(mergedData)] <- 0

file_id1 <- unlist(strsplit(opt$id1, "-"))
file_id2 <- unlist(strsplit(opt$id2, "-"))
file_id <- paste(unique(c(file_id1, file_id2)), collapse = "-")


tmp <- mergedData%>%
  mutate(id = NA)%>%
  mutate(set_val = NA)%>%
  mutate(file_id = file_id)
tmp[is.na(tmp)] <- 0

i <- 4
for(i in 1:nrow(tmp)){
  id1_list <- unlist(strsplit(tmp$id1[i], "-"))
  id2_list <- unlist(strsplit(tmp$id2[i], "-"))
  id_list <- unique(c(id1_list, id2_list))
  
  tmp$id[i] <- paste(id_list, collapse = "-")
  
  
  set_val_1 <- unlist(strsplit(tmp$set_val_1[i], "-"))
  set_val_2 <- unlist(strsplit(tmp$set_val_2[i], "-"))
  
  
  
  
  if(length(intersect(set_val_1, set_val_2)) > 0){
    set_val <- intersect(set_val_1, set_val_2)
  }else{
    set_val <- union(set_val_1, set_val_2)
  }
  tmp$set_val[i] <- paste(set_val, collapse = "-")
  
  
  
}




mergedData <- tmp[,c(1:13, (ncol(tmp) - 1):(ncol(tmp)), 18:(ncol(tmp) - 2) )]
mergedData <- mergedData%>%mutate(V1 = set_val)
colnames(mergedData)[ncol(mergedData)] <- paste(opt$out_name)


}else{

  
  if(is.null(opt$intergenic)){
  gff1Dat <- read.table(paste("~/phd/RNASeq/combined_gff_files/", opt$gff1, "_merged.gff", sep = ""), sep = "\t", header = T, as.is = T)
  gff2Dat <- read.table(paste("~/phd/RNASeq/combined_gff_files/", opt$gff2, "_merged.gff", sep = ""), sep = "\t", header = T, as.is = T)
  }else{
    gff1Dat <- read.table(paste("~/phd/RNASeq/combined_gff_files_random/", opt$gff1, "_merged.gff", sep = ""), sep = "\t", header = T, as.is = T)
    gff2Dat <- read.table(paste("~/phd/RNASeq/combined_gff_files_random/", opt$gff2, "_merged.gff", sep = ""), sep = "\t", header = T, as.is = T)
  }

  gff1Working <- gff1Dat %>% mutate(row_numbers = as.character(row_numbers))
  gff2Working <- gff2Dat %>% mutate(row_numbers = as.character(row_numbers))
  

  filenum1 <- gff1Working$file_id[1]
  filenum2 <- gff2Working$file_id[1]

  if(align == T){
   ncRNAgff <- alignAndCombine(reference = opt$alignment,
                               gff1 = gff1Working,
                               gff2 = gff2Working,
                               filenum1 = filenum1,
                               filenum2 = filenum2,
                               seqA = 1,
                               seqB = 2)

   ncRNAgff <- ncRNAgff%>%select(-changed)%>%unique()
  ncRNAgff[is.na(ncRNAgff)] <- "0"
  }else{
    ncRNAgff <- gff1Working%>%bind_rows(gff2Working)
    ncRNAgff[is.na(ncRNAgff)] <- 0
  }
  
  if(test_setup == T){
  mergeSRAData <- list(ncRNAgff = ncRNAgff, filenum1 = filenum1, filenum2 = filenum2, initial_data = initial_data, align = align)
  save(mergeSRAData, file = "~/bin/r_git/R/mergeSRAData.Rda")
  
  }
  
  mergedData <- mergeSRA(ncRNAgff = ncRNAgff,
                         filenum1 = filenum1,
                         filenum2 = filenum2,
                         align = align, 
                         initial_data = F)


  mergedData <- mergedData%>%mutate(change = ifelse(start < end, F, T))%>%
    mutate(start.tmp = end)%>%
    mutate(end.tmp = start)%>%
    mutate(start = ifelse(change == T, start.tmp, start))%>%
    mutate(end = ifelse(change == T, end.tmp, end))%>%
    select(-start.tmp, -end.tmp, -change)



  mergedData <- mergedData%>%filter(!is.na(sequence))
  mergedData[is.na(mergedData)] <- 0

  file_id1 <- unlist(strsplit(filenum1, "-"))
  file_id2 <- unlist(strsplit(filenum2, "-"))
  file_id <- paste(unique(c(file_id1, file_id2)), collapse = "-")


  tmp <- mergedData%>%
    mutate(id = NA)%>%
    mutate(set_val = NA)%>%
    mutate(file_id = file_id)
  tmp[is.na(tmp)] <- 0

  #i <- 4
for(i in 1:nrow(tmp)){
  id1_list <- unlist(strsplit(tmp$id1[i], "-"))
  id2_list <- unlist(strsplit(tmp$id2[i], "-"))
  id_list <- unique(c(id1_list, id2_list))

  tmp$id[i] <- paste(id_list, collapse = "-")


  set_val_1 <- unlist(strsplit(tmp$set_val_1[i], "-"))
  set_val_2 <- unlist(strsplit(tmp$set_val_2[i], "-"))




  if(length(intersect(set_val_1, set_val_2)) > 0){
    set_val <- intersect(set_val_1, set_val_2)
  }else{
    set_val <- union(set_val_1, set_val_2)
  }
  tmp$set_val[i] <- paste(set_val, collapse = "-")



}
  fitchTest <- tmp %>% select(id, set_val) %>% mutate(fitch = 0, prop = 0)
  
  for(i in 1:nrow(fitchTest)){
    files_1 <- c()
    files_all <- c()
    id_set <- unlist(strsplit(as.character(fitchTest$id[i]), "-"))
    for(j in 1:length(id_set)){
      number <- unlist(strsplit(as.character(id_set[j]), "_"))
      file_id <- paste(number[1:2], collapse = "_")
      number <- number[3]
      files_all <- c(files_all, file_id)
      if(number > 0){
        files_1 <- c(files_1, file_id)
      }
    }
    files_1 <- unique(files_1)
    files_all <- unique(files_all)
    if(align == F){
       fitchTest$fitch[i] <- fitchTest$set_val[i]##use for random only thgis needs removing
      #fitchTest$fitch[i] <- ifelse(length(files_1) == 1, "0-1", ifelse(length(files_1) == 0, "0", "1"))
    }else{
      fitchTest$fitch[i] <- fitchTest$set_val[i]
    }
    fitchTest$prop[i] <- paste(length(files_1), length(files_all), sep = "-")
  }
  fitchTest <- fitchTest %>% select(-set_val)
  
  
  mergedData <- tmp[,c(1:13, (ncol(tmp) - 1):(ncol(tmp)), 18:(ncol(tmp) - 2) )]
  mergedData <- mergedData %>% full_join(fitchTest, by = "id")
  
  
  if(align == F){
    mergedData <- mergedData%>%mutate(set_val = fitch)
  }
    mergedData <- mergedData%>%mutate(V1 = set_val)

  
  mergedData <- mergedData%>%mutate(V2 = prop)
  colnames(mergedData)[ncol(mergedData) - 1] <- paste(opt$out_name, "fitch", sep = "-")
  colnames(mergedData)[ncol(mergedData)] <- paste(opt$out_name, "prop", sep = "-")
  
  mergedData <- mergedData %>% select(-fitch, -prop)
  mergedData[is.na(mergedData)] <- 0
  
}

cat(paste("Writing the output to ", filePath, "/", opt$out_name, "_merged.gff\n", sep = ""))
write.table(x = mergedData, file = paste(filePath, "/", opt$out_name, "_merged.gff", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")


