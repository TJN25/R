#!/usr/bin/env Rscript
#>> This is not a working script. It contains a number of standalone functions that are called within the script
##Functions ####
assignConservationLevel <- function(ids_lookup, main_col = 7, genera_col, species_col, any_col = c(7:ncol(ids_lookup))){
  ids_lookup <- ids_lookup%>%mutate(type = "")
  for(i in 1:nrow(ids_lookup)){
    if("1" %in% ids_lookup[i, main_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "main_conserved"
    }else if("1" %in% ids_lookup[i, genera_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "genera_conserved"
    }else if("1" %in% ids_lookup[i, species_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "species_conserved"
    }else if("0-1" %in% ids_lookup[i, main_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "main_0-1"
    }else if("0-1" %in% ids_lookup[i, genera_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "genera_0-1"
    }else if("0-1" %in% ids_lookup[i, species_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "species_0-1"
    }else if("1" %in% ids_lookup[i, any_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "any_conserved"
    }else if("0-1" %in% ids_lookup[i, any_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "any_0-1"
    }

  }
  return(ids_lookup)
}

getSequences <- function(ids_lookup, keep.type, new_only = F){
  if(keep.type == "all"){
    tmp <- ids_lookup
  }else{
    tmp <- ids_lookup%>%filter(type == keep.type)
  }
  if(new_only == T){
    tmp <- tmp%>%filter(new_feature ==T)
  }

  x <- tmp%>%group_by(id)%>%summarise(count = n())
  cat(paste(nrow(x), " ids kept.\n", sep = ""))

  tmp <- tmp%>%mutate(sequence = "")
i <- 2
  for(i in 1:nrow(tmp)){
    print(i)
    #printRemaining(i = i, length = nrow(tmp), increment = 5)
    x <- system(command = paste("sed -n '", (as.numeric(tmp[i, 2]) + 1), "p'",
                                " ~/phd/RNASeq/random_calls/",
                                tmp[i, 3],
                                "_new_calls.txt", sep = ""), intern = T)
    x <- unlist(strsplit(as.character(x), split = "\t"))
    #print(x)
    genome <- x[6]
    genome <- paste(unlist(strsplit(genome, "\\."))[1:2], collapse = ".")
    genome <- paste(unlist(strsplit(genome, "_"))[1:2], collapse = "_")
    
    
    
        #print(x[3:4])


    tmp[i, ncol(tmp)] <- system(paste("test_string=`grep -v ^'>' ~/phd/RNASeq/sequences/",
                                      genome,
                                      ".fna | tr -d '\n'`; echo ${test_string:",
                                      x[3], ":", as.numeric(x[4]) - as.numeric(x[3]),
                                      "}", sep = ""), intern = T)
  }

  return(tmp)
}

proportionConserved <- function(countsDat, column_num = 2, rows = c(3:nrow(countsDat))){
  total <- sum(countsDat[,column_num])
  prop <- sum(countsDat[rows, column_num])/total
  return(prop)
}

writeSequences <- function(sequence_list, out_name = "random_sra_enterics-serratia_", mainDir = "~/phd/RNASeq/random_calls/srna_all/"){
  for(i in 1:length(unique(sequence_list$id))){
    for(k in unique(sequence_list_all$type)){
      subDir <- k
      if (file.exists(paste(mainDir,subDir, sep = "")) == F){
        dir.create(file.path(mainDir, subDir))
      }
      if(file.exists(paste(mainDir,subDir, "/new", sep = "")) == F){
        dir.create(file.path(mainDir, subDir, "/new"))
      }
      if(file.exists(paste(mainDir,subDir, "/known", sep = "")) == F){
        dir.create(file.path(mainDir, subDir, "/known"))
      }
    }


    tmp <- sequence_list%>%filter(id == unique(sequence_list$id)[i])
    tmp <- tmp%>%mutate(seq_length = nchar(sequence))
    medianLength <- median(tmp$seq_length)

    for(j in 1:nrow(tmp)){
      mat <- matrix(ncol = 1, nrow = 2)
      mat[1,1] <- paste(">", tmp[j, 3], "|", out_name, i, "|", tmp[j, ncol(tmp)], sep = "")
      mat[2,1] <- tmp[j,(ncol(tmp) -1)]
      subDir <- as.character(tmp[j, (ncol(tmp) - 2)])
      if(j == 1){
        write.table(mat, file = paste(mainDir, subDir, ifelse(tmp[j, 6], "/new/", "/known/"),
                                      out_name, i, "_",
                                      round((1 - ((max(tmp$seq_length) - min(tmp$seq_length))/max(tmp$seq_length))), 2),
                                      ".fasta", sep = ""),
                    row.names = F, col.names = F, quote = F, sep = "\t", append = F)
      }else{
        write.table(mat, file = paste(mainDir, subDir, ifelse(tmp[j, 6], "/new/", "/known/"),
                                      out_name, i, "_",
                                      round((1 - ((max(tmp$seq_length) - min(tmp$seq_length))/max(tmp$seq_length))), 2),
                                      ".fasta", sep = ""),
                    row.names = F, col.names = F, quote = F, sep = "\t", append = T)
      }
    }
  }

}

rewriteSequencesFromAlignments <- function(sequence, out_name){
  mcl_seqs <- readLines(sequence)
  mcl_dat <- data.frame(id = "", sequence = "")
  current_seq = F
  seq <- ""
  for(i in 1:length(mcl_seqs)){
    if(i == length(mcl_seqs)){
      if(substr(mcl_seqs[i], 1,1) != ">"){
        seq <- c(seq, mcl_seqs[i])
        tmp <- data.frame(id = id, sequence = paste(seq, collapse = ""))
        mcl_dat <- mcl_dat%>%bind_rows(tmp)
      }
    }
    if(substr(mcl_seqs[i], 1,1) == ">"){
      if(current_seq == T){
        tmp <- data.frame(id = id, sequence = paste(seq, collapse = ""))
        mcl_dat <- mcl_dat%>%bind_rows(tmp)
        seq <- ""
      }
      current_seq <- T
      id <- mcl_seqs[i]

    }else{
      seq <- c(seq, mcl_seqs[i])
    }
  }

  mcl_dat <- mcl_dat%>%mutate(length = nchar(sequence))

  mcl_dat <- mcl_dat%>%mutate(start = 0, end = 0)

  lenMax <- max(mcl_dat$length)
  i <- 1
  for(i in 1:nrow(mcl_dat)){
    m <- stri_locate_all(str = mcl_dat[i,2], regex = "-+")[[1]]
    missing.length <- m[ , "end"] - m[ , "start"] + 1
    tmp <- as.data.frame(cbind(m, missing.length))
    tmpStart <- tmp%>%filter(start == 1)
    tmpEnd <- tmp%>%filter(end == lenMax)
    if(nrow(tmpStart) > 0){
      mcl_dat[i, 4] <- tmpStart[1,3]

    }
    if(nrow(tmpEnd) > 0){
      mcl_dat[i, 5] <- tmpEnd[1,3]
    }
  }
  mcl_dat <- mcl_dat%>%filter(id != "")

  n <- length(mcl_dat$start)
  start_cutoff <- abs(sort(-mcl_dat$start,partial=n-1)[n-1])

  n <- length(mcl_dat$end)
  end_cutoff <- abs(sort(-mcl_dat$end,partial=n-1)[n-1])

  for(i in 1:nrow(mcl_dat)){

    mat <- matrix(ncol = 1, nrow = 2)
    mat[1,1] <- mcl_dat[i,1]
    tmp <- mcl_dat%>%mutate(seq_short = substr(x = sequence, start = start_cutoff, stop = (length - end_cutoff)))
    for(j in 1:nrow(tmp)){
      tmp[j,6] <- str_remove_all(tmp[j, 6], "-")
    }




    mat[2,1] <- tmp[i,ncol(tmp)]


    if(i == 1){
      write.table(mat, file = out_name,
                  row.names = F, col.names = F, quote = F, sep = "\t", append = F)
    }else{
      write.table(mat, file = out_name,
                  row.names = F, col.names = F, quote = F, sep = "\t", append = T)
    }

  }
}


#
##getopts setup ####
suppressMessages(library('getopt'))


spec = matrix(c(
  'gff2', 'g', 1, "character",
  'help' , 'h', 0, "logical",
  'initial_data' , 'i', 0, "logical",
  'gff1' , 'r', 1, "character",
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

suppressMessages(library(tidyverse))
suppressMessages(library(tjnFunctions))
library(stringi)
#
##working section ####


dat <- read.table("~/phd/RNASeq/random_calls/random_enterics-serratia_merged.gff", header = T, sep = "\t",
                  comment.char = "", quote = "", as.is = T )

dat[dat == "1-1"] <- "1"
dat[dat == "1-0"] <- "0-1"


#dat <- dat %>% filter(feature != "CDS")
dat <- dat %>% filter(feature == "intergenic")

dat_files <- unlist(strsplit(dat$file_id[1], "-"))
dat_ids <- dat$id

ids_lookup <- data.frame(id = "", row = "0", genome = "")
#colnames(ids_lookup) <- c("id", dat_files)
#ids_lookup[,1] <- dat_ids
#ids <- as.data.frame(ids)
i <- 1
j <- 1
for(i in 1:length(dat_ids)){
  uid <- dat_ids[i]
  uid <- unlist(strsplit(uid, "-"))
  for(j in 1:length(uid)){
    x <- paste(unlist(strsplit(uid[j], "_"))[1:3], collapse = "_")
    y <- unlist(strsplit(uid[j], "_"))[4]
    if(y != "0"){
      tmp <- data.frame(id = dat_ids[i], row = y, genome = x)
      ids_lookup <- ids_lookup%>%bind_rows(tmp)
    }

  }
}


datSmall <- dat[,c(8, 10, 12, 14, 15, 17:ncol(dat))]
ids_lookup <- ids_lookup%>%filter(row != "0")%>%left_join(datSmall, by = "id")


##classify each sRNA by the highest level of conservation found

main_col <- 7
genera_col <- c(14, 19, 20, 22)
species_col <- c(13:15, 18:22)
any_col <- c(7:ncol(ids_lookup))


ids_lookup <- assignConservationLevel(ids_lookup = ids_lookup, main_col = main_col, genera_col = genera_col,
                                      species_col = species_col, any_col = any_col)
ids_lookup %>% group_by(type) %>% summarise(count = n())



########


#Remove the CDS before further analysis

########

# sequence_list_main <- getSequences(ids_lookup, columns = c(7), new_only = F)
# sequence_list_genera_all <- getSequences(ids_lookup, columns = c(13:15, 18:22), new_only = F)
# sequence_list_all_new <- getSequences(ids_lookup, columns = c(7:ncol(ids_lookup)), new_only = T)


sequence_list_genera_conserved <- getSequences(ids_lookup, keep.type = "genera_conserved", new_only = F)
sequence_list_all <- getSequences(ids_lookup, keep.type = "all", new_only = F)


write.table(sequence_list_all, "~/phd/RNASeq/random_calls/srna_all/sequence_list.txt", sep = "\t", quote = F, row.names = F)
#ids_seq <- ids_seq%>%mutate(ids_short = paste("sra_enterics-serratia_", row_number(), sep = ""))
i <- 1
j <- 1

sequence_list_all <- read.table("~/phd/RNASeq/random_calls/srna_all/sequence_list.txt", header = T, sep = "\t", as.is = T)

sequence_list_all_srna <- sequence_list_all%>%filter(nchar(sequence) >= 30, nchar(sequence) <= 500)
write.table(sequence_list_all_srna, "~/phd/RNASeq/random_calls/srna_all/sequence_list_30-500.txt", sep = "\t", quote = F, row.names = F)


tmp <- sequence_list_all%>%mutate(seq_length = nchar(sequence))

ggplot() +
  geom_histogram(data = tmp%>%filter(seq_length < 500),
                 aes(x = seq_length, y = ..count.., fill = new_feature),
                 binwidth = 10)


count_30 <- sequence_list_all_srna%>%group_by(type)%>%summarise(count.30 = n())
count_all <- sequence_list_all%>%group_by(type)%>%summarise(count.all = n())
count_all <- count_all%>%full_join(count_30, by = "type")
count_all <- count_all%>%mutate(prop.30 = round(count.30/count.all, 2))

conservedGroupsAll <- sequence_list_all%>%group_by(id)%>%summarise(group.size.all = n())
id_list_all <- sequence_list_all%>%select(-row, -genome, -sequence)%>%unique()
conservedGroupsAll <- conservedGroupsAll%>%left_join(id_list_all, by = "id")

conservedGroupssRNA <- sequence_list_all_srna%>%group_by(id)%>%summarise(group.size.srna = n())
id_list_sRNA <- sequence_list_all_srna%>%select(-row, -genome, -sequence)%>%unique()
conservedGroupssRNA <- conservedGroupssRNA%>%left_join(id_list_all, by = "id")





count_30 <- conservedGroupssRNA%>%group_by(type)%>%summarise(count.30 = n())
count_all <- conservedGroupsAll%>%group_by(type)%>%summarise(count.all = n())
count_all <- count_all%>%full_join(count_30, by = "type")
count_all <- count_all%>%mutate(prop.30 = round(count.30/count.all, 2))


count_30_known <- conservedGroupssRNA%>%filter(new_feature ==F)%>%group_by(type)%>%summarise(count.30 = n())
count_all_known <- conservedGroupsAll%>%filter(new_feature ==F)%>%group_by(type)%>%summarise(count.all = n())
count_all_known <- count_all_known%>%full_join(count_30_known, by = "type")
count_all_known <- count_all_known%>%mutate(prop.30 = round(count.30/count.all, 2))


proportionConserved(count_all_known, column_num = 2, rows = c(3:nrow(count_all)))




#writeSequences(sequence_list = sequence_list_all_new, out_name = "new_only/sra_enterics-serratia_")
#writeSequences(sequence_list = sequence_list_main, out_name = "sra_enterics-serratia_")
#writeSequences(sequence_list = sequence_list_genera_all, out_name = "all/sra_enterics_genera_")

writeSequences(sequence_list = sequence_list_all_srna, out_name = "sra_enterics-serratia_Random-", mainDir = "~/phd/RNASeq/random_calls/srna_all/")








##take the alignments in and shorten the longer sequences to keep only the overlapping section
sequence_new_only <- "~/phd/RNASeq/sRNAs/new_only/mcl/sra_enterics-serratia_11_NEW_0.04.fasta.mcl"
out_name <- "~/phd/RNASeq/sRNAs/new_only/sra_enterics-serratia_11_NEW_0.04.fasta"




rewriteSequencesFromAlignments(sequence = sequence_new_only, out_name = out_name)



#

#####
