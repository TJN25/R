#!/usr/bin/env Rscript
suppressMessages(library('getopt'))


##Getopts setup ####
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'quiet' , 'q', 0, "logical",
  'file_path', 'p', 2, "character",
  'out_name', 'o', 2, "character",
  'name', 'i', 1, "character",
  'sequence', 's', 1, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat("combine_gff_files.R version 1.0\n")
  cat(" \n")
  cat("Use combine_gff_files.R <options> -f <files>\n")
  cat(" \n")
  cat("Options:\n")
  cat("  -f <files> The gff files\n")
  cat("  -s <stranded data> The data is stranded\n")
  cat("  -q <quiet> Do not print any updates\n")
  cat("  -p <file path> The location of the other files and the output file\n")
  cat("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input\n")
  q(status=1)
}

if ( is.null(opt$name) ) {
  cat("Error: -i <input> is required.\n")
  q(status=1)
}






suppressMessages(library(tidyverse))
suppressMessages(library(tjnFunctions))
suppressMessages(library(IRanges))
suppressMessages(library(stringi))

if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$out_name ) ) { opt$out_name = opt$name }
filePath <- opt$file_path

# Functions ---------------------------------------------------------------
##this is just copied from getRNASequences and will need to be adapted/tested


getSequencesFast <- function(ids_lookup, fasta){
  tmp <- ids_lookup
  tmp <- tmp%>%mutate(sequence = "")
  
  for(i in 1:nrow(tmp)){
    printRemaining(i = i, length = nrow(tmp), increment = 5)
    # x <- system(command = paste("sed -n '", (as.numeric(tmp[i, 2]) + 1), "p'",
    #                             " ~/phd/RNASeq/new_calls/",
    #                             tmp[i, 3],
    #                             "_new_calls.txt", sep = ""), intern = T)
    # x <- unlist(strsplit(as.character(x), split = "\t"))
    #print(x[3:4])
    
    
    try( tmp[i, ncol(tmp)] <- substr(fasta, ids_lookup[i,1], ids_lookup[i,2]))
    
  }
  
  return(tmp)
}


# working -----------------------------------------------------------------



overlaps <- read.table(file = paste(filePath, "/", opt$name, sep = ""), header = T, sep = "\t")
#overlaps <- read.table(file = "~/phd/RNASeq/random_sequences/GCA_001874505.1_random_data.txt", header = T, sep = "\t")

overlaps <- overlaps %>% filter(type == "intergenic")

#fileName <- "~/phd/RNASeq/Enterobacter/GCA_001874505.1.data/GCA_001874505.1.fna"
#fasta <- readLines(fileName)

fasta <- readLines(opt$sequence)
fasta <- data.frame(sequence = fasta)

fasta <- fasta %>% filter(grepl(">", sequence) == F) %>% 
  mutate(length = str_length(sequence)) %>% 
  mutate(left = 1) %>% 
  mutate(length.total = cumsum(length)) %>% 
  mutate(right =  length.total) %>% 
  mutate(left = right - length + 1) %>% 
  select(left, sequence, right)
fasta.list <- paste(as.vector(fasta$sequence), collapse = "")

sequences <- getSequencesFast(overlaps[,], fasta.list)

sequences <- sequences %>% 
  mutate(tmp.seq = str_replace_all(string = sequence, pattern = "A", replacement = "B")) %>% 
  mutate(tmp.seq = str_replace_all(string = tmp.seq, pattern = "C", replacement = "D")) %>% 
  mutate(tmp.seq = str_replace_all(string = tmp.seq, pattern = "G", replacement = "H")) %>% 
  mutate(tmp.seq = str_replace_all(string = tmp.seq, pattern = "T", replacement = "U")) %>% 
  mutate(tmp.seq = str_replace_all(string = tmp.seq, pattern = "B", replacement = "T")) %>% 
  mutate(tmp.seq = str_replace_all(string = tmp.seq, pattern = "D", replacement = "G")) %>% 
  mutate(tmp.seq = str_replace_all(string = tmp.seq, pattern = "H", replacement = "C")) %>% 
  mutate(tmp.seq = str_replace_all(string = tmp.seq, pattern = "U", replacement = "A")) %>% 
  mutate(tmp.seq = stri_reverse(tmp.seq))




sequences <- sequences %>% mutate(id = paste(">", opt$name, "_random_intergenic_", row_number(), sep = "")) %>% 
  select(id, sequence)

cat(paste("Writing the output to ", filePath, "/", opt$out_name, "_random_intergenic_sequences.fna\n", sep = ""))
write.table(x = sequences, file = paste(filePath, "/", opt$out_name, "_random_intergenic_sequences.fna", sep = ""), row.names = F, col.names = F, quote = F, sep = "\n")

#
