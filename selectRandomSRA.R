#!/usr/bin/env Rscript
suppressMessages(library('getopt'))


##Getopts setup ####
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'quiet' , 'q', 0, "logical",
  'file_path', 'p', 2, "character",
  'out_name', 'o', 2, "character",
  'name', 'i', 1, "character",
  'method', 'm', 1, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)
if ( !is.null(opt$help) ) {
  cat("selectRandomSRA.R version 1.1\n")
  cat(" \n")
  cat("Use selectRandomSRA. <options> -i <name>\n")
  cat(" \n")
  cat("Options:\n")
  cat("  -i <name> The GCA accession to work with. Do not inclue the gff file extension. The default is the same as the sra input\n")
  cat("  -m <method> Either 'shuffled' or 'shifted' (Shuffled in default)\n")
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
suppressMessages(library(comparativeSRA))
suppressMessages(library(IRanges))
suppressMessages(library(stringi))

if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$out_name ) ) { opt$out_name = opt$name }
if ( is.null(opt$method ) ) { opt$method = 'shuffled' }
opt$method <- tolower(opt$method)
if(opt$method != 'shuffled' && opt$method != 'shifted'){
  opt$method <- 'shuffled'
  cat("Incorrect method selected. Using 'shuffled")
}

filePath <- opt$file_path

if (file.exists(paste(filePath,"/random_sequences", sep = "")) == F){
    dir.create(file.path(filePath,"/random_sequences"))
    if(is.null(opt$quiet)){
    cat(paste(filePath,"/random_sequences made\n", sep = ""))
    }
  }else{
    if(is.null(opt$quiet)){

    cat(paste(filePath,"/random_sequences exists\n", sep = ""))
    }
  }


# Functions ---------------------------------------------------------------
##this is just copied from getRNASequences and will need to be adapted/tested
getSequences <- function(ids_lookup){
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


   try( tmp[i, ncol(tmp)] <- system(paste("test_string=`grep -v ^'>' ~/phd/RNASeq/sequences/",
                                      "GCA_000017765.1",
                                      ".fna | tr -d '\n'`; echo ${test_string:",
                                      ids_lookup[i,1], ":", as.numeric(ids_lookup[i,2]) - as.numeric(ids_lookup[i,1]),
                                      "}", sep = ""), intern = T))
  }

  return(tmp)
}


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


sra <- read.table(paste(filePath, "/", opt$name, "_new_calls.txt", sep = ""), header = T, sep = "\t", quote = "")
gff <- read.table(paste(filePath, "/", opt$name, ".gff", sep = ""), header = F, sep = "\t", quote = "")
colnames(gff) <- c("sequence", "source", "feature", "start", "end", "score", "strand", "phase", "Atrribute")

#fileName <- "~/phd/RNASeq/escherichia/GCA_000017765.1_data/GCA_000017765.1.fna"

#fasta <- readLines(fileName)

fasta <- readLines(paste(filePath, "/", opt$name, ".fna", sep = ""))
fasta <- data.frame(sequence = fasta)


fasta <- fasta %>% filter(grepl(">", sequence) == F) %>%
  mutate(length = str_length(sequence)) %>%
  mutate(left = 1) %>%
  mutate(length.total = cumsum(length)) %>%
  mutate(right =  length.total) %>%
  mutate(left = right - length + 1) %>%
  select(left, sequence, right)
fasta.list <- paste(as.vector(fasta$sequence), collapse = "")

lengths <- sra$end - sra$start
strands <- sra$strand

sdLength <- sd(lengths)
meanLength <- mean(lengths)

propStrand <- prop.table(table(strands))


if(opt$method == "shuffled"){

  if (file.exists(paste(filePath,"/random_sequences/shuffled", sep = "")) == F){
    dir.create(file.path(filePath,"/random_sequences/shuffled"))
    if(is.null(opt$quiet)){
      cat(paste(filePath,"/random_sequences/shuffled made\n", sep = ""))
    }
  }else{
    if(is.null(opt$quiet)){

      cat(paste(filePath,"/random_sequences/shuffled exists\n", sep = ""))
    }
  }



sraLength <- nrow(sra)
if(sraLength < 500){sraLength <- 500}
if(sraLength > 2000){sraLength <- 2000}

startsRandom <- round(runif(sraLength, min = 1, max = max(fasta$right)))
lengthsRandom <- round(rnorm(sraLength, mean=meanLength, sd=sdLength))
strandRandom <- sample(c("-", "+"), sraLength, replace = T, prob = c( propStrand["-"],  propStrand["+"]))


random <- data.frame(start= startsRandom, end = startsRandom + lengthsRandom, strand = strandRandom)
random <- random %>% arrange(start) %>%
  mutate(tmp.end = ifelse(start > end, start, end),
         tmp.start = ifelse(start > end, end, start)) %>%
  mutate(start = tmp.start,
         end = tmp.end) %>%
  filter(start < end)
}else{

  if (file.exists(paste(filePath,"/random_sequences/shifted", sep = "")) == F){
    dir.create(file.path(filePath,"/random_sequences/shifted"))
    if(is.null(opt$quiet)){
      cat(paste(filePath,"/random_sequences/shifted made\n", sep = ""))
    }
  }else{
    if(is.null(opt$quiet)){

      cat(paste(filePath,"/random_sequences/shifted exists\n", sep = ""))
    }
  }

  sraLength <- nrow(sra)
  if(sraLength < 50){sraLength <- 50}
  if(sraLength > 2000){sraLength <- 2000}

  shifts <- round(runif(1, min = -max(fasta$right)/2, max = max(fasta$right)/2))

  random <- sra %>% mutate(start = (start + shifts[1]) %% max(fasta$right),
                           end = (end + shifts[1])%% max(fasta$right)) %>%
    select(start,end, strand)

  
random <- random %>%
  filter(start < end)
}


ncRNA <- sra%>%
  mutate(feature = ifelse(new_feature == T, "sra_novel", "sra_known")) %>% select(start, end, strand, feature)

ncRNA <- ncRNA %>% arrange(start) %>%
  mutate(tmp.end = ifelse(start > end, start, end),
         tmp.start = ifelse(start > end, end, start)) %>%
  mutate(start = tmp.start,
         end = tmp.end) %>%
  filter(start < end)


cds <- gff %>% filter(feature == "CDS")  %>% select(start, end, strand, feature)

cds <- cds %>% arrange(start) %>%
  mutate(tmp.end = ifelse(start > end, start, end),
         tmp.start = ifelse(start > end, end, start)) %>%
  mutate(start = tmp.start,
         end = tmp.end) %>%
  filter(start < end)

features <- ncRNA %>% bind_rows(cds) %>%
  arrange(start)




query <- IRanges(random$start, random$end)
subject <- IRanges(features$start, features$end)

x <- findOverlaps(query = query, subject = subject)



randomPositions <- x@from
featuresPositions <- x@to

random <- random %>% mutate(random_row = row_number()) %>%
  select(start, end, strand, random_row) %>%
  dplyr::rename(startR = start,
                endR = end,
                strandR = strand)

features <- features %>% mutate(features_row = row_number())%>%
  select(start, end, strand, features_row, feature) %>%
  dplyr::rename(startF = start,
                endF = end,
                strandF = strand)

overlaps <- data.frame(random_row = randomPositions,
                       features_row = featuresPositions)
overlaps <- overlaps %>%
  full_join(random, by  = "random_row")%>%
  left_join(features, by  = "features_row")

overlaps <- overlaps %>% mutate(sense = ifelse(strandR == strandF, "sense", "antisense"),
                                sense = ifelse(feature == "sra_novel", "sense", sense),
                                sense = ifelse(is.na(sense), "sense", sense))

random <- overlaps %>% group_by(startR) %>% summarise(feature_list = paste(unique(feature), collapse = ","))

random <- random %>% mutate(type = ifelse(grepl("CDS", feature_list), "CDS", ifelse(grepl("known", feature_list), "sra_known", ifelse(grepl("novel", feature_list), "sra_novel", "intergenic"))))

random <- random %>% select(startR, type)

overlaps <- random %>%
  left_join(overlaps, by = "startR")




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

seqSense <- sequences %>% filter(sense == "sense", type != "intergenic") %>% mutate(id = paste(">", opt$name, "_random_", row_number(), sep = "")) %>%
  select(id, sequence)

seqAntisense <- sequences %>% filter(sense == "antisense", type != "intergenic") %>% mutate(id = paste(">", opt$name, "_random_", row_number(), sep = "")) %>%
  select(id, sequence)

sequences <- sequences %>% mutate(id = paste(">", opt$name, "_random_", row_number(), sep = "")) %>%
  select(id, sequence)

cat(paste("Writing the output to ", filePath, "/", opt$out_name, "_random_sequences.txt\n", sep = ""))
write.table(x = sequences, file = paste(filePath, "/", opt$out_name, "_random_sequences.fna", sep = ""), row.names = F, col.names = F, quote = F, sep = "\n")

write.table(x = seqSense , file = paste(filePath, "/", opt$out_name, "_random_sequences_sense.fna", sep = ""), row.names = F, col.names = F, quote = F, sep = "\n")
write.table(x = seqAntisense , file = paste(filePath, "/", opt$out_name, "_random_sequences_antisense.fna", sep = ""), row.names = F, col.names = F, quote = F, sep = "\n")


overlaps <- overlaps %>% 
  mutate(strandR = "+", sense = "sense") %>% 
  select(startR, endR, strandR, sense, type) %>% unique()

cat(paste("Writing the output to ", filePath, "/", opt$out_name, "_random_data.txt\n", sep = ""))
write.table(x = overlaps, file = paste(filePath, "/", opt$out_name, "_random_data.txt", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")

#
