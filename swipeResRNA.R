#!/usr/bin/env Rscript
suppressMessages(library('getopt'))


# getopts -----------------------------------------------------------------


spec = matrix(c(
  'file_name', 'f', 1, "character",
  'help' , 'h', 0, "logical",
  'quiet' , 'q', 0, "logical",
  'file_path', 'p', 2, "character"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat("swipeResRNA.R version 1.0\n")
  cat(" \n")
  cat("Use swipeResRNA.R <options> -f <file>\n")
  cat(" \n")
  cat("Options:\n")
  cat("  -f <file> The gff files\n")
  cat("  -q <quiet> Do not print any updates\n")
  cat("  -p <file path> The location of the other files and the output file\n")
  q(status=1)
}




# packages ----------------------------------------------------------------

if ( is.null(opt$file_name) ) {
  cat("Error: -f <file> is required.\n")
  q(status=1)
}


suppressMessages(library(tidyverse))
suppressMessages(library(comparativeSRA))

if ( is.null(opt$file_path ) ) { opt$file_path = "./" }
file_path <- opt$file_path
file_name <- opt$file_name


test <- F
if(test){
file_path <- "~/phd/RNASeq/srna_all/version_8/known_c/multiple_seqs/"
file_name <- "sra_enterics-serratia_known_version_8c_1.fasta"
}

#print(file_path)
#print(file_name)


swipe <- read.table(paste(file_path, file_name, ".swipe", sep = ""), comment.char = "#", fill = T, sep = "\t", header = F, quote = "", as.is = T)
lookup <- read.table(paste(file_path, file_name, ".lookup", sep = ""), comment.char = "#", fill = T, sep = "\t", header = F, quote = "", as.is = T)
fasta <- read.table(paste(file_path, file_name, sep = ""))



seqs <- fasta %>% filter(grepl(pattern = ">", V1) == F)
names <- fasta %>% filter(grepl(pattern = ">", V1))

fasta <- names %>% bind_cols(seqs)
fasta$Query.id <- str_remove_all(string = fasta$V1, pattern = ">")
fasta$seqs <- fasta$V11
fasta <- fasta %>% select(Query.id, seqs)
if(nrow(fasta) == 1){
  cat(paste("Single sequence for", file_name, "\n"))
  q(status=1)
}


colnames(swipe) <- c("Query.id", "Subject.id", "% identity", "alignment length", "mismatches", "gap openings", "q. start", "q. end", "s. start", "s. end","e.value", "score")
lookup <- lookup %>% filter(grepl(pattern = ">", V1)) %>% separate(col = V1, into = c("Subject.id", "rna.id"), sep = " ", remove = T, extra = "merge")
lookup$Subject.id <- str_remove_all(string = lookup$Subject.id, pattern = ">")



swipe <- swipe %>% 
  left_join(lookup, by = "Subject.id") %>% 
  left_join(fasta, by = "Query.id")

swipe <- swipe %>%  mutate(Subject.id = rna.id)

swipe <- swipe %>% filter(Query.id != Subject.id, e.value < 1)
#print(head(swipe))

fasta$new_seq <- NA

#print(head(fasta))


lengths <- swipe %>% group_by(Subject.id, Query.id) %>% summarise(length.total = sum(`alignment length`))
i <- 2
queryList <- c()

if(nrow(lengths) == 0){
  cat(paste("No conserved seqs for", file_name, "\n"))
  q(status=1)
}


for(i in 1:nrow(lengths)){
  query <- lengths$Query.id[i] 
  if(query %in% queryList == F){
  if(lengths$length.total[i] > 35){
    subject <- lengths$Subject.id[i] 
    start <- min(swipe$`q. start`[swipe$Query.id == query])
    stop <- max(swipe$`q. end`[swipe$Query.id == query])
    
    if(start <= 10){
      start <- 1
    }else{
      start <- start - 10
    }
    

    if(stop >= nchar(as.character(fasta$seqs[fasta$Query.id == query])) - 10){
      stop <- nchar(as.character(fasta$seqs[fasta$Query.id == query]))
    }else{
      stop <- stop + 10
    }
    
    
    
    seq <- substr(fasta$seqs[fasta$Query.id == query], start = start, stop = stop)
    
    fasta$new_seq[fasta$Query.id == query] <- seq
    
    mat <- data.frame(a = c(paste(">", query, sep = ""), seq))
    write.table(mat, file = paste(file_path, "new/", file_name, sep = ""),
                row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  }
  queryList <- c(queryList, query)
  }
}


