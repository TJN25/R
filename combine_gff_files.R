#!/usr/bin/env Rscript
suppressMessages(library('getopt'))


# getopts -----------------------------------------------------------------


spec = matrix(c(
  'sra', 'f', 1, "character",
  'gff', 'g', 1, 'character',
  'help' , 'h', 0, "logical",
  'stranded' , 's', 0, "logical",
  'quiet' , 'q', 0, "logical",
  'file_path', 'p', 2, "character",
  'out_name', 'o', 2, "character",
  'random_data', 'r', 1, "character"
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
  cat("  -r <random data> The file to remove CDS regions from\n")
  cat("  -q <quiet> Do not print any updates\n")
  cat("  -p <file path> The location of the other files and the output file\n")
  cat("  -o <output file name> The name of the output file. Do not inclue the gff file extension. The default is the same as the sra input\n")
  q(status=1)
}

if ( is.null(opt$sra) ) {
  cat("Error: -f <files> is required.\n")
  q(status=1)
}

if ( is.null(opt$out_name) ) {
  cat("Error: -o <output file name> is required.\n")
  q(status=1)
}


# packages ----------------------------------------------------------------


suppressMessages(library(tidyverse))
suppressMessages(library(comparativeSRA))

# defining variables ------------------------------------------------------


if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$random_data ) ) { opt$random_data = "" }
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



# gffDat  <- tryCatch({
#   suppressWarnings(gffDat <- read.table(paste(filePath, "/", gffName, ".gff", sep = ""), sep = "\t", fill = T, comment.char = "#", quote = ""))
#   gffDat
# }, error =  function(e) {
#   cat(paste("Error: ", opt$file_path, "/", opt$gff, ".gff not found.\n", sep = ""))
#   q(status=1)
# })

#####

file_path <- opt$file_path
files <- list.files(paste(file_path, opt$sra, sep = "/"), pattern = ".gff$")
# import data -------------------------------------------------------------


#print(files)
dat <- data.frame(sequence = as.character("0"), source = as.character("0"), feature = as.character("0"),
                  start = as.integer("0"), end = as.integer("0"), score = as.character("0"),
                  strand = as.character("0"), phase = as.character("0"), Atrribute = as.character("0"), file_name = as.character("start_row"), stringsAsFactors = F)
i <- 2
for(i in 1:length(files)){
  tmp  <- tryCatch({
    suppressWarnings(tmp <- read.table(paste(file_path, opt$sra, files[i], sep = "/"), comment.char = "#", quote = "", sep = "\t", as.is = T))
  }, error =  function(e) {
    cat(paste("Error: ", "row ", i, ", ", file_path, "/", opt$sra, "/", files[i], " cannot be opened.\n", sep = ""))
    cat(paste(e, "\n"))
  })

  if(class(tmp) == "NULL"){
    next
  }

  if(ncol(tmp) != 9){
    cat(paste("Error: ", "row ", i, ", ", file_path, "/", opt$sra, "/", files[i], " contains ", ncol(tmp), " columns.\n", sep = ""))
    next
  }

  colnames(tmp) <- c("sequence", "source", "feature", "start", "end", "score", "strand", "phase", "Atrribute")

  tmp <- tmp%>%mutate(file_name = files[i])%>%mutate(score = as.character(score))

  if(files[i] == opt$random_data){
    tmp <- tmp%>%
  filter(feature != "CDS", feature != "gene", feature != "pseudogene", feature != "exon", feature != "region")
  }else{
  
  dat <- dat%>%bind_rows(tmp)
}
}
if(!is.null(opt$random_data)){
   ncRNAgff <- dat%>%
     filter(feature != "gene", feature != "pseudogene", feature != "exon", feature != "region")
}else{
ncRNAgff <- dat%>%
  filter(feature != "CDS", feature != "gene", feature != "pseudogene", feature != "exon", feature != "region")
}

# main section  -------------------------------------------------------------------


ncRNAgff <- ncRNAgff%>%arrange(start) %>% filter((end - start) > 0)# %>% arrange(strand)


mergedDat <- data.frame(sequence = as.character("0"), feature = as.character("0"),
                        start = as.integer("0"), end = as.integer("0"),
                        strand = as.character("0"), file_names = as.character("start_row"),
                        row_numbers = as.character("0"), prop_overlap = as.numeric(0), new_feature = F,
                        number_of_rnaseq_files = as.integer("0"),
                        score = as.character("0"),
                        stringsAsFactors = F)

##loop through the combined gff files and combine features that overlap
i <- 3
current_feature <- F #is there a current feature being written?
new_feature <- T

for(i in 1:(nrow(ncRNAgff))){
  ##check if the feature is already known
  if(ncRNAgff$source[i] != "sraAlignedncRNAExpression"){
    new_feature <- F
  }

  ##if there is no current feature then set a new start value
  if(current_feature == F){
  start_val <- ncRNAgff$start[i]
  start_i <- i
  end_val <- ncRNAgff$end[i]
  }



  ##set the new end value
  if(ncRNAgff$end[i] > end_val){
  end_val <- ncRNAgff$end[i]
  }

  if(i == nrow(ncRNAgff)){
    
    ##check if the subsequent feature was contained within the first feature
    if(ncRNAgff$end[start_i] < end_val){
      prop_val <- (ncRNAgff$end[start_i] - ncRNAgff$start[i])/(end_val - start_val)
    }else{
      prop_val <- 1
    }
    
    tmp <- data.frame(sequence = ncRNAgff$sequence[i],
                      feature = ncRNAgff$feature[i],
                      start = start_val, end = end_val,
                      strand = ncRNAgff$strand[i],
                      file_names = paste(ncRNAgff$file_name[start_i:i], collapse = ","),
                      row_numbers = paste(c(start_i:i), collapse = ","),
                      prop_overlap = prop_val,
                      new_feature = new_feature,
                      number_of_rnaseq_files = length(start_i:i),
                      score = as.character(ncRNAgff$score[i]),
                      stringsAsFactors = F)
    mergedDat <- mergedDat%>%bind_rows(tmp)
    current_feature <- F
    new_feature <- T
  }else{
    
    
  ##check if the cuurent end value overlaps with the next starting value and update the end value if it does
  if(end_val > ncRNAgff$start[i + 1]){
    end_val <- ncRNAgff$end[i + 1]
    current_feature <- T
  }else{

    ##check if the subsequent feature was contained within the first feature
    if(ncRNAgff$end[start_i] < end_val){
    prop_val <- (ncRNAgff$end[start_i] - ncRNAgff$start[i])/(end_val - start_val)
    }else{
      prop_val <- 1
    }

    tmp <- data.frame(sequence = ncRNAgff$sequence[i],
                      feature = ncRNAgff$feature[i],
                      start = start_val, end = end_val,
                      strand = ncRNAgff$strand[i],
                      file_names = paste(ncRNAgff$file_name[start_i:i], collapse = ","),
                      row_numbers = paste(c(start_i:i), collapse = ","),
                      prop_overlap = prop_val,
                      new_feature = new_feature,
                      number_of_rnaseq_files = length(start_i:i),
                      score = as.character(ncRNAgff$score[i]),
                      stringsAsFactors = F)
    mergedDat <- mergedDat%>%bind_rows(tmp)
    current_feature <- F
    new_feature <- T
  }
  }
}





mergedDat <- mergedDat%>%filter(number_of_rnaseq_files > 0, file_names != "start_row")

# if(!is.null(opt$random_data)){
#   mergedDat <- mergedDat %>% filter(file_names != opt$gff)
# }

mergedDat <- mergedDat %>% mutate(id =  paste(opt$out_name, row_number(), sep = "_"))

cat(paste("Writing the output to ", file_path, "/", opt$out_name, "_new_calls.txt\n", sep = ""))
write.table(x = mergedDat, file = paste(file_path, "/", opt$out_name, "_new_calls.txt", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")



