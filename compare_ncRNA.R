#!/usr/bin/env Rscript

##Seems to all be working

# functions ---------------------------------------------------------------


mergeSRATest <- function(ncRNAgff, gff1, gff2, time.it = T, quiet = F, filenum1 = "1", filenum2 = "2", print_log = F, align = T){
  error_message <- "Either gff1 and gff2 or ncRNAgff are needed:\n"
  stop_val <- 0
  log_file = ""
  if(missing(gff1)){
    error_message <- paste(error_message, "\tArugment gff1 missing\n", sep = "")
    stop_val <- stop_val + 1
  }else{
    gff1 <- gff1%>%mutate(filenum = filenum1)
  }

  if(missing(gff2)){
    error_message <- paste(error_message, "\tArugment gff2 missing\n", sep = "")
    stop_val <- stop_val + 1
  }else{
    gff2b <- gff2%>%mutate(filenum = filenum2)
  }


  if(stop_val > 0){
    if(missing(ncRNAgff)){
      error_message <- paste(error_message, "\tArugment ncRNAgff missing\n", sep = "")
      stop(error_message)
    }
  }
  if(missing(ncRNAgff)){
    ncRNAgff <- gff2b%>%bind_rows(gff1)%>%unique()
  }else{
    if(quiet == F){
      cat("Using the ncRNAgff dataframe:\n")
    }
  }

  ptm <- proc.time()

  ncRNAgff <- ncRNAgff%>%arrange(start)%>%arrange(strand)


  mergedDat <- data.frame(sequence = as.character("0"), feature = as.character("0"),
                          start = as.integer("0"), end = as.integer("0"),
                          strand = as.character("0"), file_names = as.character("start_row"),
                          row_numbers = as.character("0"), prop_overlap = as.numeric(0), feature_match = F,
                          number_of_features = as.integer("0"),
                          score = as.character("0"),
                          new_feature = F,
                          number_of_rnaseq_files = as.integer("0"),
                          id1 = as.character("0"),
                          id2 = as.character("0"),
                          set_val_1 = as.character("0"),
                          set_val_2 = as.character("0"),
                          stringsAsFactors = F)

  ##loop through the combined gff files and combine features that overlap
  i <- 2
  current_feature <- F #is there a current feature being written?
  new_feature <- F
  for(i in 1:(nrow(ncRNAgff) - 1)){

    if(quiet ==F){
      printRemaining(i <- i, length = nrow(ncRNAgff) - 1, increment = 5)
    }


    if(ncRNAgff$start[i] < 0 && align == T){
      start_val <- ncRNAgff$start[i]
      start_i <- i
      end_val <- ncRNAgff$end[i]

      feature_matched <- ifelse(length(unique(ncRNAgff[start_i:i, 14])) > 1, T, F)
      prop_val <- 1
      idRows <- ncRNAgff[start_i:i,]
      id1_val <- idRows[idRows[,16] == filenum1,14]
      id2_val <- idRows[idRows[,16] == filenum2,14]
      if(is_empty(id1_val)){

        id_vals <- unlist(strsplit(filenum1, "-"))

        id_vals <- paste(id_vals, "_0", sep = "")

        id1_val <- paste(id_vals, collapse = "-")
      }
      if(is_empty(id2_val)){

        id_vals <- unlist(strsplit(filenum2, "-"))

        id_vals <- paste(id_vals, "_0", sep = "")

        id2_val <- paste(id_vals, collapse = "-")
      }


      set_val1 <- ifelse(filenum1 == ncRNAgff$file_id[i], ncRNAgff$set_val[i], 0)

      set_val2 <- ifelse(filenum2 == ncRNAgff$file_id[i], ncRNAgff$set_val[i], 0)




      if(is_empty(set_val1)){
        set_val1 <- "0"
      }
      if(is_empty(set_val2)){
        set_val2 <- "0"
      }




      #id1_val <- "1"
      #id2_val <- "2"


      tmp <- data.frame(sequence = ncRNAgff[i,1],
                        feature = ncRNAgff[i,2],
                        start = start_val, end = end_val,
                        strand = ncRNAgff[i,5],
                        file_names = paste(unique(ncRNAgff[start_i:i, 16]), collapse = ","),
                        row_numbers = paste(c(start_i:i), collapse = ","),
                        prop_overlap = prop_val,
                        feature_match = feature_matched,
                        number_of_features = sum(ncRNAgff$number_of_features[start_i:i]),
                        score = as.character(ncRNAgff[i,11]),
                        new_feature = !(F %in% ncRNAgff[start_i:i, 12]),
                        number_of_rnaseq_files = sum(as.integer(ncRNAgff[start_i:i, 10])),
                        id1 = id1_val,
                        id2 = id2_val,
                        set_val_1 = as.character(set_val1),
                        set_val_2 = as.character(set_val2),
                        stringsAsFactors = F)
      mergedDat <- mergedDat%>%bind_rows(tmp)


      next
    }



    ##if there is no current feature then set a new start value
    if(current_feature == F){
      start_val <- ncRNAgff[i,3]
      start_i <- i
      end_val <- ncRNAgff[i,4]
    }



    ##set the new end value
    if(ncRNAgff[i, 4] >= end_val){
      end_val <- ncRNAgff[i,4]
    }

    if(ncRNAgff[i, 3] <= start_val){
      start_val <- ncRNAgff[i,3]
    }
#& ncRNAgff[i,5] == ncRNAgff[i+1, 5] strand has been removed
    ##check if the current end value overlaps with the next starting value and update the end value if it does
    if(end_val >= ncRNAgff[i + 1, 3] ){
      if(ncRNAgff[i + 1, 4] >= end_val){
        end_val <- ncRNAgff[i + 1,4]
      }

      if(ncRNAgff[i + 1, 3] <= start_val){
        start_val <- ncRNAgff[i + 1,3]
      }
      current_feature <- T
    }else{

      ##check if the subsequent feature was contained within the first feature
      if(ncRNAgff[start_i, 4] <= end_val){
        prop_val <- (ncRNAgff[start_i, 4] - ncRNAgff[i, 3])/(end_val - start_val)
      }else{
        prop_val <- 1
      }
      feature_matched <- ifelse(length(unique(ncRNAgff[start_i:i, 14])) > 1, T, F)



        idRows <- ncRNAgff[start_i:i,]
        id1_val <- paste(idRows$id[idRows$file_id == filenum1], collapse = "-")
        id2_val <- paste(idRows$id[idRows$file_id == filenum2], collapse = "-")

        if(id1_val == ""){
          id_vals <- unlist(strsplit(filenum1, "-"))

          id_vals <- paste(id_vals, "_0", sep = "")

          id1_val <- paste(id_vals, collapse = "-")
        }

        if(id2_val == ""){
          id_vals <- unlist(strsplit(filenum2, "-"))

          id_vals <- paste(id_vals, "_0", sep = "")

          id2_val <- paste(id_vals, collapse = "-")
          }



        set_val_1 <- idRows$set_val[idRows$file_id == filenum1]
        set_val_2 <- idRows$set_val[idRows$file_id == filenum2]

        set_val_1 <- interset_all(set_val_1)
        set_val_2 <- interset_all(set_val_2)


        if(set_val_1 == ""){
          set_val_1 <- "0"
        }
        if(set_val_2 == ""){
          set_val_2 <- "0"
        }




        tmp <- data.frame(sequence = ncRNAgff[i,1],
                          feature = ncRNAgff[i,2],
                          start = start_val, end = end_val,
                          strand = ncRNAgff[i,5],
                          file_names = paste(unique(ncRNAgff[start_i:i, 16]), collapse = ","),
                          row_numbers = paste(c(start_i:i), collapse = ","),
                          prop_overlap = prop_val,
                          feature_match = feature_matched,
                          number_of_features = sum(ncRNAgff$number_of_features[start_i:i]),
                          score = as.character(ncRNAgff[i,11]),
                          new_feature = !(F %in% ncRNAgff[start_i:i, 12]),
                          number_of_rnaseq_files = sum(as.integer(ncRNAgff[start_i:i, 10])),
                          id1 = id1_val,
                          id2 = id2_val,
                          set_val_1 = set_val_1,
                          set_val_2 = set_val_2,
                          stringsAsFactors = F)
        mergedDat <- mergedDat%>%bind_rows(tmp)

      current_feature <- F
      new_feature <- F
    }
  }

  runningTime <- proc.time() - ptm
  if(time.it){
    if(quiet == F){
      printRunningTime(runningTime = runningTime)
    }
  }

  if(print_log){
    cat(log_file)
  }

  mergedDat <- mergedDat%>%filter(number_of_features > 0, file_names != "start_row")
  return(mergedDat)

}

reorderGFFTest <- function(ref, gff, time.it = T, quiet = F, reference.genome = F){
  ref2 <- ref%>%arrange(start.b)

  ptm <- proc.time()

  if(reference.genome == F){
    for(i in 1:nrow(gff)){
      if(quiet == F){
        printRemaining(i = i, length = nrow(gff), increment = 5)
      }
      start_val <- gff[i, 3]
      end_val <- gff[i, 4]
      for(j in 1:nrow(ref2)){
        if(start_val >= ref2[j,4] && start_val <= ref2[j,5]){
          gff[i,3] <- gff[i,3] + ref2[j, 7]
          gff[i,4] <- gff[i,4] + ref2[j, 7]
          gff[i, 17] <- T
        }
      }

    }
  }else{
    for(i in 1:nrow(gff)){
      if(quiet == F){
        printRemaining(i = i, length = nrow(gff), increment = 5)
      }
      start_val <- gff[i, 3]
      end_val <- gff[i, 4]
      for(j in 1:nrow(ref2)){
        if(start_val >= ref2[j,1] && start_val <= ref2[j,2]){
          gff[i, 17] <- T
        }
      }

    }
  }

  runningTime <- proc.time() - ptm
  if(time.it){
    if(quiet == F){
      printRunningTime(runningTime = runningTime)
    }
  }
  gff <- gff%>%filter(changed == T)
  return(gff)
}

alignAndCombineTest <- function(reference, gff1, gff2, time.it = T, quiet = F, filenum1 = "1", filenum2 = "2", seqA = 1, seqB = 2){

  referenceEsch1Serr1 <- read.table(reference, header = T, as.is = T)
  referenceEsch1Serr1Built <- buildReferenceLookup(reference = referenceEsch1Serr1,
                                                   as.numeric(seqA), seqB = as.numeric(seqB),
                                                   collapse.alignment = T,
                                                   quiet = quiet)

  esch1 <- gff1
  serr1 <- gff2

  esch1 <- esch1%>%mutate(changed = F)#%>%
    #dplyr:: mutate(id = paste(filenum1, row_number(), sep = "_"))
  serr1 <- serr1%>%mutate(changed = F)#%>%
    #dplyr::mutate(id = paste(filenum2, row_number(), sep = "_"))



  serr1b <- reorderGFFTest(ref = referenceEsch1Serr1Built, gff = serr1, time.it = time.it, quiet = quiet)
  esch1b <- reorderGFFTest(ref = referenceEsch1Serr1Built, gff = esch1, reference.genome = T, time.it = time.it, quiet = quiet)

  serr1b <- serr1b%>%mutate(filenum = filenum2)
  esch1b <- esch1b%>%mutate(filenum = filenum1)
  ncRNAgff <- esch1b%>%bind_rows(serr1b)


  serr1 <- serr1%>%mutate(filenum = filenum2)
  esch1 <- esch1%>%mutate(filenum = filenum1)
  ncRNAgffUnchanged <- esch1%>%bind_rows(serr1)
  tmp1 <- ncRNAgff %>% select(id) %>% mutate(found = T)
  tmp2 <- ncRNAgffUnchanged %>% select(id)
  tmp1 <- tmp1 %>% full_join(tmp2, by = "id") %>% filter(is.na(found)) %>% mutate(found = F)

  otherncRNA <- ncRNAgffUnchanged %>% full_join(tmp1, by = "id") %>% filter(found == F) %>% select(-found) %>%
    mutate(start = -start, end = -end)

  ncRNAgff <- ncRNAgff %>% bind_rows(otherncRNA)

  return(ncRNAgff)

}

interset_all <- function(set_val){
  setList <- list()
  for(j in 1:length(set_val)){
    value <- set_val[j]
    setList[j] <- strsplit(as.character(value), "-")
  }

  set_val <- paste(Reduce(intersect, setList), collapse = "-")
  return(set_val)
}




# getopts -----------------------------------------------------------------



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

if ( !is.null(opt$clean) ) {
  clean_data <- T
}else{
  clean_data <- F
}



if ( is.null(opt$gff2) ) {
  if(clean_data == F){
  cat("Error: -g <other gff file> is required.\n")
  q(status=1)
  }
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
opt$gff1 <- "~/phd/RNASeq/combined_gff_files/GCA_000017745.1-GCA_000017765.1_merged.gff"
opt$gff2 <- "~/phd/RNASeq/combined_gff_files/GCA_000017745.1-GCA_000017985.1_merged.gff"
opt$alignment <- "~/phd/RNASeq/escherichia/escherichia.backbone"
opt$id1 <- "GCA_000017745.1-GCA_000017765.1"
opt$id2 <- "GCA_000017745.1-GCA_000017985.1"
opt$out_name <- "esch_1-2-3"
}else{
opt$gff1 <- "~/phd/RNASeq/escherichia/GCA_000017745_data/GCA_000017745.1_new_calls.txt"
opt$gff2 <- "~/phd/RNASeq/escherichia/GCA_000017765.1_data/GCA_000017765.1_new_calls.txt"
opt$alignment <- "~/phd/RNASeq/escherichia/escherichia.backbone"
opt$id1 <- "GCA_000017745.1"
opt$id2 <- "GCA_000017765.1"
opt$out_name <- "escherichia_1-2_TEST"
}
}



# Defining variables ------------------------------------------------------



if ( is.null(opt$file_path ) ) { opt$file_path = "." }
if ( is.null(opt$out_name ) ) { opt$out_name = paste(opt$x, opt$y, sep = "-") }
if ( is.null(opt$x ) ) {  opt$x = "1" }
if ( is.null(opt$y ) ) { opt$y = "2" }
if ( is.null(opt$s ) ) {  opt$s = "1" }
if ( is.null(opt$t ) ) { opt$t= "2" }

filePath <- opt$file_path

# Main section ------------------------------------------------------------
if(initial_data == T){
  cat("Analysing initial calls (from *_new_calls.gff)\n")
  ncRNAgff <- alignAndCombine(reference = opt$alignment,
                                      gff1 = opt$gff1,
                                      gff2 = opt$gff2,
                                      filenum1 = opt$id1,
                                      filenum2 = opt$id2,
                                      seqA = opt$s,
                                      seqB = opt$t)

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



   # opt$gff1 <- "~/phd/RNASeq/combined_gff_files/esch_1-2-3-4-5_merged.gff"
   # opt$gff2 <- "~/phd/RNASeq/combined_gff_files/GCA_000017745.1-GCA_001559675.1_merged.gff"
   # align <- F

 
  
  
  gff1Dat <- read.table(opt$gff1, sep = "\t", header = T, as.is = T)
  gff2Dat <- read.table(opt$gff2, sep = "\t", header = T, as.is = T)


  gff1Working <- gff1Dat %>% mutate(row_numbers = as.character(row_numbers))
  gff2Working <- gff2Dat %>% mutate(row_numbers = as.character(row_numbers))
  

  filenum1 <- gff1Working$file_id[1]
  filenum2 <- gff2Working$file_id[1]
  
  print(opt$gff1)
  print(opt$gff2)  
  print(filenum1)
  print(filenum2)
  print(filePath)
  print(opt$out_name)
  print(align)
  print(initial_data)


  if(align == T){
   ncRNAgff <- alignAndCombineTest(reference = opt$alignment,
                               gff1 = gff1Working,
                               gff2 = gff2Working,
                               filenum1 = filenum1,
                               filenum2 = filenum2,
                               seqA = 1,
                               seqB = 2)

   ncRNAgff <- ncRNAgff%>%select(-changed, -filenum)%>%unique()

  }else{
    ncRNAgff <- gff1Working%>%bind_rows(gff2Working)
    ncRNAgff[is.na(ncRNAgff)] <- 0
}
print(nrow(ncRNAgff))
  mergedData <- mergeSRA(ncRNAgff = ncRNAgff,
                         filenum1 = filenum1,
                         filenum2 = filenum2,
                         align = align, 
                         initial_data = F)
  print(nrow(mergedData))
  

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

}

cat(paste("Writing the output to ", filePath, "/", opt$out_name, "_merged.gff\n", sep = ""))
write.table(x = mergedData, file = paste(filePath, "/", opt$out_name, "_merged.gff", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")


