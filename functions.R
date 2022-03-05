###Should list the required packages here

#take a named vector and convert to a data frame with the names turned into a column
vector_to_dataframe <- function(vec, new_col = 'V2'){
  vec <- as.data.frame(vec) 
  vec[[new_col]] <- row.names(vec)
  rownames(vec) <- NULL
  return(vec)
}

#select and order columns of a data frame using a vector of column names
select_columns_by_list <- function(dat, vec){
  dat <- dat[, colnames(dat) %in% vec]
  dat <- dat[, match(vec, colnames(dat))]
  return(dat)
}

#set first letter to capital
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#take an input, and calculate the cumulative counts of rows sorted by the target column
cumulativeCounts <- function(dists, smooth = T, target_column = 'max_dist'){
  groups <- unique(dists$group)
  for(i in groups){
    dat <- dists %>% filter(group == i)
    dat <- dat %>% mutate(count = 1) %>%
      arrange(-.data[[target_column]]) %>% group_by(group) %>%
      mutate(cumulativeCount = cumsum(count)) %>% ungroup() %>%
      group_by(group, .data[[target_column]]) %>% summarise(cumulative_prop = max(cumulativeCount)/ nrow(dat))
    if(smooth){
      dat <- as.data.frame(spline(x = dat[[target_column]],y =  dat$cumulative_prop))
    }
    dat <- dat %>% ungroup() %>% mutate(group = i)
    if(exists('combinedDat')){
      combinedDat <- combinedDat %>% bind_rows(dat)
    }else{
      combinedDat <- dat
    }
  }
  combinedDat$group <- factor(combinedDat$group, levels = c('low', 'Positive Control', 'high'))
  return(combinedDat)
}

#plot the output of the cumulative counts and run a K-S test
cumulativeDistribution <- function(dat, run.ks.test = F, alternative_pred, alternative_pc, alternative = 'two.sided', show.legend = F, target_column = 'max_dist'){
  if(missing('alternative_pred')){alternative_pred <- alternative}
  if(missing('alternative_pc')){alternative_pc <- alternative}
  distsCumulativeCount <- cumulativeCounts(dists = dat, smooth = F, target_column = target_column)
  distsCumulativeCount <- distsCumulativeCount %>% filter(group != "Predicted Known")
  
  p <- ggplot() +
    geom_line(data = distsCumulativeCount, aes(x= .data[[target_column]], y = cumulative_prop, group = group, colour = group), size = 1, show.legend = show.legend) + theme_classic()
  
  if(run.ks.test == T){
    pos <- dat %>% filter(group == "Positive Control")
    neg <- dat %>% filter(group == "Negative Control")
    pred <- dat %>% filter(group == "Predicted")
    res <- ks.test(x = pred[[target_column]], y = neg[[target_column]], alternative = alternative_pred)
    print(res)
    res <- ks.test(x = pos[[target_column]], y = neg[[target_column]], alternative = alternative_pc)
    print(res)
  }
  return(p)
}

#set up the data sets for later merging
readDepthsSetup <- function(file_path){
  readsColName <- c("mean.val", "max.val", "counts.above.threshold", "ID", "genus")
  reads.dat <- read.table(file_path, as.is = T, stringsAsFactors = F)
  colnames(reads.dat) <- readsColName
  reads.max <- reads.dat  %>% group_by(ID) %>% summarise(read.max.score  =  sum(max.val))
  reads.mean <- reads.dat %>%
    group_by(ID) %>%
    summarise(reads.mean.score  = sum(mean.val))
  reads.count <- reads.dat %>%
    group_by(ID) %>%
    filter(max.val > 0) %>%
    summarise(max.val = max(max.val), read.counts  = n()) %>% 
    select(-max.val)
  reads.dat <- reads.max %>% full_join(reads.mean, by = "ID")
  reads.dat <- reads.dat %>% full_join(reads.count, by = "ID")
  return(reads.dat)
}

rscapeCovarianceSetup <- function(file_path){
  cov.dat <- read.table(file_path, sep = "\t", comment.char = "#", as.is = T, header = F, fill = T, col.names = c("V1", "left_pos", "right_pos", "score", "e.value", "substitutions", "V2", "power", "ID"))
  cov.dat <- cov.dat %>% select(ID, score, e.value, power)
  cov.mean <- cov.dat %>% group_by(ID) %>% summarise(cov.mean.score = max(score))
  cov.count <- cov.dat %>% group_by(ID) %>% summarise(cov.count = n())
  cov.max <- cov.dat %>% group_by(ID) %>% summarise(cov.min.eval = min(e.value))
  cov.power <- cov.dat %>% group_by(ID) %>% summarise(power = sum(power))
  cov.dat <- cov.mean %>% full_join(cov.max, by = "ID") %>%
    full_join(cov.count, by = "ID") %>%
    full_join(cov.power, by = "ID")
  cov.dat <- cov.dat %>% filter(!is.na(cov.mean.score)) %>% mutate(cov.combined.score = cov.count * cov.mean.score) 
  return(cov.dat)
}

gcSetup <- function(file_path){
  dat <- read.table(file_path)
  colnames(dat) <- c("counts", "letter", "ID")
  
  datTotals <- dat %>% group_by(ID) %>% summarise(total = sum(counts))
  
  dat <- dat %>% filter(letter %in% c("C", "G")) %>% group_by(ID) %>% summarise(gc.count = sum(counts)) %>% left_join(datTotals, by = "ID") %>% 
    mutate(gc.score = (gc.count/total)*100) %>% select(ID, gc.score)
  
  return(dat)
}

alifoldSetup <- function(file_path){
  dat<- read.table(file_path, header = F, comment.char = "#", quote = "", sep = "",   fill = T, as.is = T, col.names = c( "From",      "To",    "Strand",    "Native.MFE",    "Mean.MFE",     "STDV",        "Z", "ID"))
  
  dat <- dat %>% filter(ID != "")
  
  datSD <- dat %>% select(ID, STDV)
  datMean <- dat %>% group_by(ID) %>% summarise(z_mean = mean(as.numeric(Z), na.rm = T))
  datMax <- dat %>% group_by(ID) %>% summarise(z_max = max(as.numeric(Z), na.rm = T))
  dat <- datMean %>% full_join(datMax, by = "ID") %>% full_join(datSD, by = "ID") %>% 
    filter(!is.nan(STDV)) %>% select(-STDV) %>% unique()
  
  return(dat)
  
}

motifSetup <- function(file_path){
  dat <- read.table(file_path, sep = "", comment.char = "#", as.is = T, header = F, fill = T)
  
  colnames(dat) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute", "ID")
  dat <- dat %>% group_by(feature, start, end, strand, ID) %>% summarise(score = max(score))
  
  
  
  datMean <- dat %>% group_by(ID) %>% summarise(motif.mean.score = mean(score))
  datMax <- dat %>% group_by(ID) %>% summarise(motif.max.score = max(score))
  datCount <- dat %>% group_by(ID) %>% summarise(motif_count = n())
  
  dat <- datMean %>% full_join(datMax, by = "ID") %>% 
    full_join(datCount, by = "ID")
  return(dat)
}

#replace the NAs with suitable null values to keep rows for random forest
replaceNAs <- function(dat){
  
  dat$mfe.score[is.na(dat$mfe.score)] <- 0
  dat$gc.score[is.na(dat$gc.score)] <- 50
  dat$distance[is.na(dat$distance)] <- 0
  dat$reads.mean.score[is.na(dat$reads.mean.score)] <- 0
  dat$read.max.score[is.na(dat$read.max.score)] <- 0
  dat$cov.mean.score[is.na(dat$cov.mean.score)] <- 0
  dat$cov.min.eval[is.na(dat$cov.min.eval)] <- 10
  # dat$cov.combined.score[is.na(dat$cov.combined.score)] <- 0
  dat$cov.count[is.na(dat$cov.count)] <- 1
  dat$read.counts[is.na(dat$read.counts)] <- 0
  dat$motif.mean.score[is.na(dat$motif.mean.score)] <- 0
  dat$motif.max.score[is.na(dat$motif.max.score)] <- 0
  dat$motif_count[is.na(dat$motif_count)] <- 0
  dat$z_mean[is.na(dat$z_mean)] <- 10
  dat$z_max[is.na(dat$z_max)] <- 10
  dat <- dat[dat$z_max != -Inf,]
  dat$alifold.score[is.na(dat$alifold.score)] <- 0
  dat$alifold_cov_score[is.na(dat$alifold_cov_score)] <- 0
  # dat$alifold_cov_score[dat$alifold_cov_score > 0] <- 0
  return(dat)
}

#calculate statistics for a dataset
scoreProbabities <- function(probDat, threshold, target_column = 'max_dist'){
  
  tp <- probDat %>% filter(group == "Positive Control", .data[[target_column]] > threshold) %>% nrow()
  fp <- probDat %>% filter(group == "Negative Control", .data[[target_column]] > threshold) %>% nrow()
  pos <- probDat %>% filter(group == "Positive Control") %>% nrow()
  pred_pos <- probDat %>% filter(group == "Predicted", .data[[target_column]] > threshold) %>% nrow()
  pred_neg <- probDat %>% filter(group == "Predicted", .data[[target_column]] <= threshold) %>% nrow()
  pred <- probDat %>% filter(group == "Predicted") %>% nrow()
  
  
  tn <- probDat %>% filter(group == "Negative Control", .data[[target_column]] <= threshold) %>% nrow()
  fn <- probDat %>% filter(group == "Positive Control", .data[[target_column]] <= threshold) %>% nrow()
  neg <- probDat %>% filter(group == "Negative Control") %>% nrow()
  
  sens <- tp/pos ##sensitivity
  spec <- tn/neg ##specificity
  ppv <- tp/(tp+fp) ##positive predictive value
  fpr <- fp/neg
  fnr <- fn/pos
  pred_pct <- pred_pos/pred
  pc_pct <- tp/pos
  nc_pct <- fp/neg
  
  res <- probDat %>% group_by(group) %>% summarise(positive_res = round(mean(.data[[target_column]] > threshold), 3)) %>% mutate(negative_res = 1 - positive_res)
  
  outList <- list()
  outList$sensitivity <- sens
  outList$specificity <- spec
  outList$ppv <- ppv
  outList$fpr <- fpr
  outList$fnr <- fnr
  outList$tp <- tp
  outList$fp <- fp
  outList$pred_pos <- pred_pos
  outList$pred_neg <- pred_neg
  outList$tn <- tn
  outList$fn <- fn
  outList$neg <- neg
  outList$pred <- pred
  outList$pred_pct <- pred_pct
  outList$pc_pct <- pc_pct
  outList$nc_pct <- nc_pct
  outList$res <- res
  return(outList)
}


printListSubset <- function(inList, vec = NA, startText = NA, getNames = F, round_val = 10){
  if(getNames){
    outNames <- "c('"
    for(value in names(scores)){
      outNames <- paste0(outNames, value, "', '")
    }
    outNames <- paste0(outNames, "')\n")
    cat(outNames)
  }else{
  outputPrint = ''
  if(!is.na(startText)){
    outputPrint = paste0(startText, '\n  ')
  }
  
  for(value in vec){
      if(value %in% names(inList)){
        outputPrint = paste0(outputPrint, value, ': ', round(inList[[value]], round_val), '\n  ')
      }else{
        print(paste0(value, 'not found'))
      }
  }
  cat(outputPrint)
  }
}


allScores <- function(dat, steps, target_column = 'max_dist'){
  valuesDat <- data.frame(threshold = seq((0-steps),max(dat[[target_column]]), by=steps),
                          ppv = NA, 
                          fpr = NA,
                          fnr = NA,
                          pred_pct = NA,
                          pc_pct = NA,
                          pred_pos = NA,
                          sensitivity = NA,
                          specificity = NA
                          )
  counter <- 0
  for(i in seq((0-steps),max(dat[[target_column]]), by=steps)){
    counter <-  counter + 1
    scores <- scoreProbabities(dat, i, target_column)
    valuesDat$ppv[counter] <- scores$ppv 
    valuesDat$fpr[counter] <- scores$fpr 
    valuesDat$fnr[counter] <- scores$fnr 
    valuesDat$pred_pct[counter] <- scores$pred_pct 
    valuesDat$pc_pct[counter] <- scores$pc_pct 
    valuesDat$pred_pos[counter]<- scores$pred_pos 
    valuesDat$sensitivity[counter]<- scores$sensitivity 
    valuesDat$specificity[counter]<- scores$specificity 
  }
  return(valuesDat)
}

#take input from the contigs and the coordinates and 
#rearrange so that start < stop
reformatContigPositionData <- function(dat){
  colnames(dat) <- c("contig", "start", "stop", "srna")
  dat <- dat  %>% select(contig, srna, start, stop) %>% mutate(strand = "+")
  dat <- dat %>% mutate(start = as.numeric(start), stop = as.numeric(stop))
  dat <- dat %>% mutate(tmpstart = ifelse(start < stop, start, stop),
                        tmpend = ifelse(start > stop, start, stop))
  return(dat)
}

#check for overlapping coordinates on the same contig for two given data sets
getOverlapIDs <- function(queryDat, targetDat){
  queryDat <- queryDat %>% arrange(start)
  targetDat <- targetDat %>% arrange(start)
  query <- GRanges(seqnames = queryDat$contig,
                   ranges = IRanges(start = queryDat$tmpstart, end = queryDat$tmpend),
                   strand = queryDat$strand, query_name = queryDat$srna)
  
  target <- GRanges(seqnames = targetDat$contig,
                    ranges = IRanges(start = targetDat$tmpstart, end = targetDat$tmpend),
                    strand = targetDat$strand, query_name = targetDat$srna)
  
  
  lookup1NC <- data.frame(id1 = queryDat$srna, queryHits = c(1:length(queryDat$srna)))
  lookup2PC <- data.frame(id2 = targetDat$srna, subjectHits = c(1:length(targetDat$srna)))
  
  tmp <- GenomicRanges::findOverlaps(query, target, type = 'any')
  
  tmp <- as.data.frame(tmp)
  
  tmp <- tmp %>% left_join(lookup1NC) %>% left_join(lookup2PC)
  tmp <- tmp %>% mutate_all(as.character)
  tmp <- tmp %>% filter(id1 != id2)
  
  smallDat <- tmp %>% select(id1, id2) %>% unique()
  return(smallDat)
}





#