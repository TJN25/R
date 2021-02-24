
groupOverlapItems <- function(smallDat, item, current_ids){
  # print(length(current_ids))
  # print(item)
  if(length(current_ids) == 0){
    current_ids <- item
  }
  df <- smallDat %>% filter(id1 == item)
  
  if(nrow(df) == 0){
    return(current_ids)
  }
  df$seen <- df$id2 %in% current_ids
  df <- df %>% filter(seen == F)
  
  if(nrow(df) == 0){
    return(current_ids)
  }
  for(item2 in df$id2){
    current_ids <- unique(c(current_ids, item2))
    current_ids <- groupOverlapItems(smallDat, item2, current_ids)
  }
  return(current_ids)
}


dat <- read.table("~/phd/RNASeq/srnas/known_predicted/predicted_genomic_sequence_matches.txt", sep = " ", fill = T, as.is = T)

dat <- dat %>% select(V1, V4)
dat <- dat %>% separate(col = V1, into = c("contig", "coordinates"), sep = "\\/", remove = F)
dat <- dat %>% separate(col = coordinates, into = c("start", "stop"), sep = "-", remove = F)
dat <- dat %>% dplyr::rename(srna = V4) %>% select(contig, srna, start, stop) %>% mutate(strand = "+")
colnames(dat) <- c("contig","srna", "start", "stop", "strand") 

dat <- dat %>% mutate(start = as.numeric(start), stop = as.numeric(stop))

dat <- dat %>% mutate(tmpstart = ifelse(start < stop, start, stop),
                      tmpend = ifelse(start > stop, start, stop))

query <- GRanges(seqnames = dat$contig,
                 ranges = IRanges(start = dat$tmpstart, end = dat$tmpend),
                 strand = dat$strand, query_name = dat$srna)


lookup1 <- data.frame(id1 = dat$srna, queryHits = c(1:length(dat$srna)))
lookup2 <- data.frame(id2 = dat$srna, subjectHits = c(1:length(dat$srna)))

save(query, file="~/phd/RNASeq/r_files/query.Rda")
save(lookup1, file="~/phd/RNASeq/r_files/lookup1.Rda")
save(lookup2, file="~/phd/RNASeq/r_files/lookup1.Rda")

range_out <- GenomicRanges::findOverlaps(query, query, type = 'any')
save(range_out, file="~/phd/RNASeq/r_files/range_out.Rda")

tmpDat <- as.data.frame(range_out)

rangesDat <- tmpDat %>% left_join(lookup1) %>% left_join(lookup2)

rangesDat <- rangesDat %>% filter(id1 != id2)

smallDat <- rangesDat %>% select(id1, id2) %>% unique()

save(smallDat, file="~/phd/RNASeq/r_files/smallDat.Rda")
save(rangesDat, file="~/phd/RNASeq/r_files/rangesDat.Rda")


smallDat <- smallDat %>% mutate_all(as.character)

s2 <- smallDat
s2$id1 <- smallDat$id2
s2$id2 <- smallDat$id1

smallDat <- smallDat %>% bind_rows(s2) %>% unique()

lookup1$id1[lookup1$id1 == "GCA_000015425.1_1002"]

ids <- as.character(unique(dat$srna))
item <- ids[68]
item2 <- "alignments_GCA_000017765.1_689"

ids[1:10]
idsDat <- data.frame(ids = ids)

item <- ids[ids == "GCA_000015425.1_1002"]
tmp <- groupOverlapItems(smallDat, item, current_ids = c())

checked_ids <- c()
for(item in ids){
  if(item %in% checked_ids){
    next
  }
  print(item)
  current_ids <- groupOverlapItems(smallDat, item, current_ids = c())
  
  checked_ids <- c(checked_ids, current_ids)
  
  write.table(x = current_ids, file = paste("~/phd/RNASeq/srnas/known_predicted/combined_alignments_ids/", current_ids[1], "_combined_list.txt", sep = ""), append = T, quote = F, row.names = F, col.names = F)
  
  # mat <- as.data.frame(t(as.matrix(current_ids)))
  # write.table(mat, file = "~/phd/RNASeq/tmp/nc_overlaps", append = T, quote = F, col.names = F, row.names = F, sep = "\t")
}