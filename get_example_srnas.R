#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library(tidyverse, quietly = T, warn.conflicts = F))
suppressMessages(library(stringi, quietly = T, warn.conflicts = F))
suppressMessages(library(ape, quietly = T, warn.conflicts = F))
suppressMessages(library(genoPlotR, quietly = T, warn.conflicts = F))
# library(dplyr, quietly = T, warn.conflicts = F)

args = commandArgs(trailingOnly=TRUE)
# print(args)
id <- "GCA_000213655.1_267"
outname <- "test"
filepath <- "~/phd/RNASeq/examples"

id <- args[1]
outname <- args[2]
filepath <- args[3]




plot_region <- function(i, adj.range = 5000, counter, dat = dat, write_data = F){
  
  genome <- as.character(dat$genome[i])
  # print(genome)
  genus <- dat$genus[i]
  gff  <- tryCatch({
    gff <- read.gff(paste("~/phd/RNASeq/annotation_files/", genome,".gff", sep = ""))
    
    gff
  }, error =  function(e) {
    print(paste("Error: ~/phd/RNASeq/annotation_files/", genome,".gff not found", sep = ""))
    return(F)
  })
  if(gff == F){
    return(counter)
  }
  
  start_val <- as.numeric(dat$start[i])
  end_val <- as.numeric(dat$end[i])
  contig <- dat$query.name[i]
  strand_val <- ifelse(start_val < end_val, '+', '-')
  
  subgff <- gff %>% filter(start > (start_val - adj.range), start < (start_val + adj.range), type == 'CDS', seqid == contig) %>% 
    separate(attributes, into = c(NA, "t1"), extra = 'merge', remove = F, sep = "product=") %>% separate(t1, into = "Name", sep = ";")
  
  
  names1 <- c("start", "end", subgff$Name, "rna_1")
  starts1 <- c((start_val - 5000), (start_val + 5000 - 1), subgff$start, min(start_val, end_val))
  ends1 <- c((start_val - 5000 + 1), (start_val + 5000),subgff$end, max(start_val, end_val))
  
  starts1[starts1 < (start_val - adj.range)] <- start_val - adj.range
  ends1[ends1 > (start_val + adj.range)] <- start_val + adj.range
  
  strands1 <- c("+", "-", as.character(subgff$strand), strand_val)
  
  # names1 <- subgff$t2
  # starts1 <- subgff$start
  # ends1 <- subgff$end
  # strands1 <- as.character(subgff$strand)
  
  df1 <- data.frame(name=names1, start=starts1, end=ends1,
                    strand=strands1, col="Blue")
  df1 <- df1 %>% mutate(col = ifelse(name == "rna_1", "Red", "Blue"),
                        fill = ifelse(name == "rna_1", "Red", "Blue"))
  
  
  dna_seg1 <- dna_seg(df1)
  is.dna_seg(dna_seg1)
  
  # if(i == 1){
  dna_segs <- list(dna_seg1)
  # }else{
  #   dna_segs[[i]] <- dna_seg1
  # }
  
  mid_pos <- middle(dna_segs[[1]])
  
  # comparison1 <- as.comparison(data.frame(start1=dna_segs[[1]]$start, end1=dna_segs[[1]]$end,
  #                                         start2=dna_segs[[2]]$start,
  #                                         end2=dna_segs[[2]]$end))
  
  annot2 <- genoPlotR::annotation(x1=c(mid_pos),
                                  x2=c(dna_segs[[1]]$end),
                                  text=c(dna_segs[[1]]$name),
                                  rot=c(90), col=c("black"))
  
  
  if(write_data){
    svg(filename=paste(filepath, "/", outname, "_example_files/", genome, "_", start_val, ".svg", sep = ""),
        width=20,
        height=20,
        pointsize=12)
    plot_gene_map(dna_segs=dna_segs,
                  annotations=annot2, annotation_height=1.3)
    
    dev.off()
  }else{
    plot_gene_map(dna_segs=dna_segs,
                  annotations=annot2, annotation_height=1.3)
  }
  
  counter <- counter + 1
  return(counter)
  
}

# print(id)
# print(outname)
# print(filepath)

contig_labels <- read.table("~/phd/RNASeq/genome_contig_pairs.txt")
colnames(contig_labels) <- c("query.name", "genome")

dat <- read.table(paste("~/phd/RNASeq/srna_seqs/version_1/predicted/large_alignments/alignments_rnaalifold/alignments_", id, ".stk", sep = ""), comment.char = "#", fill = T)

dat <- dat %>% separate(V1, into = c("query.name", "location"), sep = "/", remove = T, extra = 'drop') %>%  select(query.name, location) %>% 
  mutate(found = T) %>% mutate(rev.query.name = stri_reverse(query.name)) %>%  separate(rev.query.name, into = c("t1"), sep = "\\|", remove = F) %>% mutate(query.name = stri_reverse(t1))

dat <- dat %>% left_join(contig_labels, by = "query.name")

genus_labels <- read.table("~/phd/RNASeq/contig_genus_lables.txt")

colnames(genus_labels) <- c("query.name", "genus")

genus_labels$genus[genus_labels$genus == "Methylococcaceae"] <- "Moraxella"

dat <- dat %>% left_join(genus_labels, by = "query.name") %>%  unique() %>%  separate(location, into = c("start", "end"), sep = "-", remove = F) %>% filter(!is.na(genome)) %>% mutate(srna = id)

summaryDat <- dat %>% select(genome, genus) %>% unique() %>% filter(!is.na(genome))

subsetDat <- dat %>% filter(grepl(pattern = "GCA", x = genome) == T)
write.csv(x = subsetDat, file = paste(filepath, "/", outname, "_example_files/", outname, "_list.csv", sep = ""), row.names = F, quote = F)


counter <- 1
i <- 1
for(i in 1:nrow(subsetDat)){
  # print(i)
  counter <- plot_region(i = i, counter = counter, dat = subsetDat, write_data = T, adj.range = 5000)
  if(counter > 20){
    break
  }
}
