library(tidyverse)
library(VennDiagram)
library(shiny)
library(ggplot2)
library(viridis)
library(RColorBrewer)
#library(stringr)
#library(plyr)
library(devtools)
#library(tidyr)
library(shinyjs)
library(shinyWidgets)
library(DT)
library(lubridate)
library(dplyr)
library(svglite)
library(genoPlotR)
library(drake)
library(ape)
library(Biostrings)
library(ggtree)
library(treeio)
library(geiger)
library(ROSE)
library(reshape2)
library(igraph)
library("viridis")



####

load(file = "~/bin/r_git/R/allMerged_version_8b.Rda")
load(file = "~/bin/r_git/R/allRandom_version_8b.Rda")

summaryCountsRNA <- function(dat){
  

  # dat <- dat %>% separate(col = target.name, into = c("target.file", "target.seq"), remove = F, extra = "merge", sep =":")
  # 
  # dat$target.file <- str_replace_all(dat$target.file, ".tbl", ".fasta")
  
 # dat <- dat %>% filter(E.value < 0.001)
  
  dat %>% filter(target.name != query.name) %>% nrow()
  

  
  
  
  
  threshold <- 1e-5
  current_data <- dat %>% arrange(E.value) %>% group_by(target.name, query.name) %>% top_n(n = 1, wt = -E.value)

  summaryCount <- current_data %>% group_by(query.name) %>% select(query.name, target.name) %>% unique %>% summarise(count = n())

  summaryCount <- summaryCount %>% group_by(count) %>% summarise(number.of.hits = n()) %>% mutate(unique.hits = number.of.hits/count)
  return(summaryCount)
 
}
contigLookup <- read.table("~/phd/RNASeq/sequences/contig_ids.lookup")
positiveLookup <- read.table("~/phd/RNASeq/srna_seqs/version_1/positive_control/single_seqs/srna_positive_ids.lookup")

head(contigLookup)
head(positiveLookup)


colnames(contigLookup) <- c("target.name", "target.genome")
colnames(positiveLookup) <- c("query.name", "query.description", "query.genome")



datPositive <- read.table("~/phd/RNASeq/srna_seqs/version_1/positive_control.tbl", comment.char = "#", fill = T, sep = "", header = F, quote = "", as.is = T)
datPredicted <- read.table("~/phd/RNASeq/srna_seqs/version_1/predicted.tbl", comment.char = "#", fill = T, sep = "", header = F, quote = "", as.is = T)
datNegative <- read.table("~/phd/RNASeq/srna_seqs/version_1/negative_control.tbl", comment.char = "#", fill = T, sep = "", header = F, quote = "", as.is = T)

# datPredictedKnown <- read.table("~/phd/RNASeq/srna_seqs/version_1/predicted_known.tbl", comment.char = "#", fill = T, sep = "", header = F, quote = "", as.is = T)
# datNegativeKnown <- read.table("~/phd/RNASeq/srna_seqs/version_1/negative_control_positive_control.tbl", comment.char = "#", fill = T, sep = "", header = F, quote = "", as.is = T)

datPositive <- datPositive %>% filter(V2 == "-")
datPredicted <- datPredicted %>% filter(V2 == "-")
datNegative <- datNegative %>% filter(V2 == "-")

datPositive <- datPositive[,c(1:16)]
datPredicted <- datPredicted[,c(1:16)]
datNegative <- datNegative[,c(1:16)]

colnames(datPositive) <- c("target.name", "accession", "query.name", "accession.2", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sq.len", "strand", "E.value", "score", "bias", "description.of.target")
colnames(datPredicted) <- c("target.name", "accession", "query.name", "accession.2", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sq.len", "strand", "E.value", "score", "bias", "description.of.target")
# colnames(datPredictedKnown) <- c("target.name", "accession", "query.name", "accession.2", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sq.len", "strand", "E.value", "score", "bias", "description.of.target")
colnames(datNegative) <- c("target.name", "accession", "query.name", "accession.2", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sq.len", "strand", "E.value", "score", "bias", "description.of.target")
# colnames(datNegativeKnown) <- c("target.name", "accession", "query.name", "accession.2", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sq.len", "strand", "E.value", "score", "bias", "description.of.target")

summaryCountPositive <- summaryCountsRNA(dat = datPositive)
summaryCountPredicted <- summaryCountsRNA(dat = datPredicted)
# summaryCountPredictedKnown <- summaryCountsRNA(dat = datPredictedKnown)
summaryCountNegative <-  summaryCountsRNA(dat = datNegative)
# summaryCountNegativeKnown <- summaryCountsRNA(dat = datNegativeKnown)





# tmp <- datPredictedKnown %>% select(query.name) %>% unique() %>% mutate(known = T)
# datPredictedAll <- datPredicted %>% full_join(tmp, by = "query.name")
# datPredictedAll <- datPredictedAll %>% mutate(known = ifelse(is.na(known), F, T))
# datPredictedNovel <- datPredictedAll %>% filter(known == F) %>% select(-known)
# datPredictedKnown <- datPredictedAll %>% filter(known == T) %>% select(-known)
# summaryCountPredictedNovel <- summaryCountsRNA(dat = datPredictedNovel)
# summaryCountPredictedKnown <- summaryCountsRNA(dat = datPredictedKnown)


datPositive <- datPositive %>% left_join(contigLookup, by = "target.name") %>% left_join(positiveLookup, by = "query.name")

datPositivencRNA <- datPositive %>% filter(grepl(pattern = "tRNA", query.description) == F, 
                                           grepl(pattern = "repeat_region", query.description) == F, 
                                           grepl(pattern = "tRNA", query.description) == F, 
                                           grepl(pattern = "repeat_region", query.description) == F, 
                                           grepl(pattern = "rRNA", query.description) == F, 
                                           grepl(pattern = "rRNA", query.description) == F)
summaryCountPositivencRNA <- summaryCountsRNA(dat = datPositivencRNA)



summaryPredictedMerged <- allMerged %>% filter(total > 25, new_feature ==T)%>% group_by(count) %>% summarise(unique.hits = n())
summaryPositiveControlMerged <- allMerged %>% filter(total > 25, new_feature ==F)%>% group_by(count) %>% summarise(unique.hits = n())
summaryRandomMerged <- allRandom %>% filter(total > 25)%>% group_by(count) %>% summarise(unique.hits = n())



summaryCountPositive <- summaryCountPositive %>% mutate(Group = "Positive Control") %>% mutate(count = as.numeric(count))
summaryCountPositivencRNA <- summaryCountPositivencRNA %>% mutate(Group = "Positive Control (tRNA removed)") %>% mutate(count = as.numeric(count))
summaryCountNegative <-summaryCountNegative %>% mutate(Group = "Negative Control")%>% mutate(count = as.numeric(count))
summaryCountPredicted <- summaryCountPredicted %>% mutate(Group = "Predicted sRNA")%>% mutate(count = as.numeric(count))
summaryCountPredictedKnown <- summaryCountPredictedKnown %>% mutate(Group = "Predicted sRNA (matching known sRNA)")%>% mutate(count = as.numeric(count))
summaryCountNegativeKnown <-summaryCountNegativeKnown %>% mutate(Group = "Negative Control (matching known sRNA)")%>% mutate(count = as.numeric(count))
summaryCountPredictedNovel <- summaryCountPredictedNovel %>% mutate(Group = "Predicted sRNA (not matching known sRNA)")%>% mutate(count = as.numeric(count))
summaryPredictedMerged <- summaryPredictedMerged %>% mutate(Group = "Positive Control (whole genome alignment)")%>% mutate(count = as.numeric(count))
summaryPositiveControlMerged <- summaryPositiveControlMerged %>% mutate(Group = "Predicted sRNA (whole genome alignment)")%>% mutate(count = as.numeric(count))
summaryRandomMerged <- summaryRandomMerged %>% mutate(Group = "Negative Control (whole genome alignment)")%>% mutate(count = as.numeric(count))


summaryCount <- summaryCountPositive %>% 
  bind_rows(summaryCountPositivencRNA)%>% 
  bind_rows(summaryCountNegative)%>% 
  bind_rows(summaryCountPredicted)#%>% 
  # bind_rows(summaryCountPredictedKnown)%>% 
  # bind_rows(summaryCountNegativeKnown)%>% 
  # bind_rows(summaryCountPredictedNovel)%>% 
  # bind_rows(summaryPredictedMerged)%>% 
  # bind_rows(summaryPositiveControlMerged)%>% 
  # bind_rows(summaryRandomMerged)

ggplot() +
  geom_point(data = summaryCount %>% filter(Group == "Positive Control" | Group == "Negative Control" | Group == "Predicted sRNA" | Group == "Negative Control" | Group == "Positive Control (whole genome alignment)") , aes(x = as.numeric(count),y = log(unique.hits), group = Group, fill = Group, colour = Group)) +
  xlab("Number of copies of sRNAs") +
  ylab("Number of different sRNAs with a given copy number")

# ggplot() +
#   geom_point(data = summaryCount %>% filter(Group == "Positive Control" | Group == "Negative Control (matching known sRNA)" | Group == "Predicted sRNA (not matching known sRNA)") , aes(x = as.numeric(count),y = log(unique.hits), group = Group, fill = Group, colour = Group)) +
#   xlab("Number of copies of sRNAs") +
#   ylab("Number of different sRNAs with a given copy number")


ggplot() +
  geom_point(data = summaryCount %>% filter(Group == "Positive Control (tRNA removed)" | Group == "Predicted sRNA (matching known sRNA)" | Group == "Predicted sRNA (not matching known sRNA)" ) , aes(x = as.numeric(count),y = log(unique.hits), group = Group, fill = Group, colour = Group)) +
  xlab("Number of copies of sRNAs") +
  ylab("Number of different sRNAs with a given copy number (log)")

ggplot() +
  geom_point(data = summaryCount %>% filter( Group == "Predicted sRNA (matching known sRNA)" | Group == "Predicted sRNA (not matching known sRNA)" | Group == "Negative Control") , aes(x = as.numeric(count),y = log(unique.hits), group = Group, fill = Group, colour = Group)) +
  xlab("Number of copies of sRNAs") +
  ylab("Number of different sRNAs with a given copy number (log)")

ggplot() +
  geom_point(data = summaryCount %>% filter( Group == "Positive Control (whole genome alignment)" | Group == "Predicted sRNA (whole genome alignment)" | Group == "Negative Control (whole genome alignment)") , aes(x = as.numeric(count),y = log(unique.hits), group = Group, fill = Group, colour = Group)) +
  xlab("Number of copies of sRNAs") +
  ylab("Number of different sRNAs with a given copy number (log)")

ggplot() +
  geom_point(data = summaryCount %>% filter( Group == "Positive Control (whole genome alignment)" | Group == "Positive Control") , aes(x = as.numeric(count),y = log(unique.hits), group = Group, fill = Group, colour = Group)) +
  xlab("Number of copies of sRNAs") +
  ylab("Number of different sRNAs with a given copy number (log)")




allMerged$count <- as.numeric(allMerged$count)

