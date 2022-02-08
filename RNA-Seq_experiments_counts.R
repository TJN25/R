library(tidyverse)

##working out where the numbers in Figure 1 came from as a repeatable approach
##The flow diagram states:
#   strain of genome in experiment is reported. Was this manually checked?
#   >4 experiments are available for each strain
#   each genus contained at least 2 genomes after filtering (this was ignored 
#     for a couple of genomes to fill in a large taxonomic gap)

load('r_files/SRA_bacteria_RNASeq_v2_all.Rda')

filterData <- function(dat){
  dat <- dat %>% filter(grepl(pattern = 'RNA-Seq', x = Assay.Type, ignore.case = T),
                        grepl(pattern = 'lumina', x= Instrument, ignore.case = T), 
                        LibraryLayout == 'PAIRED')
  return(dat)
}

getGenusCounts <- function(dat){
  dat <- dat %>% separate(col = Organism, into = c("genus"), sep = " ", extra = 'drop', remove = F)
  dat <- dat %>% group_by(Organism, genus) %>% summarise(experiments.count = n()) %>% 
    ungroup() %>% group_by(genus) %>% 
    mutate(genomes.count = sum(experiments.count))
  genera <- dat %>% select(genus, genomes.count) %>% unique()
  returnList <- list(dat, genera)
  return(returnList)
}


sraFiltered <- filterData(SRA_bacteria_RNASeq_v2_all)
returnList <- getGenusCounts(sraFiltered)

experimentsCounts <- returnList[[1]]
genomeCounts <- returnList[[2]]
