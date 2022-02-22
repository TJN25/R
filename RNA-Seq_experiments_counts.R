library(tidyverse)

##working out where the numbers in Figure 1 came from as a repeatable approach
##The flow diagram states:
#   strain of genome in experiment is reported. Was this manually checked?
#   >4 experiments are available for each strain
#   each genus contained at least 2 genomes after filtering (this was ignored 
#     for a couple of genomes to fill in a large taxonomic gap)



#take the input data and filter based on criteria
filterData <- function(dat){
  dat <- dat %>% filter(grepl(pattern = 'RNA-Seq', x = Assay.Type, ignore.case = T),
                        grepl(pattern = 'lumina', x= Instrument, ignore.case = T), 
                        LibraryLayout == 'PAIRED')
  return(dat)
}

#count the number of experiments and genomes found
getGenusCounts <- function(dat){
  dat <- dat %>% separate(col = Organism, into = c("genus"), sep = " ", extra = 'drop', remove = F)
  dat <- dat %>% group_by(Organism, genus) %>% summarise(experiments.count = n()) %>% 
    ungroup() %>% group_by(genus) %>% 
    mutate(genomes.count = sum(experiments.count))
  genera <- dat %>% select(genus, genomes.count) %>% unique()
  returnList <- list(dat, genera)
  return(returnList)
}

#the columns
getGenusColumns <- function(dat, genus_names){
  genus_vals <- list(rows = rep(NA, nrow(dat)), 
                     cols = rep(NA, nrow(dat)), 
                     organism = rep(NA, nrow(dat)),
                     run = rep(NA, nrow(dat)),
                     layout = rep(NA, nrow(dat)),
                     instrument = rep(NA, nrow(dat)))
  dat <- dat %>% mutate_all(as.character)
  for(row in 1:nrow(dat)){
    print(row)
    for(name in genus_names){
      for(col in 1:ncol(dat)){
        if(is.na(dat[row, col])){
          next
        }
        layout_found <-  'PAIRED' == dat[row, col] | 'SINGLE' == dat[row, col]
        if(layout_found){
          genus_vals$layout[row] <- dat[row, col]
          next
        }
        instrument_found <- grepl(pattern = "Illumina", x = dat[row, col], ignore.case = T)
        if(instrument_found){
          genus_vals$instrument[row] <- dat[row, col]
          next
        }
        name_found <- grepl(pattern = name, x = dat[row, col])
        if(name_found){
          genus_vals$rows[row] <- row
          genus_vals$cols[row] <- col
          genus_vals$run[row] <- dat[row, 1]
          genus_vals$organism[row] <- name
          next
        }
        
      }
    }
  }
  
  return(genus_vals)
}

load('r_files/SRA_bacteria_RNASeq_v2_all.Rda')

genus_names <- c("Escherichia", "Shigella", "Salmonella", "Enterobacter", "Klebsiella",
                 "Plautia", 'Citrobacter', 'Serratia', 'Wigglesworthia', 'Buchnera', 
                 'Sodalis', 'Lonsdelea', 'Brenneria', 'Dickeya', 'Edwardsiella', 
                 'Pantoea', 'Erwinia', 'Proteus', 'Providencia', 'Xenorhabdus', 
                 'Photorhabdus', 'Mannheimia', 'Aggregatibacter', 'Alishewanella', 
                 'Acinetobacillus', 'Yersinia', 'Vibrio', 'Alteromonas', 'Agarivoranas', 
                 'Pseudoalteromonas', 'Moritella', 'Shewanella', 'Psychrobacter', 
                 'Moraxella', 'Acinetobacter', 'Pseudomonas', 'Marinobacter', 'Methylomicrobium', 
                 'Methyloccocus', 'Francisella', 'Pseudoxanthomonas', 'Stenotrophomonas', 
                 'Xanthomonas', 'Xylella', 'Lysobacter')


if(FALSE){

sraFiltered <- filterData(SRA_bacteria_RNASeq_v2_all)
returnList <- getGenusCounts(sraFiltered)

experimentsCounts <- returnList[[1]]
genomeCounts <- returnList[[2]]

#this has been taken out of a function
dat <- SRA_bacteria_RNASeq_v2_all
genus_vals <- list(rows = rep(NA, nrow(dat)), 
                   cols = rep(NA, nrow(dat)), 
                   organism = rep(NA, nrow(dat)),
                   run = rep(NA, nrow(dat)),
                   layout = rep(NA, nrow(dat)),
                   instrument = rep(NA, nrow(dat)))
dat <- dat %>% mutate_all(as.character)
for(row in 1:nrow(dat)){
  genus_vals$run[row] <- dat[row, 1]
  for(name in genus_names){
    for(col in 1:ncol(dat)){
      if(is.na(dat[row, col])){
        next
      }
      layout_found <-  'PAIRED' == dat[row, col] | 'SINGLE' == dat[row, col]
      if(layout_found){
        genus_vals$layout[row] <- dat[row, col]
        next
      }
      instrument_found <- grepl(pattern = "Illumina", x = dat[row, col], ignore.case = T)
      if(instrument_found){
        genus_vals$instrument[row] <- dat[row, col]
        next
      }
      name_found <- grepl(pattern = name, x = dat[row, col])
      if(name_found){
        genus_vals$rows[row] <- row
        genus_vals$cols[row] <- col
        genus_vals$organism[row] <- dat[row, col]
        next
      }
      
    }
  }
}






check_rows <- as.numeric(row.names(genus_dat[is.na(genus_dat$rows),]))



for(i in check_rows){
  print(i)
  genus_vals$run[row] <- dat[row, 1]
  for(name in genus_names){
    for(col in 1:ncol(dat)){
      if(is.na(dat[row, col])){
        next
      }
      layout_found <-  'PAIRED' == dat[row, col] | 'SINGLE' == dat[row, col]
      if(layout_found){
        genus_vals$layout[row] <- dat[row, col]
        next
      }
      instrument_found <- grepl(pattern = "Illumina", x = dat[row, col], ignore.case = T)
      if(instrument_found){
        genus_vals$instrument[row] <- dat[row, col]
        next
      }
      name_found <- grepl(pattern = name, x = dat[row, col])
      if(name_found){
        genus_vals$rows[row] <- row
        genus_vals$cols[row] <- col
        genus_vals$organism[row] <- dat[row, col]
        next
      }
      
    }
  }
}

save(genus_vals, file = 'r_files/rnaseq_experiments_list_all.Rda')
}
## Working from here!

load('r_files/rnaseq_experiments_list_all.Rda')
load('r_files/SRA_bacteria_RNASeq_v2_all.Rda')
dat <- SRA_bacteria_RNASeq_v2_all %>% mutate_all(as.character)

# originally had the function call genus_vals <- getGenusColumns(SRA_bacteria_RNASeq_v2_all, genus_names)
# this failed so took out the loop so that I can restart mid loop
# should really learn a better way of writing this.


plot_files <- read.table('~/bin/PhD/Chapter_4/chapter_4_files/plot_file_names.txt')
colnames(plot_files) <- c('genus.name', 'genome.accession', 'run')

genus_dat <- as.data.frame(genus_vals)

genus_dat <-  genus_dat %>% full_join(plot_files, by = 'run')

genus_dat <- genus_dat[!is.na(genus_dat$rows) | !is.na(genus_dat$genome.accession),]

for(i in 1:nrow(genus_dat)){
  if(is.na(genus_dat$rows[i]) | is.na(genus_dat$cols[i])){
    next
  }
  genus_dat$organism[i] <- dat[genus_dat$rows[i], genus_dat$cols[i]]
}




filtered_dat <- genus_dat[genus_dat$layout == 'PAIRED' | is.na(genus_dat$layout),]
filtered_dat <- filtered_dat[!is.na(filtered_dat$rows) | !is.na(filtered_dat$genome.accession),]


# next step ---------------------------------------------------------------

#Moraxella not working

#there appear to be a number of runs missing
#some are labelled as experiments, 
#some will be from the original dataset
#need to track them all down


x  <- c(1:10)
y  <- x^2

ggplot() +
  geom_point(aes(x, y), color = "Red")

