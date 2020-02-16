library(ggplot2)
library(comparativeSRA)
library(tidyverse)




allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/escherichia/GCA_000017745_data/", gff_names = "gff_files", sra = "GCA_000017745.1")
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/escherichia/GCA_000017765.1_data/", gff_names = "gff_files", sra = "GCA_000017765.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/escherichia/GCA_000017985.1.data/", gff_names = "gff_files", sra = "GCA_000017985.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/escherichia/GCA_001559675.1.data/", gff_names = "gff_files", sra = "GCA_001559675.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/escherichia/GCA_002843685.1.data/", gff_names = "gff_files", sra = "GCA_002843685.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/escherichia/GCA_900186905.1.data/", gff_names = "gff_files", sra = "GCA_900186905.1", allData = allData)


allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Shigella/GCA_000007405.1.data/", gff_names = "gff_files", sra = "GCA_000007405.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Shigella/GCA_000283715.1.data/", gff_names = "gff_files", sra = "GCA_000283715.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Shigella/GCA_000497505.1.data/", gff_names = "gff_files", sra = "GCA_000497505.1", allData = allData)


allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Salmonella/GCA_000006945.2.data/", gff_names = "gff_files", sra = "GCA_000006945.2", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Salmonella/GCA_000210855.2.data/", gff_names = "gff_files", sra = "GCA_000210855.2", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Salmonella/GCA_002504125.1.data/", gff_names = "gff_files", sra = "GCA_002504125.1", allData = allData)

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Salmonella/GCA_002813995.1.data/", gff_names = "gff_files", sra = "GCA_002813995.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Salmonella/GCA_900184385.1.data/", gff_names = "gff_files", sra = "GCA_900184385.1", allData = allData)

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Enterobacter/GCA_001874505.1.data/", gff_names = "gff_files", sra = "GCA_001874505.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Enterobacter/GCA_002303275.1.data/", gff_names = "gff_files", sra = "GCA_002303275.1", allData = allData)

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Klebsiella/GCA_000220485.1.data/", gff_names = "gff_files", sra = "GCA_000220485.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/Klebsiella/GCA_002848605.1.data/", gff_names = "gff_files", sra = "GCA_002848605.1", allData = allData)

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/serratia/GCA_002220715.1.data/", gff_names = "gff_files", sra = "GCA_002220715.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/serratia/GCA_000438825.1_data/", gff_names = "gff_files", sra = "GCA_000438825.1", allData = allData)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/serratia/GCA_000747565.1_data/", gff_names = "gff_files", sra = "GCA_000747565.1", allData = allData) 

allReadDepths <- allData
save(allReadDepths, file = "allReadDepths.Rda")
ggplot() +
  geom_density(data = allData , aes(x = score_2, group = group, colour = group)) +
  xlim(c(0,10)) 

ggplot() +
  geom_density(data = allData %>% filter(group != "3") , aes(x = score_2, group = group, colour = group)) +
  xlim(c(0,10)) 
