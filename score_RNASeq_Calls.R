library(ggplot2)
library(comparativeSRA)
library(tidyverse)




allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Escherichia/GCA_000017745.1.data/", gff_names = "gff_files", sra = "GCA_000017745.1", mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Escherichia/GCA_000017765.1.data/", gff_names = "gff_files", sra = "GCA_000017765.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Escherichia/GCA_000017985.1.data/", gff_names = "gff_files", sra = "GCA_000017985.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Escherichia/GCA_001559675.1.data/", gff_names = "gff_files", sra = "GCA_001559675.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Escherichia/GCA_002843685.1.data/", gff_names = "gff_files", sra = "GCA_002843685.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Escherichia/GCA_900186905.1.data/", gff_names = "gff_files", sra = "GCA_900186905.1", allData = allData, mean_only = F)


allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Shigella/GCA_000007405.1.data/", gff_names = "gff_files", sra = "GCA_000007405.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Shigella/GCA_000283715.1.data/", gff_names = "gff_files", sra = "GCA_000283715.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Shigella/GCA_000497505.1.data/", gff_names = "gff_files", sra = "GCA_000497505.1", allData = allData, mean_only = F)


allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Salmonella/GCA_000006945.2.data/", gff_names = "gff_files", sra = "GCA_000006945.2", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Salmonella/GCA_000210855.2.data/", gff_names = "gff_files", sra = "GCA_000210855.2", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Salmonella/GCA_002504125.1.data/", gff_names = "gff_files", sra = "GCA_002504125.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Salmonella/GCA_002813995.1.data/", gff_names = "gff_files", sra = "GCA_002813995.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Salmonella/GCA_900184385.1.data/", gff_names = "gff_files", sra = "GCA_900184385.1", allData = allData, mean_only = F)

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Enterobacter/GCA_001874505.1.data/", gff_names = "gff_files", sra = "GCA_001874505.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Enterobacter/GCA_002303275.1.data/", gff_names = "gff_files", sra = "GCA_002303275.1", allData = allData, mean_only = F)

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Klebsiella/GCA_000220485.1.data/", gff_names = "gff_files", sra = "GCA_000220485.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Klebsiella/GCA_002848605.1.data/", gff_names = "gff_files", sra = "GCA_002848605.1", allData = allData, mean_only = F)

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Serratia/GCA_002220715.1.data/", gff_names = "gff_files", sra = "GCA_002220715.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Serratia/GCA_000438825.1_data/", gff_names = "gff_files", sra = "GCA_000438825.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Serratia/GCA_000747565.1_data/", gff_names = "gff_files", sra = "GCA_000747565.1", allData = allData, mean_only = F) 

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Acinetobacter/GCA_900088705.1.data/", gff_names = "gff_files", sra = "GCA_900088705.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Acinetobacter/GCA_000015425.1.data/", gff_names = "gff_files", sra = "GCA_000015425.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Acinetobacter/GCA_000046845.1.data/", gff_names = "gff_files", sra = "GCA_000046845.1", allData = allData, mean_only = F) 

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Pseudomonas/GCA_900243355.1.data/", gff_names = "gff_files", sra = "GCA_900243355.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Pseudomonas/GCA_000006765.1.data/", gff_names = "gff_files", sra = "GCA_000006765.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Pseudomonas/GCA_000281215.1.data/", gff_names = "gff_files", sra = "GCA_000281215.1", allData = allData, mean_only = F) 
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Pseudomonas/GCA_000007565.2.data/", gff_names = "gff_files", sra = "GCA_000007565.2", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Pseudomonas/GCA_001648195.1.data/", gff_names = "gff_files", sra = "GCA_001648195.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Pseudomonas/GCA_002208745.1.data/", gff_names = "gff_files", sra = "GCA_002208745.1", allData = allData, mean_only = F) 

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Xanthomonas/GCA_002850215.1.data/", gff_names = "gff_files", sra = "GCA_002850215.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Xanthomonas/GCA_001042875.1.data/", gff_names = "gff_files", sra = "GCA_001042875.1", allData = allData, mean_only = F)

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Lysobacter/GCA_002355295.1.data/", gff_names = "gff_files", sra = "GCA_002355295.1", allData = allData, mean_only = F) 

allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Alteromonas/GCA_002849875.1.data/", gff_names = "gff_files", sra = "GCA_002849875.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Alteromonas/GCA_000310085.1.data/", gff_names = "gff_files", sra = "GCA_000310085.1", allData = allData, mean_only = F)
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Alteromonas/GCA_000213655.1.data/", gff_names = "gff_files", sra = "GCA_000213655.1", allData = allData, mean_only = F) 
allData <- scoreRNAReadDepth(file_path = "~/phd/RNASeq/genera/Alteromonas/GCA_001886455.1.data/", gff_names = "gff_files", sra = "GCA_001886455.1", allData = allData, mean_only = F)



allReadDepthValues <- allData
save(allReadDepthValues, file = "~/bin/r_git/R/allReadDepthValues.Rda")
allReadDepths <- allData
save(allReadDepths, file = "~/bin/r_git/R/allReadDepths.Rda")
ggplot() +
  geom_density(data = allData , aes(x = score_2, group = group, colour = group)) +
  xlim(c(0,10)) 

ggplot() +
  geom_density(data = allData %>% filter(group != "3") , aes(x = score_2, group = group, colour = group)) +
  xlim(c(0,10)) 
