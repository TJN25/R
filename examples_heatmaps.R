#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library(tidyverse, quietly = T, warn.conflicts = F))
suppressMessages(library(ggplot2, quietly = T, warn.conflicts = F))
suppressMessages(library(zoo, quietly = T, warn.conflicts = F))
suppressMessages(library(reshape2, quietly = T, warn.conflicts = F))

args = commandArgs(trailingOnly=TRUE)
# print(args)
id <- "GCA_000213655.1_267"
outname <- "GCA_000222975.1_486"
filepath <- "~/phd/RNASeq/examples"


id <- args[1]
outname <- args[2]
filepath <- args[3]

# for(outname in c("GCA_000015425.1_257", "GCA_000213655.1_267", "GCA_000196795.1_104", "GCA_000015425.1_172", "GCA_000438825.1_154",
#                  "GCA_000438825.1_15", "GCA_000222975.1_486", "GCA_002850215.1_141")){

dat <- read.table(paste(filepath, "/", outname, "_example_files/", outname,"_read_values.txt", sep = ""), as.is = T, stringsAsFactors = F, fill =T)

dat[is.na(dat)] <- 0
cols <- data.frame(names = colnames(dat), genomes = as.character(dat[1,]))
dat <- dat[-1,]
dat <- dat %>% mutate_all(as.numeric)
col_sum <- colSums(dat)
scaled.dat <- dat
scaled.dat <- as.data.frame(scaled.dat)
i <- 2
for(i in 1:ncol(scaled.dat)){
  zooDat <- zoo::zoo(as.numeric(as.character(scaled.dat[,i])))
  smoothDat <- zoo::rollapply(zooDat, width = 25, by = 1, FUN = mean, align = "center", partial = T) 
  scaled.dat[,i] <- smoothDat
}

scaled.dat <- scaled.dat %>% mutate_all(as.numeric) %>% mutate(row_num = row_number())
scaled.dat2 <- scaled.dat
scaled.dat2[is.na(scaled.dat2)] <- 0

# i <- 2
# for(i in 1:(ncol(scaled.dat2) -1)){
#   scaled.dat2[,i] <- scaled.dat2[,i]/sum(scaled.dat2[,i])
#   maxVal <- max(scaled.dat2[,i])
#   adjVal <- 1/maxVal
#   # scaled.dat2[,i] <- scaled.dat2[,i]*adjVal
# }



mat <- scaled.dat2 %>% as.matrix() 
meltDat <- melt(data = mat)
meltDat <- meltDat %>% mutate(value = ifelse(value > 5000, 5000, ifelse(value < 1, 1, value)))



# meltDat <- meltDat %>% filter(Var1 > 2648, Var1 < (10000 - 2363))
tmpDat <- data.frame(Var1 = (max(meltDat$Var1) + 1), Var2  ="V1", value = 5000)
meltDat <- meltDat %>% bind_rows(tmpDat)

p <- ggplot() +
  geom_tile(data = meltDat[meltDat$Var2 != "row_num",], aes(x = Var1, y = Var2, fill = log10(value))) +
  scale_fill_gradient(low = 'white', high = 'red') +
  ggtitle(outname)
# print(p)
ggsave(filename = paste(filepath, "/heatmaps/", outname,"_read_depths_heatmap_log.svg", sep = ""), plot = p, width = 450, height = 450, units = "mm")
# print(outname)
print(cols)
# }
# 
# df <- data.frame(a = c(1,5,10,50,100,500,1000,5000), b = c(1,2,3,4,5,6,7,8), c = c(1,1,1,1,1,1,1,1))
# ggplot() +
#   geom_tile(data =df, aes(x = b, y = c, fill = log10(a))) +
#   scale_fill_gradient(low = 'white', high = 'red') 
# 
# log10(df$a)
