#!/usr/bin/env Rscript
library(dplyr, quietly = T, warn.conflicts = F)

args = commandArgs(trailingOnly=TRUE)
print(args)
dat <- read.table(args[1])

valuesDat <- dat %>%
  summarise_all(list(max)) %>% as.matrix() %>% t()


values <- valuesDat[,1]
outDat <- data.frame(mean = mean(values), max  = max(values), count = length(values[values > 0]))

write.table(x = outDat, file = "~/phd/RNASeq/tmp/test.values", row.names = F, col.names = F, quote = F)



