#!/usr/bin/env Rscript
library(dplyr, quietly = T, warn.conflicts = F)

args = commandArgs(trailingOnly=TRUE)
print(args)
try(dat <- read.table(args[1], fill = T))

try(valuesDat <- dat %>% summarise_all(list(max)) %>% as.matrix() %>% t())

try(values <- valuesDat[,1])
try(outDat <- data.frame(mean = mean(values), max  = max(values), count = length(values[values > 0])))

try(write.table(x = outDat, file = "~/phd/RNASeq/tmp/test.values", row.names = F, col.names = F, quote = F))



