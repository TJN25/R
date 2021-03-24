#!/usr/bin/env Rscript
options(warn = -1)
library(dplyr, quietly = T, warn.conflicts = F)

args = commandArgs(trailingOnly=TRUE)
# print(args)

id <- args[1]
filepath <- args[2]

dat <- read.table(paste(filepath, "/tmp/test.plot", sep = ""))

valuesDat <- dat %>% as.matrix() %>% t() %>% as.data.frame() %>%
  summarise_all(list(max)) %>% as.matrix() %>% t()

colnames(valuesDat)[1] <- id


# values <- valuesDat[,1]
# outDat <- data.frame(mean = mean(values), max  = max(values), count = length(values[values > 0]))

write.table(x = valuesDat, file = paste(filepath, "/tmp/tmp.values", sep = ""), row.names = F, col.names = T, quote = F)


