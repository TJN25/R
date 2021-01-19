#!/usr/bin/env Rscript
library(dplyr, quietly = T, warn.conflicts = F)
dat <- read.table("~/phd/RNASeq/tmp/test.plot")

valuesDat <- dat %>%
  summarise_all(list(max)) %>% as.matrix() %>% t()

values <- valuesDat[,1]
outDat <- data.frame(mean = mean(values), max  = max(values))

write.table(x = outDat, file = "~/phd/RNASeq/tmp/test.values", row.names = F, col.names = F, quote = F)



