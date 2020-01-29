#!/usr/bin/env Rscript
library('getopt')
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help' , 'h', 0, "logical",
  'count' , 'c', 1, "integer",
  'mean' , 'm', 1, "double",
  'sd' , 's', 1, "double"
), byrow=TRUE, ncol=4)


opt = getopt(spec)

print(opt)
