library(tidyverse)
setwd("~/phd/pred_virus_host/")
hmm <- read.table("arVOG_archaeal_1")


lookup <- readLines("archaeal_virus_lookup.txt")
lookup <- data.frame(text = lookup)
lookup <- lookup%>%
  separate(col = text, into = c("id", "description"), sep = " ", remove = F, extra = "merge")%>%
  separate(col = description, into = c("protein", "genome"), sep = "\\[", remove = F, extra = "merge")





models <- readLines("archaeal_models_1_renamed.txt")
models <- data.frame(text = models, stringsAsFactors = F)

models <- models%>%
  separate(col = text, into = c("id", "i1"), sep = " ", remove = F, extra = "merge")%>%
  separate(col = i1, into = c("protein", "i2"), sep = "\\[", remove = T, extra = "merge")%>%
  separate(col = i2, into = c("genome", "model"), sep = "\\]", remove = T, extra = "merge")

hypotheicals <- models%>%filter(grepl("hypothetical", protein) == T)
ORFs <- models%>%filter(grepl("ORF", protein) == T)
VPs <- models%>%filter(grepl("VP", protein) == T)
other <- models%>%filter(grepl("hypothetical", protein) == F)%>%filter(grepl("ORF", protein) == F)%>%filter(grepl("VP", protein) == F)


hypotheicals <- hypotheicals%>%bind_rows(ORFs)%>%bind_rows(VPs)

hypotheicals$protein <- rep(" ", nrow(hypotheicals))

models <- hypotheicals%>%bind_rows(other)


proteins <- models%>%group_by(model)%>%summarise(protein.list = paste(unique(protein), collapse = ","))
proteinCounts <- models%>%group_by(model)%>%select(model, protein)%>%unique()%>%summarise(protein.count = n())

genomes <- models%>%group_by(model)%>%summarise(genome.list = paste(unique(genome), collapse = ","))
genomeCounts <- models%>%group_by(model)%>%select(model, genome)%>%unique()%>%summarise(genome.count = n())


modelSummary <- proteins%>%left_join(proteinCounts)%>%left_join(genomes)%>%left_join(genomeCounts)%>%
  arrange(model)

write.table(modelSummary, "archaeal_1_model_summary.txt", quote = F, col.names = F, row.names = F, sep = "\t")

#
hmm <- readLines("arVOG_archaeal_1_full.txt")
hmm <- data.frame(text = hmm)

hmm <- hmm%>%separate(text, into = c("i1", "i2"), sep = "\\[", remove = F, extra = "merge")%>%
  separate(i1, into = c("target.name", "accession",  "query.name", "accession.2", "E.value" , "score" , "bias" ,  "E.value.2" , "score.2",  "bias.2" ,  "exp", "reg", "clu",  "ov", "env", "dom", "rep", "inc", "description.of.target", "protein.name"),
           sep = " ", remove = T, extra = "merge")%>%
  separate(col = i2, into = c("genome", "i3"), remove = T, sep = "\\]", extra = "merge")%>%
  separate(col = i3, into =c("i4", "protein.list", "protein.count", "genome.list", "genome.count"), sep = " ", remove = T, extra = "merge")

#
