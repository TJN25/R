library(tidyverse)

hmmRes <- read.table("~/phd/pred_virus_host/arVOG_fungi_formatted.tab", sep = "", comment.char = "#", quote = "", fill = T, as.is = T)
trainingRes <- read.table("~/phd/pred_virus_host/arVOG_archaeal_train.txt", sep = "", comment.char = "#", quote = "", fill = T, as.is = T)

colnames(hmmRes) <- c("target.name", "accession", "query.name", "accession.2", "E.value", "score", "bias", "E.value.2", "score.2", "bias.2", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description.of.target")
colnames(trainingRes) <- c("target.name", "accession", "query.name", "accession.2", "E.value", "score", "bias", "E.value.2", "score.2", "bias.2", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description.of.target")

suppressWarnings(hmmRes <- hmmRes%>%filter(query.name != "")%>%mutate(score = as.numeric(score))%>%filter(!is.na(score)))
suppressWarnings(trainingRes <- trainingRes%>%filter(query.name != "")%>%mutate(score = as.numeric(score))%>%filter(!is.na(score)))


modelMean <- hmmRes%>%group_by(query.name)%>%summarise(mean = mean(as.numeric(score)))
modelSD <- hmmRes%>%group_by(query.name)%>%summarise(sd = sd(as.numeric(score)))
modelCount <- hmmRes%>%group_by(query.name)%>%summarise(count = n())

modelSummary <- modelMean%>%full_join(modelSD, by = "query.name")%>%full_join(modelCount, by = "query.name")

ggplot() +
  geom_histogram(data = modelSummary, aes(x = mean, y = ..density..), binwidth = 1) +
  geom_histogram(data = modelSummary, aes(x = upper, y = ..density..), binwidth = 1) +
  geom_histogram(data = modelSummary, aes(x = lower, y = ..density..), binwidth = 1)


ggplot() +
  geom_histogram(data = hmmRes, aes(x = score, y = ..density..), binwidth = 1)




trainingMean <- trainingRes%>%group_by(query.name)%>%summarise(mean.train = mean(as.numeric(score)))
trainingSD <- trainingRes%>%group_by(query.name)%>%summarise(sd.train = sd(as.numeric(score)))
trainingCount <- trainingRes%>%group_by(query.name)%>%summarise(count.train = n())

trainingSummary <- trainingMean%>%full_join(trainingSD, by = "query.name")%>%full_join(trainingCount, by = "query.name")


summary <- trainingSummary%>%full_join(modelSummary, by = "query.name")%>%
  mutate(difference = (mean.train - mean)/mean.train)%>%
  mutate(x = mean - mean.train)

ggplot() +
  geom_histogram(data = trainingRes, aes(x = score, y = ..density..), binwidth = 1) +
  geom_histogram(data = hmmRes, aes(x = score, y = ..density..), binwidth = 1, fill = "Red", alpha = 0.3) +
  xlim(0, 200)

ggplot() +
  geom_histogram(data = summary, aes(x = x, y = ..density..), binwidth = 10) 



median(hmmRes$score)
sd(hmmRes$score)

