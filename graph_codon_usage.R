print("Starting.")

library(tidyverse)
library(ggplot2)
setwd("/Volumes/2TTJN/Virus_Host_Prediction/")

options(warn = -1)


args <- commandArgs(trailingOnly = T)


# args <- c("/Volumes/2TTJN/Virus_Host_Prediction/Synechococcus_codon_usage.txt",
#           "/Volumes/2TTJN/Virus_Host_Prediction/phage_codon_usage.txt",
#           "Synechococcus",
#           "/Volumes/2TTJN/Virus_Host_Prediction/Synechococcus_codons_test_18-06-2018.png",
#           "test_Synechococcus.codons.txt")

print(args)

setwd(args[6])

codonLookup <- read.table("/Volumes/2TTJN/Virus_Host_Prediction/codon_lookup.txt", header = F, sep = "\t", comment.char = "", quote = "", as.is = T, fill = T)

colnames(codonLookup) <- c("codon", "amino.acid")

bacCodons <- read.table(args[1], header = F, sep = "\t", comment.char = "", quote = "", as.is = T, fill = T)

phageCodons <- read.table(args[2], header = F, sep = "\t", comment.char = "", quote = "", as.is = T, fill = T)

aminoAcidFreq <- read.table(args[5], header = T, sep = "\t", comment.char = "", quote = "", as.is = T, fill = T)

aminoAcidFreq <- aminoAcidFreq%>%
  filter(Genome != "Genome")%>%
  filter(Genome != "A")%>%select(-X.1)

aminoSummary <- data.frame(amino.acid = colnames(aminoAcidFreq)[2:ncol(aminoAcidFreq)],
                           freq = rep(NA, ncol(aminoAcidFreq) - 1))


for(i in 2:ncol(aminoAcidFreq)){
  aminoSummary[(i - 1), 2] <- mean(as.numeric(aminoAcidFreq[,i]))
}

aminoSummary <- aminoSummary%>%arrange(-freq)

aminoSummary <- aminoSummary%>%left_join(codonLookup, by = "amino.acid")


##check that each codon is in the same column
for(i in 1:ncol(bacCodons)){
  if(i %% 2 == 0){
    xx <- length(unique(bacCodons[, i]))
    if(xx > 1){
      print(paste("Column V", i, " has ", xx, " codons.", sep = ""))
      
    }
  }
}
for(i in 1:ncol(phageCodons)){
  if(i %% 2 == 0){
    xx <- length(unique(phageCodons[, i + 1]))
    if(xx > 1){
      print(paste("Column V", i, " has ", xx, " codons.", sep = ""))
    }
    }
}


##get genera infor from the phageCodons

phageCodons <- phageCodons%>%separate(col = V2, into = c("Genera", "Species", "Extra"), sep = " ", extra = "merge", remove = F)

SynechococcusCodons <- phageCodons%>%filter(Genera == args[3])

phageLength <- nrow(SynechococcusCodons)



##format the codons for plotting
phageDat <- data.frame(genera = rep(NA, 1), codon = rep(NA, 1), prop = rep(NA, 1))
for(i in 6:(ncol(SynechococcusCodons) - 1)){
  if(i %% 2 == 0){
    genera = SynechococcusCodons[,2]
    codon = SynechococcusCodons[,i]
    j <- i + 1
    prop = SynechococcusCodons[,j]

    df <- data.frame(genera = genera, codon = codon, prop = prop)

    phageDat <- phageDat%>%bind_rows(df)


  }
}

phageDat <- phageDat%>%
  mutate(type = rep("phage", nrow(phageDat)))%>%
  filter(!is.na(codon))


##format the codons for plotting
bacDat <- data.frame(genera = rep(NA, 1), codon = rep(NA, 1), prop = rep(NA, 1))
for(i in 1:(ncol(bacCodons) - 1)){
  if(i %% 2 == 0){
    genera = bacCodons[,1]
    codon = bacCodons[,i]
    j <- i + 1
    prop = bacCodons[,j]

    df <- data.frame(genera = genera, codon = codon, prop = prop)

    bacDat <- bacDat%>%bind_rows(df)


  }
}

bacDat <- bacDat%>%
  mutate(type = rep("bacteria", nrow(bacDat)))%>%
  filter(!is.na(codon))


bacSummary <- bacDat%>%group_by(codon)%>%summarise(bac_prop = mean(prop))

dat <- phageDat%>%left_join(bacSummary, by = "codon")

dat <- dat%>%mutate(log_odds = log2(prop/bac_prop))


logOddsSummary <- dat%>%select(codon, prop, bac_prop, log_odds)%>%unique()%>%
  group_by(codon)%>%summarise(score = sum(abs(log_odds)))

dat <- dat%>%left_join(aminoSummary, by = "codon")

dat <- dat%>%arrange(-freq)


tmp <- dat
tmp <- tmp%>%left_join(logOddsSummary, by = "codon")


tmp$codon <- factor(tmp$codon, levels = unique(tmp$codon))

levels(tmp$codon)



#log_odds_sum <- sum(abs(logOddsSummary$score))

tmp2 <- tmp%>%select(score, codon)%>%unique()


p <- ggplot(data = tmp,aes(x = codon, y = log_odds)) +
  geom_boxplot( outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
   geom_text(data = tmp2, aes(label = round(score), x = codon, y = (2.5 + score/50)))

p 

print("Saving png.")
ggsave(filename = args[4], plot = p, device = "png", width = 15, height = 10)

print("Done.")
