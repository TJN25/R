##set up and data import #####
options(warn=-1)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tjnFunctions)
library(zoo)

mround <- function(x,base){
  base*round(x/base)
}


simulateSmallContigs <- function(arVOGres, baPOGres, euVOGres, lookup,
                                 arVOGModels = archaealModels,
                                 baPOGModels = phageModels,
                                 euVOGModels = eukaryoticModels,
                                 discriminant_models_only = F){

  ptm1 <- proc.time()
  aa <- arVOGres%>%full_join(lookup, by = "target.name")
  pa <- baPOGres%>%full_join(lookup, by = "target.name")
  ea <- euVOGres%>%full_join(lookup, by = "target.name")

  aa <- aa%>%full_join(arVOGModels, by = "query.name")
  pa <- pa%>%full_join(baPOGModels, by = "query.name")
  ea <- ea%>%full_join(euVOGModels, by = "query.name")


  ##adjust the model scores
  aa <- aa%>%mutate(score = score*score.percentage)%>%mutate(score = ifelse(is.na(score), 0, score))%>%mutate(query.name = ifelse(is.na(query.name), "arVOG_0", query.name))
  pa <- pa%>%mutate(score = score*score.percentage)%>%mutate(score = ifelse(is.na(score), 0, score))%>%mutate(query.name = ifelse(is.na(query.name), "baPOG_0", query.name))
  ea <- ea%>%mutate(score = score*score.percentage)%>%mutate(score = ifelse(is.na(score), 0, score))%>%mutate(query.name = ifelse(is.na(query.name), "euVOG_0", query.name))

  if(discriminant_models_only == T){
    aa <- aa%>%mutate(score = ifelse(score.percentage < 0.8, 0, score*score.percentage))%>%mutate(score = ifelse(is.na(score), 0, score))%>%mutate(query.name = ifelse(is.na(query.name), "arVOG_0", query.name))
    pa <- pa%>%mutate(score = ifelse(score.percentage < 0.8, 0, score*score.percentage))%>%mutate(score = ifelse(is.na(score), 0, score))%>%mutate(query.name = ifelse(is.na(query.name), "baPOG_0", query.name))
    ea <- ea%>%mutate(score = ifelse(score.percentage < 0.8, 0, score*score.percentage))%>%mutate(score = ifelse(is.na(score), 0, score))%>%mutate(query.name = ifelse(is.na(query.name), "euVOG_0", query.name))
    
  }else{
    aa <- aa%>%mutate(score = score*score.percentage)%>%mutate(score = ifelse(is.na(score), 0, score))%>%mutate(query.name = ifelse(is.na(query.name), "arVOG_0", query.name))
    pa <- pa%>%mutate(score = score*score.percentage)%>%mutate(score = ifelse(is.na(score), 0, score))%>%mutate(query.name = ifelse(is.na(query.name), "baPOG_0", query.name))
    ea <- ea%>%mutate(score = score*score.percentage)%>%mutate(score = ifelse(is.na(score), 0, score))%>%mutate(query.name = ifelse(is.na(query.name), "euVOG_0", query.name))
    
  }
  
  ##keep only the top matching model

  aa1 <- aa%>%arrange(-score)%>%group_by(target.name)%>%top_n(n = 1, wt = score)
  pa1 <- pa%>%arrange(-score)%>%group_by(target.name)%>%top_n(n = 1, wt = score)
  ea1 <- ea%>%arrange(-score)%>%group_by(target.name)%>%top_n(n = 1, wt = score)



  ## score the archaeal results ##

  ##create a data frame to store the output
  res <- data.frame(genome = NA, archaeal_score = NA, phage_score = NA, eukaryotic_score = NA, call = NA, protein_count = NA)


  ##simulate varying numbers of proteins
  for(j in 1:10){
    for(i in 1:25){
      printRemaining(i = ((j-1)*25 + i), length = 250, increment = 5)
      ##select i proteins from genomes
      aaTmp <- aa1%>%group_by(genome)%>%mutate(protein_count = n())%>%filter(protein.count.full >= i)%>%sample_n(i)
      paTmp <- pa1%>%group_by(genome)%>%mutate(protein_count = n())%>%filter(protein.count.full >= i)%>%sample_n(i)
      eaTmp <- ea1%>%group_by(genome)%>%mutate(protein_count = n())%>%filter(protein.count.full >= i)%>%sample_n(i)


      ##calucute contig scores
      aaRes <- aaTmp%>%group_by(genome)%>%summarise(archaeal_score = sum(score))%>%
        mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))

      paRes <- paTmp%>%group_by(genome)%>%summarise(phage_score = sum(score))%>%
        mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))

      eaRes <- eaTmp%>%group_by(genome)%>%summarise(eukaryotic_score = sum(score))%>%
        mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))


      tmpRes <- aaRes%>%full_join(paRes, by = "genome")%>%full_join(eaRes, by = "genome")

      tmpRes <- tmpRes%>%
        mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))%>%
        mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))%>%
        mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))%>%mutate(call = ifelse(archaeal_score >= phage_score,
                                                                                                            ifelse(archaeal_score >= eukaryotic_score,
                                                                                                                   ifelse(archaeal_score > 60, "archaeal", "none"),
                                                                                                                   ifelse(eukaryotic_score > 150, "eukaryotic", "none")),
                                                                                                            ifelse(phage_score >= eukaryotic_score, ifelse(phage_score > 110, "phage", "none"), ifelse(eukaryotic_score > 150, "eukaryotic", "none"))))
      tmpRes <- tmpRes%>%mutate(protein_count = i)


      res <- res%>%bind_rows(tmpRes)

    }
  }
  for(j in 1:3){
    for(i in 26:50){
      printRemaining(i = ((j-1)*25 + i-25), length = 75, increment = 5)

      ##select i proteins from genomes
      aaTmp <- aa1%>%group_by(genome)%>%mutate(protein_count = n())%>%filter(protein.count.full >= i)%>%sample_n(i)
      paTmp <- pa1%>%group_by(genome)%>%mutate(protein_count = n())%>%filter(protein.count.full >= i)%>%sample_n(i)
      eaTmp <- ea1%>%group_by(genome)%>%mutate(protein_count = n())%>%filter(protein.count.full >= i)%>%sample_n(i)


      ##calucute contig scores
      aaRes <- aaTmp%>%group_by(genome)%>%summarise(archaeal_score = sum(score))%>%
        mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))

      paRes <- paTmp%>%group_by(genome)%>%summarise(phage_score = sum(score))%>%
        mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))

      eaRes <- eaTmp%>%group_by(genome)%>%summarise(eukaryotic_score = sum(score))%>%
        mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))


      tmpRes <- aaRes%>%full_join(paRes, by = "genome")%>%full_join(eaRes, by = "genome")

      tmpRes <- tmpRes%>%
        mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))%>%
        mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))%>%
        mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))%>%
        mutate(call = ifelse(archaeal_score >= phage_score, ifelse(archaeal_score >= eukaryotic_score, ifelse(archaeal_score > 0, "archaeal", "none"), "eukaryotic"), ifelse(phage_score >= eukaryotic_score, "phage", "eukaryotic")))

      tmpRes <- tmpRes%>%mutate(protein_count = i)


      res <- res%>%bind_rows(tmpRes)

    }
  }




  runningTime <- proc.time() - ptm1
  printRunningTime(runningTime = runningTime, type = "The simulation")
  return(res)
}


getwd()
setwd("/Volumes/userdata/student_users/thomasnicholson/PredVirusHost/")
#####

##read in the hmm output ####
aa <- read.table("arVOG_archaeal_test.txt")
ap <- read.table("arVOG_phage_test.txt")
ae <- read.table("arVOG_eukaryotic_test.txt")

pa <- read.table("baPOG_archaeal_test.txt")
pp <- read.table("baPOG_phage_test.txt")
pe <- read.table("baPOG_eukaryotic_test.txt")

ea <- read.table("euVOG_archaeal_test.txt")
ep <- read.table("euVOG_phage_test.txt")
ee <- read.table("euVOG_eukaryotic_test.txt")

hmmColnames <- c("target.name", "accession", "query.name", "accession.2", "E.value", "score", "bias", "E.value.2", "score.2", "bias.2", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description.of.target")


colnames(aa) <- hmmColnames
colnames(ae) <- hmmColnames
colnames(ap) <- hmmColnames

colnames(pa) <- hmmColnames
colnames(pe) <- hmmColnames
colnames(pp) <- hmmColnames

colnames(ea) <- hmmColnames
colnames(ep) <- hmmColnames
colnames(ee) <- hmmColnames
##get genome names matching to protein names ####
archaeal_lookup <- read.table("~/DB/viral/archaeal_virus_lookup.txt", sep = "\t", header = F, comment.char = "", quote = "")
phage_lookup <- read.table("~/DB/viral/phage_lookup.txt", sep = "\t", header = F, comment.char = "", quote = "")

#>I dont think this is correct
#eukaryotic_lookup <- read.table("~/DB/viral/eukaryotic_virus_proteinname_lookup.txt", sep = "\t", header = F, comment.char = "", quote = "")
#>> this is probably the right one
eukaryotic_lookup <- read.table("~/DB/viral/eukaryotic_lookup.txt", sep = "\t", header = F, comment.char = "", quote = "", fill = T)

phage_lookup <- phage_lookup%>%
  separate(col = V2, into = c("protein_description", "genome"), extra = "merge", remove = T, sep = "\\[")%>%
  separate(col = genome, into = c("i1", "i2"), extra = "warn", remove = F, sep = "\\[")%>%
  mutate(genome = ifelse(!is.na(i2), i2, genome))%>%
  select(-i1, -i2)

##get protein count for the test set ####
archaealProteins <- read.table("~/PredVirusHost/archaeal_virus_test_protein_list.txt")
phageProteins <- read.table("~/PredVirusHost/phage_test_protein_list.txt")
eukaryoticProteins <- read.table("~/PredVirusHost/eukaryotic_virus_test_protein_list.txt")

archaealProteins <- archaealProteins%>%mutate(test = T)%>%dplyr::rename(target.name = V1)
phageProteins <- phageProteins%>%mutate(test = T)%>%dplyr::rename(target.name = V1)
eukaryoticProteins <- eukaryoticProteins%>%mutate(test = T)%>%dplyr::rename(target.name = V1)


archaeal_lookup <- archaeal_lookup%>%dplyr::rename(target.name = V1, genome = V3)%>%group_by(genome)%>%full_join(archaealProteins)%>%filter(!is.na(test))%>%mutate(protein.count.full = n())
eukaryotic_lookup <- eukaryotic_lookup%>%dplyr::rename(target.name = V2, genome = V1)%>%group_by(genome)%>%full_join(eukaryoticProteins)%>%filter(!is.na(test))%>%mutate(protein.count.full = n())
phage_lookup <- phage_lookup%>%dplyr::rename(target.name = V1)%>%group_by(genome)%>%group_by(genome)%>%full_join(phageProteins)%>%filter(!is.na(test))%>%mutate(protein.count.full = n())


##get model scores ####
archaealModels <- read.table("/Volumes/userdata/student_users/thomasnicholson/complete_packages/predvirushost/archaeal_model_scores.txt", sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)
phageModels <- read.table("/Volumes/userdata/student_users/thomasnicholson/complete_packages/predvirushost/phage_model_scores.txt", sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)
eukaryoticModels <- read.table("/Volumes/userdata/student_users/thomasnicholson/complete_packages/predvirushost/eukaryotic_model_scores.txt", sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)

##simulate small contigs ####
##include the models scores

archaealRes <- simulateSmallContigs(arVOGres = aa, baPOGres = pa, euVOGres = ea, lookup = archaeal_lookup, discriminant_models_only = F)
table(archaealRes$call)
archaealRes <- read.table("~/PredVirusHost/archaeal_simulation.txt", sep = "", header = T, as.is = T)


archaealRes <- archaealRes%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))%>%mutate(call = ifelse(archaeal_score >= phage_score,
                                                                                                      ifelse(archaeal_score >= eukaryotic_score,
                                                                                                             ifelse(archaeal_score > 60, "archaeal", "none"),
                                                                                                             ifelse(eukaryotic_score > 150, "eukaryotic", "none")),
                                                                                                      ifelse(phage_score >= eukaryotic_score, ifelse(phage_score > 110, "phage", "none"), ifelse(eukaryotic_score > 150, "eukaryotic", "none"))))


write.table(x = archaealRes, "~/PredVirusHost/archaeal_simulation_discriminant_only.txt", quote = F, row.names = F, sep = "\t")

phageRes <- simulateSmallContigs(arVOGres = ap, baPOGres = pp, euVOGres = ep, lookup = phage_lookup, discriminant_models_only = T)
table(phagelRes$call)
phageRes <- read.table("~/PredVirusHost/phage_simulation.txt", sep = "", header = T, as.is = T)

phageRes <- phageRes%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))%>%mutate(call = ifelse(archaeal_score >= phage_score,
                                                                                                         ifelse(archaeal_score >= eukaryotic_score,
                                                                                                                ifelse(archaeal_score > 60, "archaeal", "none"),
                                                                                                                ifelse(eukaryotic_score > 150, "eukaryotic", "none")),
                                                                                                         ifelse(phage_score >= eukaryotic_score, ifelse(phage_score > 110, "phage", "none"), ifelse(eukaryotic_score > 150, "eukaryotic", "none"))))

write.table(x = phagelRes, "~/PredVirusHost/phage_simulation.txt", quote = F, row.names = F, sep = "\t")

eukaryoticRes <- simulateSmallContigs(arVOGres = ae, baPOGres = pe, euVOGres = ee, lookup = eukaryotic_lookup, discriminant_models_only = T)
table(eukaryoticRes$call)
eukaryoticRes <- read.table("~/PredVirusHost/eukaryotic_simulation.txt", sep = " ", header = T, as.is = T, fill = T)
eukaryoticRes <- eukaryoticRes%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))%>%mutate(call = ifelse(archaeal_score >= phage_score,
                                                                                                         ifelse(archaeal_score >= eukaryotic_score,
                                                                                                                ifelse(archaeal_score > 60, "archaeal", "none"),
                                                                                                                ifelse(eukaryotic_score > 150, "eukaryotic", "none")),
                                                                                                         ifelse(phage_score >= eukaryotic_score, ifelse(phage_score > 110, "phage", "none"), ifelse(eukaryotic_score > 150, "eukaryotic", "none"))))


write.table(x = eukaryoticRes, "~/PredVirusHost/eukaryotic_simulation.txt", quote = F, row.names = F, sep = "\t")








##roc curve ####

eukaryoticRes <- read.table(file = "~/PredVirusHost/eukaryotic_simulation.txt", sep = "\t", header = T)
archaealRes <- read.table(file = "~/PredVirusHost/archaeal_simulation.txt", sep = "\t", header = T)
phagelRes <- read.table(file = "~/PredVirusHost/phage_simulation.txt", sep = "\t", header = T)

## set max score
archaeal <- archaealRes%>%mutate(max.score = ifelse(call == "archaeal", archaeal_score, ifelse(call == "phage", phage_score, eukaryotic_score)))%>% 
  mutate(group = "archaeal")
phage <- phagelRes%>%mutate(max.score = ifelse(call == "archaeal", archaeal_score, ifelse(call == "phage", phage_score, eukaryotic_score)))%>% 
  mutate(group = "phage")
eukaryotic <- eukaryoticRes%>%mutate(max.score = ifelse(call == "archaeal", archaeal_score, ifelse(call == "phage", phage_score, eukaryotic_score))) %>% 
  mutate(group = "eukaryotic")

df <- data.frame(score = seq(0,10000, 10),
                 aTP = rep(NA, 1001), aFP = rep(NA, 1001),
                 pTP = rep(NA, 1001), pFP = rep(NA, 1001),
                 eTP = rep(NA, 1001), eFP = rep(NA, 1001),
                 eaTP = rep(NA, 1001), eaFP = rep(NA, 1001))

i <- 0
##loop through scores
##archaeal results
for(i in 0:10000){
 # printRemaining(i = i, length = 10000, increment = 5)
  ##work on increments of 10
  if( i %% 10 == 0 ){

    ##keep results where the max score is above i
    aLocal <- archaeal%>%filter(archaeal_score >= i)#%>%filter(max.score < i + 10)
    pLocal <- phage%>%filter(archaeal_score >= i)#%>%filter(max.score < i + 100)
    eLocal <- eukaryotic%>%filter(archaeal_score >= i)#%>%filter(max.score < i + 100)

    ##calculate and store TP and FP for each host
    df[(i/10 +1),2] <- nrow(aLocal)/nrow(archaeal)
    df[(i/10 +1),3] <- (nrow(pLocal) + nrow(eLocal))/(nrow(phage) + nrow(eukaryotic))

  }
}
##phage
for(i in 0:10000){
 # printRemaining(i = i, length = 10000, increment = 5)

  ##work on increments of 10
  if( i %% 10 == 0 ){

    ##keep results where the max score is above i
    aLocal <- archaeal%>%filter(phage_score >= i)#%>%filter(max.score < i + 10)
    pLocal <- phage%>%filter(phage_score >= i)#%>%filter(max.score < i + 100)
    eLocal <- eukaryotic%>%filter(phage_score >= i)#%>%filter(max.score < i + 100)

    ##calculate and store TP and FP for each host
    df[(i/10 +1),4] <- nrow(pLocal)/nrow(phage)
    df[(i/10 +1),5] <- (nrow(aLocal) + nrow(eLocal))/(nrow(archaeal) + nrow(eukaryotic))

  }
}
##eukaryotic
for(i in 0:10000){
  #printRemaining(i = i, length = 10000, increment = 5)

  ##work on increments of 10
  if( i %% 10 == 0 ){

    ##keep results where the max score is above i
    aLocal <- archaeal%>%filter(eukaryotic_score >= i)#%>%filter(max.score < i + 10)
    pLocal <- phage%>%filter(eukaryotic_score >= i)#%>%filter(max.score < i + 100)
    eLocal <- eukaryotic%>%filter(eukaryotic_score >= i)#%>%filter(max.score < i + 100)

    ##calculate and store TP and FP for each host
    df[(i/10 +1),6] <- nrow(eLocal)/nrow(eukaryotic)
    df[(i/10 +1),7] <- (nrow(aLocal) + nrow(pLocal))/(nrow(archaeal) + nrow(phage))

  }
}
##eukaryotic vs archaeal only
for(i in 0:10000){
 # printRemaining(i = i, length = 10000, increment = 5)

  ##work on increments of 10
  if( i %% 10 == 0 ){

    ##keep results where the max score is above i
    aLocal <- archaeal%>%filter(eukaryotic_score >= i)#%>%filter(max.score < i + 10)
    eLocal <- eukaryotic%>%filter(eukaryotic_score >= i)#%>%filter(max.score < i + 100)

    ##calculate and store TP and FP for each host
    df[(i/10 +1),8] <- nrow(eLocal)/nrow(eukaryotic)
    df[(i/10 +1),9] <- nrow(aLocal)/nrow(archaeal)

  }
}

all <- archaeal %>% bind_rows(phage) %>% bind_rows(eukaryotic)



rocData <- archaeal %>% mutate(response = ifelse(group == "Random Data", 0, ifelse(group == "Known SRA predicted", 1, NA))) %>% 
  filter(!is.na(response))

roc.curve(response = rocData$response, predicted = rocData$score_2, 
          main="ROC curve for Read Depths Scores")


##plot FP vs TP
ggplot() +
  geom_line(data = df, aes(x = pFP, y = pTP), colour = "Blue") +
  geom_line(data = df, aes(x = aFP, y = aTP), colour = "Red") +
  geom_line(data = df, aes(x = eFP, y = eTP), colour = "Green") +
  geom_line(data = df, aes(x = eaFP, y = eaTP), colour = "Yellow") +
#  geom_vline(xintercept = 0.01) +
  geom_vline(xintercept = 0.02) +
#  geom_vline(xintercept = 0.03) +
#  geom_vline(xintercept = 0.04) +
#  geom_vline(xintercept = 0.05) +
  xlim(0,1) +
  ylim(0,1)


df2 <- df%>%mutate(x = ifelse(eaFP >= 0.02, 0.02, 0))





## Results ####
archaealTab <- archaealRes%>%group_by(call)%>%summarise(archaea = n())%>%filter(!is.na(call))
phageTab <- phageRes%>%group_by(call)%>%summarise(phage = n())%>%filter(!is.na(call))
eukaryoticTab <- eukaryoticRes%>%group_by(call)%>%summarise(eukarya = n())%>%filter(!is.na(call))

combinedTab <- archaealTab%>%full_join(phageTab)%>%full_join(eukaryoticTab)

rownames(combinedTab) <- combinedTab$call




##This is specificity
spec <- c(sum(combinedTab[2:4, 3:4], na.rm = T)/(sum(combinedTab[1,3:4], na.rm = T) + sum(combinedTab[2:4, 3:4], na.rm = T)),
          sum(combinedTab[c(1,3,4), c(2,4)], na.rm = T)/(sum(combinedTab[2,c(2,4)], na.rm = T) + sum(combinedTab[c(1,3,4), c(2,4)], na.rm = T)),
          sum(combinedTab[c(1,2,4), 2:3], na.rm = T)/(sum(combinedTab[3,2:3], na.rm = T) + sum(combinedTab[c(1,2,4), 2:3], na.rm = T)))

##This is sensitivity
sens <- c(archaealTab[archaealTab$call == "archaeal",2]/sum(archaealTab$archaea),
          phageTab[phageTab$call == "phage",2]/sum(phageTab$phage),
          eukaryoticTab[eukaryoticTab$call == "eukaryotic",2]/sum(eukaryoticTab$eukarya))

tab <- data.frame(host = c("Archaeal Virus", "Bacteriophage", "Eukaryotic Virus"),
                  sensitivity = t(as.data.frame(sens)),
                  specificity = spec)

tt <- tab%>%select(host, sensitivity)%>%mutate(Result = rep("Sensitivity", 3))%>%dplyr::rename(score = sensitivity)%>%bind_rows(
  tab%>%select(host, specificity)%>%mutate(Result = rep("Specificity", 3))%>%dplyr::rename(score = specificity)
)

ggplot(data = tt, aes(x = host, y = score, fill = Result)) +
  geom_col(position = "dodge")

## Analysis of model specificity #####

##archaeal
aa <- read.table("arVOG_archaeal_train.txt")
ap <- read.table("arVOG_phage_train.txt")
ae <- read.table("arVOG_eukaryotic_train.txt")

hmmColnames <- c("target.name", "accession", "query.name", "accession.2", "E.value", "score", "bias", "E.value.2", "score.2", "bias.2", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description.of.target")


colnames(aa) <- hmmColnames
colnames(ae) <- hmmColnames
colnames(ap) <- hmmColnames

aa <- aa%>%filter(as.numeric(E.value) < 1e-5)
ap <- ap%>%filter(as.numeric(E.value) < 1e-5)
ae <- ae%>%filter(as.numeric(E.value) < 1e-5)

aaSum <- aa%>%group_by(query.name)%>%summarise(archaeal_bit_score_mean = mean(score))
apSum <- ap%>%group_by(query.name)%>%summarise(phage_bit_score_mean = mean(score))
aeSum <- ae%>%group_by(query.name)%>%summarise(eukaryotic_bit_score_mean = mean(score))

archaealModels <- aaSum%>%full_join(apSum)%>%full_join(aeSum)

archaealModels <- archaealModels%>%
  mutate(archaeal_bit_score_mean = ifelse(is.na(archaeal_bit_score_mean), 0, archaeal_bit_score_mean))%>%
  mutate(phage_bit_score_mean = ifelse(is.na(phage_bit_score_mean), 0, phage_bit_score_mean))%>%
  mutate(eukaryotic_bit_score_mean = ifelse(is.na(eukaryotic_bit_score_mean), 0, eukaryotic_bit_score_mean))%>%
  mutate(score.percentage = ifelse(phage_bit_score_mean >= eukaryotic_bit_score_mean,
                                   round((archaeal_bit_score_mean - phage_bit_score_mean)/archaeal_bit_score_mean, 2),
                                   round((archaeal_bit_score_mean - eukaryotic_bit_score_mean)/archaeal_bit_score_mean, 2)))%>%
  mutate(score.percentage = ifelse(score.percentage < 0, 0, score.percentage))


write.table(archaealModels%>%select(query.name,score.percentage), "archaeal_model_scores.txt", quote = F, col.names = T, row.names = F, sep = "\t")



##phage
pa <- read.table("baPOG_archaeal_train.txt")
pp <- read.table("baPOG_phage_train.txt")
pe <- read.table("baPOG_eukaryotic_train.txt")

hmmColnames <- c("target.name", "accession", "query.name", "accession.2", "E.value", "score", "bias", "E.value.2", "score.2", "bias.2", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description.of.target")


colnames(pa) <- hmmColnames
colnames(pe) <- hmmColnames
colnames(pp) <- hmmColnames

pa <- pa%>%filter(as.numeric(E.value) < 1e-5)
pp <- pp%>%filter(as.numeric(E.value) < 1e-5)
pe <- pe%>%filter(as.numeric(E.value) < 1e-5)

paSum <- pa%>%group_by(query.name)%>%summarise(archaeal_bit_score_mean = mean(score))
ppSum <- pp%>%group_by(query.name)%>%summarise(phage_bit_score_mean = mean(score))
peSum <- pe%>%group_by(query.name)%>%summarise(eukaryotic_bit_score_mean = mean(score))

phageModels <- paSum%>%full_join(ppSum)%>%full_join(peSum)

phageModels <- phageModels%>%
  mutate(archaeal_bit_score_mean = ifelse(is.na(archaeal_bit_score_mean), 0, archaeal_bit_score_mean))%>%
  mutate(phage_bit_score_mean = ifelse(is.na(phage_bit_score_mean), 0, phage_bit_score_mean))%>%
  mutate(eukaryotic_bit_score_mean = ifelse(is.na(eukaryotic_bit_score_mean), 0, eukaryotic_bit_score_mean))%>%
  mutate(score.percentage = ifelse(archaeal_bit_score_mean >= eukaryotic_bit_score_mean,
                                   round((phage_bit_score_mean - archaeal_bit_score_mean)/phage_bit_score_mean, 2),
                                   round((phage_bit_score_mean - eukaryotic_bit_score_mean)/phage_bit_score_mean, 2)))%>%
  mutate(score.percentage = ifelse(score.percentage < 0, 0, score.percentage))


write.table(phageModels%>%select(query.name,score.percentage), "phage_model_scores.txt", quote = F, col.names = T, row.names = F, sep = "\t")



##eukaryotic
ea <- read.table("euVOG_archaeal_train.txt")
ep <- read.table("euVOG_phage_train.txt")
ee <- read.table("euVOG_eukaryotic_train.txt")

hmmColnames <- c("target.name", "accession", "query.name", "accession.2", "E.value", "score", "bias", "E.value.2", "score.2", "bias.2", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description.of.target")


colnames(ea) <- hmmColnames
colnames(ee) <- hmmColnames
colnames(ep) <- hmmColnames

ea <- ea%>%filter(as.numeric(E.value) < 1e-5)
ep <- ep%>%filter(as.numeric(E.value) < 1e-5)
ee <- ee%>%filter(as.numeric(E.value) < 1e-5)

eaSum <- ea%>%group_by(query.name)%>%summarise(archaeal_bit_score_mean = mean(score))
epSum <- ep%>%group_by(query.name)%>%summarise(phage_bit_score_mean = mean(score))
eeSum <- ee%>%group_by(query.name)%>%summarise(eukaryotic_bit_score_mean = mean(score))

eukaryoticModels <- eaSum%>%full_join(epSum)%>%full_join(eeSum)

eukaryoticModels <- eukaryoticModels%>%
  mutate(archaeal_bit_score_mean = ifelse(is.na(archaeal_bit_score_mean), 0, archaeal_bit_score_mean))%>%
  mutate(phage_bit_score_mean = ifelse(is.na(phage_bit_score_mean), 0, phage_bit_score_mean))%>%
  mutate(eukaryotic_bit_score_mean = ifelse(is.na(eukaryotic_bit_score_mean), 0, eukaryotic_bit_score_mean))%>%
  mutate(score.percentage = ifelse(phage_bit_score_mean >= archaeal_bit_score_mean,
                                   round((eukaryotic_bit_score_mean - phage_bit_score_mean)/eukaryotic_bit_score_mean, 2),
                                   round((eukaryotic_bit_score_mean - archaeal_bit_score_mean)/eukaryotic_bit_score_mean, 2)))%>%
  mutate(score.percentage = ifelse(score.percentage < 0, 0, score.percentage))


write.table(eukaryoticModels%>%select(query.name,score.percentage), "eukaryotic_model_scores.txt", quote = F, col.names = T, row.names = F, sep = "\t")


##protein results#####

proteinSens <- data.frame(protein_count = 0, sensitivity = 0, count = 0)
tmpRes <- archaealRes#%>%mutate(protein_count = mround(protein_count, 5))%>%mutate(protein_count = ifelse(is.na(protein_count), 0, protein_count))
i <- 20
for(i in 1:(max(tmpRes$protein_count, na.rm = T))){
  print(i)
  tmp <- as.data.frame(tmpRes%>%filter(protein_count == i))#%>%filter(!is.na(call))
  tmp2 <- as.data.frame(tmp%>%group_by(call)%>%summarise(count = n())%>%filter(!is.na(call)))


  sens  <- (nrow(tmp) - sum(as.numeric(tmp2[tmp2$call != "archaeal", 2])))/nrow(tmp)*100



  if(length(sens) == 1){
    proteinSens <- proteinSens%>%bind_rows(data.frame(protein_count = i, sensitivity = sens, count = nrow(tmp%>%filter(!is.na(call)))))
  }else(
    proteinSens <- proteinSens%>%bind_rows(data.frame(protein_count = i, sensitivity = NA, count = nrow(tmp%>%filter(!is.na(call)))))
  )
}

#proteinSens[3,2] <- 94.56237


proteinSpec <- data.frame(protein_count = 0, specificity = 0, count = 0)
tmpRes <- phageRes#%>%mutate(protein_count = mround(protein_count, 5))%>%mutate(protein_count = ifelse(is.na(protein_count), 0, protein_count))

i <- 4
for(i in 1:(max(tmpRes$protein_count, na.rm = T))){
  print(i)
  tmp <- as.data.frame(tmpRes%>%filter(protein_count == i))%>%filter(!is.na(call))
  tmp2 <- as.data.frame(tmp%>%group_by(call)%>%summarise(count = n())%>%filter(!is.na(call)))
  spec  <- (nrow(tmp) - as.numeric(tmp2[tmp2$call == "archaeal", 2]))/nrow(tmp)*100
  if(length(spec) == 1){
    proteinSpec <- proteinSpec%>%bind_rows(data.frame(protein_count = i, specificity = spec, count = nrow(tmp%>%filter(!is.na(call)))))
  }else(
    proteinSpec <- proteinSpec%>%bind_rows(data.frame(protein_count = i, specificity = NA, count = nrow(tmp%>%filter(!is.na(call)))))
  )
}

proteinSens <- proteinSens%>%mutate(sensitivity = ifelse(is.na(sensitivity), 100, sensitivity))
proteinSpec <- proteinSpec%>%mutate(specificity = ifelse(is.na(specificity), 100, specificity))


ggplot() +
  #geom_point(data = proteinSens, aes(x = protein_count, y = sensitivity, colour="Sensitivity")) +
  geom_line(data = proteinSens, aes(x = protein_count, y = sensitivity, colour = "Sensitivity"), size = 1.5) +
  #geom_point(data = proteinSpec, aes(x = protein_count, y = specificity, colour="Specificity")) +
  geom_line(data = proteinSpec, aes(x = protein_count, y = specificity, colour = "Specificity"), size = 1.5) +
  coord_cartesian(ylim = c(0, 105), xlim = c(0, 28)) +
  labs(x = "Number of Proteins", y = "Percentage (Sensitivity and Specificity)") +
  theme_bw()

#####



scoreSens <- data.frame(score = 0, sensitivity = 0, count = 0)

tmpRes <- archaealRes%>%mutate(archaeal_score = mround(archaeal_score, 50))%>%mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))




i <- 2
for(i in 1:(max(tmpRes$archaeal_score)/50)){
  print(i*50)
  tmp <- as.data.frame(tmpRes%>%filter(archaeal_score == i*50))
  tmp2 <- as.data.frame(tmp%>%group_by(call)%>%summarise(count = n())%>%filter(!is.na(call)))
  sens  <- as.numeric(tmp2[tmp2$call == "archaeal", 2])/nrow(tmp)*100
  if(length(sens) == 1){
    scoreSens <- scoreSens%>%bind_rows(data.frame(score = i*50, sensitivity = sens, count = nrow(tmp)))
  }else(
    scoreSens <- scoreSens%>%bind_rows(data.frame(score = i*50, sensitivity = NA, count = nrow(tmp)))
  )
}


scoreSpec <- data.frame(score = 0, specificity = 0, count = 0)
tmpRes <- phageRes%>%mutate(phage_score = mround(phage_score, 10))%>%mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))

i <- 3
for(i in 1:(max(tmpRes$phage_score)/10)){
  print(i*10)
  tmp <- as.data.frame(tmpRes%>%filter(phage_score == i*10))
  tmp2 <- as.data.frame(tmp%>%group_by(call)%>%summarise(count = n())%>%filter(!is.na(call)))
  spec  <- (nrow(tmp) - as.numeric(tmp2[tmp2$call == "archaeal", 2]))/nrow(tmp)*100


  if(length(spec) == 1){
    scoreSpec <- scoreSpec%>%bind_rows(data.frame(score = i*10, specificity = spec, count = nrow(tmp)))
  }else(
    scoreSpec <- scoreSpec%>%bind_rows(data.frame(score = i*10, specificity = 100, count = nrow(tmp)))
  )
}




ggplot() +
  geom_histogram(data = archaealRes, aes(x = archaeal_score, y = ..count..), fill='grey', binwidth = 50) +
  #geom_point(data = scoreSens, aes(x = score, y = sensitivity, colour="Sensitivity")) +
  geom_line(data = scoreSens, aes(x = score, y = sensitivity, colour = "Sensitivity"), size = 1.5) +
  #geom_point(data = scoreSpec, aes(x = score, y = specificity, colour="Specificity")) +
  geom_line(data = scoreSpec, aes(x = score, y = specificity, colour = "Specificity"), size = 1.5) +
  coord_cartesian(ylim = c(0, 105), xlim = c(0, 5000)) +
  labs(x = "Score", y = "Percentage (Sensitivity and Specificity)") +
  theme_bw()


##Protein count data Metagenomes #####

imgvr <- read.table("~/PredVirusHost/IMGVR/IMGVR.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
bovine <- read.table("~/DB/metagenomes/bovine/bovine.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
mgm4441095 <- read.table("~/DB/metagenomes/mgm4441095.3.350.genecalling.coding.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
mgm4441096 <- read.table("~/DB/metagenomes/mgm4441096.3.350.genecalling.coding.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
mgm4442583 <- read.table("~/DB/metagenomes/mgm4442583.3.350.genecalling.coding.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
mgm4449206 <- read.table("~/DB/metagenomes/mgm4449206.3.350.genecalling.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
mgm4522044 <- read.table("~/DB/metagenomes/mgm4522044.3.350.genecalling.coding_bovine_rumen.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
mgm4522044Cluster <- read.table("~/DB/metagenomes/mgm4522044.3.550.cluster.aa90.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
mgm4529716 <- read.table("~/DB/metagenomes/mgm4529716.3.350.genecalling.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
mgm4529719 <- read.table("~/DB/metagenomes/mgm4529719.3.350.genecalling.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
mgm4530143 <- read.table("~/DB/metagenomes/mgm4530143.3.350.genecalling.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
mgm4544453 <- read.table("~/DB/metagenomes/mgm4544453.3.350.genecalling.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")
HG_ORFans <- read.table("~/DB/metagenomes/HG_ORFans.faa.tmp.folder/scores.txt", header = T, sep = "\t", fill = T, comment.char = "")

bovine <- bovine%>%mutate(genome = as.character(genome))
mgm4441095 <- mgm4441095%>%mutate(genome = as.character(genome))
mgm4441096 <- mgm4441096%>%mutate(genome = as.character(genome))
HG_ORFans <- HG_ORFans%>%mutate(genome = as.character(genome))
mgm4442583 <- mgm4442583%>%mutate(genome = as.character(genome))
mgm4449206 <- mgm4449206%>%mutate(genome = as.character(genome))
mgm4522044 <- mgm4522044%>%mutate(genome = as.character(genome))
mgm4522044Cluster <- mgm4522044Cluster%>%mutate(genome = as.character(genome))
mgm4529716 <- mgm4529716%>%mutate(genome = as.character(genome))
mgm4529719 <- mgm4529719%>%mutate(genome = as.character(genome))
mgm4530143 <- mgm4530143%>%mutate(genome = as.character(genome))
mgm4544453 <- mgm4544453%>%mutate(genome = as.character(genome))

metagenomes <- imgvr%>%
  bind_rows(bovine)%>%
  bind_rows(mgm4441095)%>%
  bind_rows(mgm4441096)%>%
  bind_rows(HG_ORFans)%>%
  bind_rows(mgm4442583)%>%
  bind_rows(mgm4449206)%>%
  bind_rows(mgm4522044)%>%
  bind_rows(mgm4522044Cluster)%>%
  bind_rows(mgm4529716)%>%
  bind_rows(mgm4529719)%>%
  bind_rows(mgm4530143)%>%
  bind_rows(mgm4544453)

hotpools <-  mgm4441095%>%
  bind_rows(mgm4441096)%>%
  bind_rows(mgm4442583)%>%
  bind_rows(mgm4449206)%>%
  bind_rows(mgm4529716)%>%
  bind_rows(mgm4529719)%>%
  bind_rows(mgm4530143)%>%
  bind_rows(mgm4544453)

dat <- data.frame(protein.count = c(1:2200), prop = rep(NA, 2200))
i <- 2
for(i in 1:2200){
  printRemaining(i, 2200, increment = 10)
  tmpMetagenomes <- metagenomes%>%
    filter(protein.counts >= i)%>%
    filter(protein.counts <= i + 2)
  if(nrow(tmpIMGVR) > 0){
    dat[i,2] <- 1- nrow(tmpMetagenomes%>%filter(call == "none"))/nrow(tmpMetagenomes)
  }else{
    dat[i,2] <- 0
  }
}


ggplot(data = dat) +
  geom_line( aes(x = protein.count, y = prop)) +
#  stat_smooth(aes(x = protein.count, y = prop), method = lm, formula = y ~ poly(x,3), se = FALSE) +
  xlim(0, 50)


metagenomes <- metagenomes%>%mutate(max.score = ifelse(call == "archaeal", archaeal_score, ifelse(call == "phage", phage_score, ifelse(call == "eukaryotic", eukaryotic_score, 0))))

meanMetagenomes <- metagenomes%>%group_by(protein.counts)%>%summarise(mean = median(max.score))
sdMetagenomes <- metagenomes%>%group_by(protein.counts)%>%summarise(sd = sd(max.score))
countsMetagenomes <- metagenomes%>%group_by(protein.counts)%>%summarise(count = n())

summaryMetagenomes <- meanMetagenomes%>%
  full_join(sdMetagenomes, by = "protein.counts")%>%
  full_join(countsMetagenomes, by = "protein.counts")

TS <- zoo::zoo(summaryMetagenomes$mean)
mean <- zoo::rollapply(TS, width = 5, by = 1, FUN = mean, align = "left")
tmp1 <- as.data.frame(mean)
df <- data.frame(mean = c(0,0,0,0))
tmp1 <- tmp1%>%bind_rows(df)

TS <- zoo::zoo(summaryMetagenomes$sd)
sd <- zoo::rollapply(TS, width = 5, by = 1, FUN = mean, align = "left")
tmp2 <- as.data.frame(sd)
df <- data.frame(sd = c(0,0,0,0))
tmp2 <- tmp2%>%bind_rows(df)




summaryMetagenomes <- summaryMetagenomes%>%select(-mean, -sd)%>%bind_cols(tmp1, tmp2)


summaryMetagenomes <- summaryMetagenomes%>%mutate(upper = mean + sd, lower = mean - sd)




ggplot(data = summaryMetagenomes%>%filter(protein.counts <= 50)) +
  geom_line(aes(x = protein.counts, y = mean)) +
  geom_line(aes(x = protein.counts, y = upper)) +
  geom_line(aes(x = protein.counts, y = lower))


##Protein count data Simulation #####
archaeal <- read.table("~/PredVirusHost/archaeal_simulation.txt", header = T, sep = "\t", fill = T, comment.char = "")
phage <- read.table("~/PredVirusHost/phage_simulation.txt", header = T, sep = "\t", fill = T, comment.char = "")
eukaryotic <- read.table("~/PredVirusHost/eukaryotic_simulation.txt", header = T, sep = "\t", fill = T, comment.char = "")


archaeal <- archaeal%>%mutate(genome = as.character(genome))
phage <- phage%>%mutate(genome = as.character(genome))
eukaryotic <- eukaryotic%>%mutate(genome = as.character(genome))


all <- archaeal%>%
  bind_rows(phage)%>%
  bind_rows(eukaryotic)



dat <- data.frame(protein.count = c(1:50), prop = rep(NA, 50))
i <- 2
for(i in 1:50){
  printRemaining(i, 50, increment = 10)
  tmpAll <- all%>%
    filter(protein_count >= i)%>%
    filter(protein_count <= i + 2)
  if(nrow(tmpIMGVR) > 0){
    dat[i,2] <- 1- nrow(tmpAll%>%filter(call == "none"))/nrow(tmpAll)
  }else{
    dat[i,2] <- 0
  }
}


ggplot(data = dat) +
  geom_line( aes(x = protein.count, y = prop)) +
  #  stat_smooth(aes(x = protein.count, y = prop), method = lm, formula = y ~ poly(x,3), se = FALSE) +
  xlim(0, 50)


all <- all%>%mutate(max.score = ifelse(call == "archaeal", archaeal_score, ifelse(call == "phage", phage_score, ifelse(call == "eukaryotic", eukaryotic_score, 0))))

meanAll <- all%>%group_by(protein_count)%>%summarise(mean = median(max.score))
sdAll <- all%>%group_by(protein_count)%>%summarise(sd = sd(max.score))
countsAll <- all%>%group_by(protein_count)%>%summarise(count = n())

summaryAll <- meanAll%>%
  full_join(sdAll, by = "protein_count")%>%
  full_join(countsAll, by = "protein_count")

summaryAll <- summaryAll%>%mutate(upper = mean + sd, lower = mean - sd)




ggplot(data = summaryAll%>%filter(protein_count <= 50)) +
  geom_line(aes(x = protein_count, y = mean)) +
  geom_line(aes(x = protein_count, y = upper)) +
  geom_line(aes(x = protein_count, y = lower))
