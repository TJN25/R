library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
mround <- function(x,base){ 
  base*round(x/base) 
} 

## Base scoring script on this section #####
getwd()
setwd("/Volumes/userdata/student_users/thomasnicholson/PredVirusHost/")


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


archaeal_genome_info <- read.table("~/DB/viral/archaeal_virus_test_genome_info.txt", sep = " ", header = F, comment.char = "", quote = "", fill = T)

archaeal_genome_info <- archaeal_genome_info%>%
  dplyr::rename(protein_count = V1, genome = V2)%>%
  filter(!is.na(protein_count))

archaeal_lookup <- read.table("~/DB/viral/archaeal_virus_lookup.txt", sep = "\t", header = F, comment.char = "", quote = "")
phage_lookup <- read.table("~/DB/viral/phage_lookup.txt", sep = "\t", header = F, comment.char = "", quote = "")
eukaryotic_lookup <- read.table("~/DB/viral/eukaryotic_virus_proteinname_lookup.txt", sep = "\t", header = F, comment.char = "", quote = "")

phage_lookup <- phage_lookup%>%
  separate(col = V2, into = c("protein_description", "genome"), extra = "merge", remove = T, sep = "\\[")%>%
  separate(col = genome, into = c("i1", "i2"), extra = "warn", remove = F, sep = "\\[")%>%
  mutate(genome = ifelse(!is.na(i2), i2, genome))%>%
  select(-i1, -i2)

archaeal_lookup <- archaeal_lookup%>%dplyr::rename(target.name = V1, genome = V3)
eukaryotic_lookup <- eukaryotic_lookup%>%dplyr::rename(target.name = V1, genome = V2)
phage_lookup <- phage_lookup%>%dplyr::rename(target.name = V1)


archaealModels <- read.table("archaeal_model_scores.txt", sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)
phageModels <- read.table("phage_model_scores.txt", sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)
eukaryoticModels <- read.table("eukaryotic_model_scores.txt", sep = "\t", comment.char = "", quote = "", fill = T, as.is = T, header = T)

##archaeal ####
aa <- aa%>%left_join(archaeal_lookup)%>%filter(as.numeric(score) > 30)
pa <- pa%>%left_join(archaeal_lookup)%>%filter(as.numeric(score) > 30)
ea <- ea%>%left_join(archaeal_lookup)%>%filter(as.numeric(score) > 30)

aa <- aa%>%left_join(archaealModels)
pa <- pa%>%left_join(phageModels)
ea <- ea%>%left_join(eukaryoticModels)

aa <- aa%>%mutate(score = score*score.percentage)
pa <- pa%>%mutate(score = score*score.percentage)
ea <- ea%>%mutate(score = score*score.percentage)

aaRes <- aa%>%group_by(genome)%>%summarise(archaeal_score = sum(score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))

paRes <- pa%>%group_by(genome)%>%summarise(phage_score = sum(score))%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))

eaRes <- ea%>%group_by(genome)%>%summarise(eukaryotic_score = sum(score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))

archaealRes <- aaRes%>%full_join(paRes)%>%full_join(eaRes)

archaealRes <- archaealRes%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))%>%
  mutate(call = ifelse(archaeal_score >= phage_score, ifelse(archaeal_score >= eukaryotic_score, "archaeal", "eukaryotic"), ifelse(phage_score >= eukaryotic_score, "phage", "eukaryotic")))%>%
  full_join(archaeal_genome_info)

##phage ####
ap <- ap%>%left_join(phage_lookup)%>%filter(as.numeric(score) > 30)
pp <- pp%>%left_join(phage_lookup)%>%filter(as.numeric(score) > 30)
ep <- ep%>%left_join(phage_lookup)%>%filter(as.numeric(score) > 30)


ap <- ap%>%left_join(archaealModels)
pp <- pp%>%left_join(phageModels)
ep <- ep%>%left_join(eukaryoticModels)

ap <- ap%>%mutate(score = score*score.percentage)
pp <- pp%>%mutate(score = score*score.percentage)
ep <- ep%>%mutate(score = score*score.percentage)

apRes <- ap%>%group_by(genome)%>%summarise(archaeal_score = sum(score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))

ppRes <- pp%>%group_by(genome)%>%summarise(phage_score = sum(score))%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))

epRes <- ep%>%group_by(genome)%>%summarise(eukaryotic_score = sum(score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))

phageRes <- apRes%>%full_join(ppRes)%>%full_join(epRes)

phageRes <- phageRes%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))%>%
  mutate(call = ifelse(archaeal_score >= phage_score, ifelse(archaeal_score >= eukaryotic_score, "archaeal", "eukaryotic"), ifelse(phage_score >= eukaryotic_score, "phage", "eukaryotic")))

##eukaryotic ####
ae <- ae%>%separate(col = target.name, into = c("protein", "genome"), sep = "\\|", remove = F, extra = "drop")%>%filter(as.numeric(score) > 30)
pe <- pe%>%separate(col = target.name, into = c("protein", "genome"), sep = "\\|", remove = F, extra = "drop")%>%filter(as.numeric(score) > 30)
ee <- ee%>%separate(col = target.name, into = c("protein", "genome"), sep = "\\|", remove = F, extra = "drop")%>%filter(as.numeric(score) > 30)

ae <- ae%>%left_join(archaealModels)
pe <- pe%>%left_join(phageModels)
ee <- ee%>%left_join(eukaryoticModels)

ae <- ae%>%mutate(score = score*score.percentage)
pe <- pe%>%mutate(score = score*score.percentage)
ee <- ee%>%mutate(score = score*score.percentage)

aeRes <- ae%>%group_by(genome)%>%summarise(archaeal_score = sum(score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))

peRes <- pe%>%group_by(genome)%>%summarise(phage_score = sum(score))%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))

eeRes <- ee%>%group_by(genome)%>%summarise(eukaryotic_score = sum(score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))

eukaryoticRes <- aeRes%>%full_join(peRes)%>%full_join(eeRes)

eukaryoticRes <- eukaryoticRes%>%
  mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))%>%
  mutate(archaeal_score = ifelse(is.na(archaeal_score), 0, archaeal_score))%>%
  mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))%>%
  mutate(call = ifelse(archaeal_score >= phage_score, ifelse(archaeal_score >= eukaryotic_score, "archaeal", "eukaryotic"), ifelse(phage_score >= eukaryotic_score, "phage", "eukaryotic")))



##


archaealRes <- archaealRes%>%
  mutate(difference.percentage = ifelse(call == "archaeal", 
                                        ifelse(phage_score >= eukaryotic_score, 
                                               round((archaeal_score - phage_score)/archaeal_score*100, 2), 
                                               round((archaeal_score - eukaryotic_score)/archaeal_score*100, 2)
                                               ),
                                        ifelse(call == "phage",ifelse(archaeal_score >= eukaryotic_score, 
                                                                     round((phage_score - archaeal_score)/phage_score*100, 2), 
                                                                     round((phage_score - eukaryotic_score)/phage_score*100, 2)
                                        ), ifelse(archaeal_score >= phage_score, 
                                                  round((eukaryotic_score - archaeal_score)/eukaryotic_score*100, 2), 
                                                  round((eukaryotic_score - phage_score)/eukaryotic_score*100, 2)
                                        ))
                                        
                                        ))%>%
  mutate(call = ifelse(archaeal_score == 0, ifelse(phage_score == 0, ifelse(eukaryotic_score == 0, "none", call),call),call))





phageRes <- phageRes%>%
  mutate(difference.percentage = ifelse(call == "archaeal", 
                                        ifelse(phage_score >= eukaryotic_score, 
                                               round((archaeal_score - phage_score)/archaeal_score*100, 2), 
                                               round((archaeal_score - eukaryotic_score)/archaeal_score*100, 2)
                                        ),
                                        ifelse(call == "phage",ifelse(archaeal_score >= eukaryotic_score, 
                                                                      round((phage_score - archaeal_score)/phage_score*100, 2), 
                                                                      round((phage_score - eukaryotic_score)/phage_score*100, 2)
                                        ), ifelse(archaeal_score >= phage_score, 
                                                  round((eukaryotic_score - archaeal_score)/eukaryotic_score*100, 2), 
                                                  round((eukaryotic_score - phage_score)/eukaryotic_score*100, 2)
                                        ))
                                        
  ))%>%
  mutate(call = ifelse(archaeal_score == 0, ifelse(phage_score == 0, ifelse(eukaryotic_score == 0, "none", call),call),call))


eukaryoticRes <- eukaryoticRes%>%
  mutate(difference.percentage = ifelse(call == "archaeal", 
                                        ifelse(phage_score >= eukaryotic_score, 
                                               round((archaeal_score - phage_score)/archaeal_score*100, 2), 
                                               round((archaeal_score - eukaryotic_score)/archaeal_score*100, 2)
                                        ),
                                        ifelse(call == "phage",ifelse(archaeal_score >= eukaryotic_score, 
                                                                      round((phage_score - archaeal_score)/phage_score*100, 2), 
                                                                      round((phage_score - eukaryotic_score)/phage_score*100, 2)
                                        ), ifelse(archaeal_score >= phage_score, 
                                                  round((eukaryotic_score - archaeal_score)/eukaryotic_score*100, 2), 
                                                  round((eukaryotic_score - phage_score)/eukaryotic_score*100, 2)
                                        ))
                                        
  ))%>%
  mutate(call = ifelse(archaeal_score == 0, ifelse(phage_score == 0, ifelse(eukaryotic_score == 0, "none", call),call),call))

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

proteinSens <- data.frame(protein_count = 0, sensitivity = 35, count = 0)
tmpRes <- archaealRes%>%mutate(protein_count = mround(protein_count, 5))%>%mutate(protein_count = ifelse(is.na(protein_count), 0, protein_count))

for(i in 1:(max(tmpRes$protein_count)/5)){
  print(i) 
  tmp <- as.data.frame(tmpRes%>%filter(protein_count == i*5))
  tmp2 <- as.data.frame(tmp%>%group_by(call)%>%summarise(count = n())%>%filter(!is.na(call)))
  sens  <- as.numeric(tmp2[tmp2$call == "archaeal", 2])/nrow(tmp)*100
  if(length(sens) == 1){
    proteinSens <- proteinSens%>%bind_rows(data.frame(protein_count = i*5, sensitivity = sens, count = nrow(tmp))) 
  }else(
    proteinSens <- proteinSens%>%bind_rows(data.frame(protein_count = i*5, sensitivity = NA, count = nrow(tmp))) 
  )
}

#proteinSens[3,2] <- 94.56237


proteinSpec <- data.frame(protein_count = 0, specificity = 35, count = 0)
tmpRes <- phageRes%>%mutate(protein_count = mround(protein_count, 5))%>%mutate(protein_count = ifelse(is.na(protein_count), 0, protein_count))

for(i in 1:(max(tmpRes$protein_count)/5)){
  print(i) 
  tmp <- as.data.frame(tmpRes%>%filter(protein_count == i*5))
  tmp2 <- as.data.frame(tmp%>%group_by(call)%>%summarise(count = n())%>%filter(!is.na(call)))
  spec  <- as.numeric(tmp2[tmp2$call == "phage", 2])/nrow(tmp)*100
  if(length(spec) == 1){
    proteinSpec <- proteinSpec%>%bind_rows(data.frame(protein_count = i*5, specificity = spec, count = nrow(tmp))) 
  }else(
    proteinSpec <- proteinSpec%>%bind_rows(data.frame(protein_count = i*5, specificity = NA, count = nrow(tmp))) 
  )
}

proteinSens <- proteinSens%>%mutate(sensitivity = ifelse(is.na(sensitivity), 100, sensitivity))


ggplot() +
  geom_point(data = proteinSens, aes(x = protein_count, y = sensitivity, colour="Sensitivity")) +
  geom_line(data = proteinSens, aes(x = protein_count, y = sensitivity, colour = "Sensitivity")) +
  geom_point(data = proteinSpec, aes(x = protein_count, y = specificity, colour="Specificity")) +
  geom_line(data = proteinSpec, aes(x = protein_count, y = specificity, colour = "Specificity")) +
  coord_cartesian(ylim = c(40, 105), xlim = c(0, 70)) +
  labs(x = "Number of Proteins", y = "Percentage (Sensitivity and Specificity)") +
  theme_bw()

#####



scoreSens <- data.frame(score = 0, sensitivity = 35, count = 0)

tmpRes <- phageRes%>%mutate(phage_score = mround(phage_score, 200))%>%mutate(phage_score = ifelse(is.na(phage_score), 0, phage_score))





for(i in 1:(max(tmpRes$phage_score)/200)){
  print(i*200) 
  tmp <- as.data.frame(tmpRes%>%filter(phage_score == i*200))
  tmp2 <- as.data.frame(tmp%>%group_by(call)%>%summarise(count = n())%>%filter(!is.na(call)))
  sens  <- as.numeric(tmp2[tmp2$call == "phage", 2])/nrow(tmp)*100
  if(length(sens) == 1){
    scoreSens <- scoreSens%>%bind_rows(data.frame(score = i*200, sensitivity = sens, count = nrow(tmp)))
  }else(
    scoreSens <- scoreSens%>%bind_rows(data.frame(score = i*200, sensitivity = NA, count = nrow(tmp)))
  )
}


scoreSpec <- data.frame(score = 0, specificity = 35, count = 0)
tmpRes <- eukaryoticRes%>%mutate(eukaryotic_score = mround(eukaryotic_score, 200))%>%mutate(eukaryotic_score = ifelse(is.na(eukaryotic_score), 0, eukaryotic_score))


for(i in 1:(max(tmpRes$eukaryotic_score)/200)){
  print(i*200) 
  tmp <- as.data.frame(tmpRes%>%filter(eukaryotic_score == i*200))
  tmp2 <- as.data.frame(tmp%>%group_by(call)%>%summarise(count = n())%>%filter(!is.na(call)))
  spec  <- as.numeric(tmp2[tmp2$call != "phage", 2])/nrow(tmp)*100
  if(length(spec) == 1){
    scoreSpec <- scoreSpec%>%bind_rows(data.frame(score = i*200, specificity = spec, count = nrow(tmp)))
  }else(
    scoreSpec <- scoreSpec%>%bind_rows(data.frame(score = i*200, specificity = NA, count = nrow(tmp)))
  )
}




ggplot() +
  geom_point(data = scoreSens, aes(x = score, y = sensitivity, colour="Sensitivity")) +
  geom_line(data = scoreSens, aes(x = score, y = sensitivity, colour = "Sensitivity")) +
  geom_point(data = scoreSpec, aes(x = score, y = specificity, colour="Specificity")) +
  geom_line(data = scoreSpec, aes(x = score, y = specificity, colour = "Specificity")) +
  coord_cartesian(ylim = c(40, 105), xlim = c(0, 3000)) +
  labs(x = "Score", y = "Percentage (Sensitivity and Specificity)") +
  theme_bw()


#####