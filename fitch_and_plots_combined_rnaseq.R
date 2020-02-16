
# packages ----------------------------------------------------------------

suppressMessages(library(tidyverse))
library(dplyr)
suppressMessages(library(tjnFunctions))
library(VennDiagram)

# functions ---------------------------------------------------------------


plotKnownvsConserved <- function(dat, columns, not_zero = F){
  dat <- dat%>%mutate(conserved = F)
if(not_zero){
  for(i in 1:nrow(dat)){
    dat[i, ncol(dat)] <- ("1" %in% dat[i, columns])
    if(dat[i, ncol(dat)] == F){
    dat[i, ncol(dat)] <- ("0-1" %in% dat[i, columns])
    }

  }
}else{
  for(i in 1:nrow(dat)){
    dat[i, ncol(dat)] <- ("1" %in% dat[i, columns])
  }
}


  conservedSet <- dat%>%filter(conserved)
  knownSet <- dat%>%filter(new_feature == F)

  vennSet <- conservedSet%>%bind_rows(knownSet)%>%unique()



  area1 <- nrow(subset(vennSet, conserved == T))
  area2 <- nrow(subset(vennSet, new_feature == F))
  cross.area <- nrow(subset(vennSet, new_feature == F & conserved == T))

  grid.newpage()
  draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                     scaled = T,
                     #cat.default.pos= "text",
                     #cat.pos = c(-50, 50),
                     #category = c("Conserved and Expressed", "Known")
                     category = c("", "")
  )
}



# old work ----------------------------------------------------------------


#----#### Need to manually specifiy all the files and the taxonomy
setwd("~/phd/RNASeq/combined_gff_files/")
genusOrder <- c("GCA_000017745-GCA_000017765_merged.gff",
                "GCA_000438825-GCA_000747565_merged.gff",
                "GCA_000007385-GCA_001042875_merged.gff",
                "GCA_000017745-GCA_000747565_merged.gff",
                "GCA_000017745-GCA_000438825_merged.gff",
                "GCA_000017745-GCA_000007385_merged.gff",
                "GCA_000017745-GCA_001042875.1_merged.gff")

#Merge all the files
##escherichia

mergedData <- read.table(genusOrder[1], header = T, sep = "\t", comment.char = "", quote = "", as.is = T)

workingDat <- mergedData%>%
  mutate(merge1 = ifelse(nchar(id1) > 0, 1, 0))%>%
  mutate(merge2 = ifelse(nchar(id2) > 0, ifelse(prop_overlap > 0.5, 1 , 0), 0))%>%
  rename(E1 = id1, E2 = id2, E1found = merge1, E2found = merge2, new_feature_E1 = new_feature)%>%
select(E1, E2, E1found, E2found, new_feature_E1)


refDat1 <- workingDat%>%mutate(E1 = ifelse(nchar(E1) > 0, E1, paste("E1_1", row_number(), sep = "_")),
                               E2 = ifelse(nchar(E2) > 0, E2, paste("E2_1", row_number(), sep = "_")))





##serratia

mergedData <- read.table(genusOrder[2], header = T, sep = "\t", comment.char = "", quote = "", as.is = T)

workingDat <- mergedData%>%
  mutate(merge1 = ifelse(nchar(id1) > 0, 1, 0))%>%
  mutate(merge2 = ifelse(nchar(id2) > 0, ifelse(prop_overlap > 0.5, 1 , 0), 0))%>%
  rename(S1 = id1, S2 = id2, S1found = merge1, S2found = merge2, new_feature_S1 = new_feature)%>%
select(S1, S2, S1found, S2found, new_feature_S1)



refDat2 <- workingDat%>%mutate(S1 = ifelse(nchar(S1) > 0, S1, paste("S1_1", row_number(), sep = "_")),
                               S2 = ifelse(nchar(S2) > 0, S2, paste("S2_1", row_number(), sep = "_")))



##xanthomonas

mergedData <- read.table(genusOrder[3], header = T, sep = "\t", comment.char = "", quote = "", as.is = T)

workingDat <- mergedData%>%
  mutate(merge1 = ifelse(nchar(id1) > 0, 1, 0))%>%
  mutate(merge2 = ifelse(nchar(id2) > 0, ifelse(prop_overlap > 0.5, 1 , 0), 0))%>%
  rename(X1 = id1, X2 = id2, X1found = merge1, X2found = merge2, new_feature_X1 = new_feature)%>%
  select(X1, X2, X1found, X2found, new_feature_X1)


refDat4 <- workingDat%>%mutate(X1 = ifelse(nchar(X1) > 0, X1, paste("X1_1", row_number(), sep = "_")),
                               X2 = ifelse(nchar(X2) > 0, X2, paste("X2_1", row_number(), sep = "_")))







##trying to join escherichia and serratia
mergedData <- read.table(genusOrder[5], header = T, sep = "\t", comment.char = "", quote = "", as.is = T)

workingDat <- mergedData%>%
  mutate(merge1 = ifelse(nchar(id1) > 0, 1, 0))%>%
  mutate(merge2 = ifelse(nchar(id2) > 0, ifelse(prop_overlap > 0.5, 1 , 0), 0))%>%
  rename(E1 = id1, S1 = id2, E1found = merge1, E1S1 = merge2, new_feature_E2 = new_feature)%>%
select(E1, S1, E1found, E1S1, new_feature_E2)


refDat3 <- workingDat%>%mutate(E1 = ifelse(nchar(E1) > 0, E1, paste("E1_2", row_number(), sep = "_")),
                               S1 = ifelse(nchar(S1) > 0, S1, paste("S1_2", row_number(), sep = "_")))




##trying to join escherichia and xanthomonas
mergedData <- read.table(genusOrder[6], header = T, sep = "\t", comment.char = "", quote = "", as.is = T)

workingDat <- mergedData%>%
  mutate(merge1 = ifelse(nchar(id1) > 0, 1, 0))%>%
  mutate(merge2 = ifelse(nchar(id2) > 0, ifelse(prop_overlap > 0.5, 1 , 0), 0))%>%
  rename(E1 = id1, X1 = id2, E1found = merge1, E1X1 = merge2, new_feature_E3 = new_feature)%>%
  select(E1, X1, E1found, E1X1, new_feature_E3)


refDat5 <- workingDat%>%mutate(E1 = ifelse(nchar(E1) > 0, E1, paste("E1_3", row_number(), sep = "_")),
                               X1 = ifelse(nchar(X1) > 0, X1, paste("X1_2", row_number(), sep = "_")))






dat <- refDat1%>%
  full_join(refDat3%>%select(-E1found), by = "E1")%>%
  full_join(refDat2, by = "S1")%>%
  full_join(refDat5%>%select(-E1found), by = "E1")%>%
  full_join(refDat4, by = "X1")

dat <- dat%>%
  mutate(E1S1 = ifelse(is.na(E1S1), 0, E1S1),
         E2found = ifelse(is.na(E2found), 0, E2found),
         E1found = ifelse(is.na(E1found), 0, E1found),
         S1found = ifelse(is.na(S1found), 0, S1found),
         X1found = ifelse(is.na(X1found), 0, X1found),
         S2found = ifelse(is.na(S2found), 0, S2found),
         X2found = ifelse(is.na(X2found), 0, X2found),
         E1X1 = ifelse(is.na(E1X1), 0, E1X1))

dat <- dat%>%
  mutate(E = "")%>%
  mutate(S = "")%>%
  mutate(X = "")%>%
  mutate(ES = "")%>%
  mutate(ESX = "")

dat <- dat%>%select(E1, E2, E1found, E2found, S1, E1S1, S2, S1found, S2found, X1, E1X1, X2, X1found, X2found, E, S, X, ES, ESX,
             new_feature_E1, new_feature_E2, new_feature_S1, new_feature_X1, new_feature_E3)

i <- 1
library(UpSetR)

dat <- dat%>%
  mutate(E1 = ifelse(nchar(E1) > 0, E1, paste("E1_4", row_number(), sep = "_")),
         E2 = ifelse(nchar(E2) > 0, E2, paste("E2_2", row_number(), sep = "_")),
         S1 = ifelse(nchar(S1) > 0, S1, paste("S1_3", row_number(), sep = "_")),
         S2 = ifelse(nchar(S2) > 0, S2, paste("S2_2", row_number(), sep = "_")),
         X1 = ifelse(nchar(X1) > 0, X1, paste("X1_3", row_number(), sep = "_")),
         X2 = ifelse(nchar(X2) > 0, X2, paste("X2_2", row_number(), sep = "_")))%>%
  mutate(E1 = ifelse(!is.na(E1), E1, paste("E1_4", row_number(), sep = "_")),
         E2 = ifelse(!is.na(E2), E2, paste("E2_2", row_number(), sep = "_")),
         S1 = ifelse(!is.na(S1), S1, paste("S1_3", row_number(), sep = "_")),
         S2 = ifelse(!is.na(S2), S2, paste("S2_2", row_number(), sep = "_")),
         X1 = ifelse(!is.na(X1), X1, paste("X1_3", row_number(), sep = "_")),
         X2 = ifelse(!is.na(X2), X2, paste("X2_2", row_number(), sep = "_")))

dat <- dat%>%mutate(new_feature = ifelse(!is.na(new_feature_E1), new_feature_E1,
                                         ifelse(!is.na(new_feature_E2), new_feature_E2,
                                                ifelse(!is.na(new_feature_S1), new_feature_S1,
                                                       ifelse(!is.na(new_feature_X1), new_feature_X1,
                                                              ifelse(!is.na(new_feature_E3), new_feature_E3, NA))))))


tmp <- dat%>%select(E1found, E2found, S1found, S2found, X1found, X2found)%>%
  rename(E1 = E1found, E2 = E2found, S1 = S1found, S2 = S2found, X1 = X1found, X2 = X2found)
upset(tmp, sets = c("E1", "E2", "S1", "S2", "X1",
                    "X2"), mb.ratio = c(0.55, 0.45), order.by = "freq")

dat <- dat%>%mutate(conserved = ifelse(E == 1, T,
                                       ifelse(S == 1, T,
                                              ifelse(ES == 1, T,
                                                     ifelse(X == 1, T,
                                                            ifelse(ESX == 1, T, F))))))

library(VennDiagram)
conservedSet <- dat%>%filter(conserved)
knownSet <- dat%>%filter(new_feature == F)

vennSet <- conservedSet%>%bind_rows(knownSet)%>%unique()

nrow(subset(vennSet, conserved == T))
nrow(subset(vennSet, new_feature == F))
nrow(subset(vennSet, new_feature == F & conserved == T))

grid.newpage()
draw.pairwise.venn(area1 = 640, area2 = 1122, cross.area = 412, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)
draw.pairwise.venn(area1 = 640, area2 = 1122, cross.area = 412, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)
draw.pairwise.venn(area1 = 640, area2 = 1122, cross.area = 412, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)



conservedSetE <- dat%>%filter(conserved, E != "0")
knownSetE <- dat%>%filter(new_feature == F, E != "0")

vennSetE <- conservedSetE%>%bind_rows(knownSetE)%>%unique()

area1 <- nrow(subset(vennSetE, conserved == T))
area2 <- nrow(subset(vennSetE, new_feature == F))
cross.area <- nrow(subset(vennSetE, new_feature == F & conserved == T))

grid.newpage()
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)

conservedSetS <- dat%>%filter(conserved, S  != "0")
knownSetS <- dat%>%filter(new_feature == F, S != "0")

vennSetS <- conservedSetS%>%bind_rows(knownSetS)%>%unique()

area1 <- nrow(subset(vennSetS, conserved == T))
area2 <- nrow(subset(vennSetS, new_feature == F))
cross.area <- nrow(subset(vennSetS, new_feature == F & conserved == T))

grid.newpage()
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)

conservedSetX <- dat%>%filter(conserved, X  != "0")
knownSetX <- dat%>%filter(new_feature == F, X != "0")

vennSetX <- conservedSetX%>%bind_rows(knownSetX)%>%unique()

area1 <- nrow(subset(vennSetX, conserved == T))
area2 <- nrow(subset(vennSetX, new_feature == F))
cross.area <- nrow(subset(vennSetX, new_feature == F & conserved == T))

grid.newpage()
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)

conservedSetES <- dat%>%filter(conserved, ES  == "1")
knownSetES <- dat%>%filter(new_feature == F, ES != "0")

vennSetES <- conservedSetES%>%bind_rows(knownSetES)%>%unique()%>%mutate(conserved = ifelse(ES == 1, T, F))

area1 <- nrow(subset(vennSetES, conserved == T))
area2 <- nrow(subset(vennSetES, new_feature == F))
cross.area <- nrow(subset(vennSetES, new_feature == F & conserved == T))

grid.newpage()
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)
draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross.area, fill = c("blue", "red"),
                   scaled = T,
                   #cat.default.pos= "text",
                   #cat.pos = c(-50, 50),
                   #category = c("Conserved and Expressed", "Known")
                   category = c("", "")
)





#Fitch algorithm calcuations ####
for(i in 1:nrow(dat)){

  if(length(intersect(dat[i,3], dat[i,4])) > 0){
  E <- intersect(dat[i,3], dat[i,4])
  }else{
    E <- union(dat[i,3], dat[i,4])
  }
  if(length(intersect(dat[i,8], dat[i,9])) > 0){
    S <- intersect(dat[i,8], dat[i,9])
  }else{
    S <- union(dat[i,8], dat[i,9])
  }

  if(length(intersect(dat[i,13], dat[i,14])) > 0){
    X <- intersect(dat[i,13], dat[i,14])
  }else{
    X <- union(dat[i,13], dat[i,14])
  }



  if(length(intersect(E, S)) > 0){
      ES <- intersect(E, S)
  }else{
    ES <- union(E, S)
  }

  if(length(intersect(ES, X)) > 0){
    ESX <- intersect(ES, X)
  }else{
    ESX <- union(ES, X)
  }
  dat[i, 15] <- paste(E, collapse = ", ")
  dat[i, 16] <- paste(S, collapse = ", ")
  dat[i, 17] <- paste(X, collapse = ", ")
  dat[i, 18] <- paste(ES, collapse = ", ")
  dat[i, 19] <- paste(ESX, collapse = ", ")
  }


# UpsetR plots ------------------------------------------------------------

dat <- read.table("~/phd/RNASeq/combined_gff_files/escherichia-shigella-salmonella-enterobacter-klebsiella-serratia_merged.gff", header = T, sep = "\t",
                  comment.char = "", quote = "", as.is = T )

##get list of files and individual ids
dat_files <- unlist(strsplit(dat$file_id[1], "-"))
dat_ids <- dat$id

##create matrix listing all ids in rows and the files for columns
ids <- matrix(nrow = length(dat_ids), ncol = (length(dat_files)+ 1))
colnames(ids) <- c("id", dat_files)
ids[,1] <- dat_ids
#ids <- as.data.frame(ids)
i <- 1
j <- 1

##test yes/no (1/NA) if an ID was found in a file
for(i in 1:nrow(ids)){
printRemaining(length = nrow(ids),i = i)
    uid <- ids[i,1]
  uid <- unlist(strsplit(uid, "-"))
  for(j in 1:length(uid)){
    x <- paste(unlist(strsplit(uid[j], "_"))[1:2], collapse = "_")
    y <- unlist(strsplit(uid[j], "_"))[3]
    if(y != "0"){
    val <- match(x , colnames(ids))
    ids[i, val] <- 1
    }
  }
}

ids[is.na(ids)] <- 0

ids <- as.data.frame(ids)

##set the values as numeric (as.data.frame messed this up)
for(i in 2:ncol(ids)){
  ids[,i] <- as.numeric(as.character(ids[,i]))
}

colnames(ids)

##make an upsetR plot using all columns
UpSetR::upset(ids, sets = colnames(ids)[2:ncol(ids)], mb.ratio = c(0.55, 0.45), order.by = 'freq',
              keep.order = T)


##select the genera level column from the dat dataframe
genera <- dat%>%select(serratia.fitch, shigella.fitch, salmonella.fitch,
                       escherichia.fitch, klebsiella, enterobacter)
#genera[genera == "1-1"] <- 1

##set any value that is not 0 to 1
genera[genera != 0] <- 1

colnames(genera)[1:2] <- c("serratia", "enterobacter")

##set the values as numeric
for(i in 1:ncol(genera)){
  genera[,i] <- as.numeric(as.character(genera[,i]))
}

##select rows where atleast one of the genera contained a value
generaSmall <- genera%>%
  mutate(count =  serratia + enterobacter + salmonella + escherichia + klebsiella)%>%
  filter(count > 1)%>%
  select(-count)


##make an upsetR plot using all genera values
UpSetR::upset(genera, sets = colnames(genera), mb.ratio = c(0.55, 0.45), order.by = "freq")

##make an upsetR plot using genera values
UpSetR::upset(generaSmall, sets = colnames(generaSmall), mb.ratio = c(0.55, 0.45), order.by = "freq")

#

# Get sequences -----------------------------------------------------------
dat <- read.table("~/phd/RNASeq/enterics-serratia_merged.gff", header = T, sep = "\t",
                  comment.char = "", quote = "", as.is = T )


dat_files <- unlist(strsplit(dat$file_id[1], "-"))
dat_ids <- dat$id

ids_lookup <- matrix(nrow = length(dat_ids), ncol = (length(dat_files)+ 1))
colnames(ids_lookup) <- c("id", dat_files)
ids_lookup[,1] <- dat_ids
#ids <- as.data.frame(ids)
i <- 1
j <- 1
for(i in 1:nrow(ids_lookup)){
  uid <- ids_lookup[i,1]
  uid <- unlist(strsplit(uid, "-"))
  for(j in 1:length(uid)){
    x <- paste(unlist(strsplit(uid[j], "_"))[1:2], collapse = "_")
    y <- unlist(strsplit(uid[j], "_"))[3]
    if(y != "0"){
      val <- match(x , colnames(ids))
      ids_lookup[i, val] <- y
    }

  }
}


ids_lookup[is.na(ids_lookup)] <- ""

ids_lookup <- as.data.frame(ids_lookup)

for(i in 2:ncol(ids_lookup)){
  ids_lookup[,i] <- as.character(ids_lookup[,i])
}

i <- 1
j <- 5
ids_seq <- ids_lookup%>%full_join(dat%>%select(set_val, id), by = "id")%>%filter(set_val != "0")
for(i in 1:nrow(ids_seq)){
  for(j in 2:(ncol(ids_seq) - 1)){


    if(ids_seq[i,j] != ""){
      x <- system(command = paste("sed -n '", (as.numeric(ids_seq[i, j]) + 1), "p'",
                         " ~/phd/RNASeq/new_calls/",
                         colnames(ids_seq)[j],
                         "_new_calls.txt", sep = ""), intern = T)
      x <- unlist(strsplit(as.character(x), split = "\t"))
      #print(x[3:4])


      ids_seq[i, j] <- system(paste("test_string=`grep -v ^'>' ~/phd/RNASeq/sequences/",
                         colnames(ids_seq)[j],
                         ".fna | tr -d '\n'`; echo ${test_string:",
                         x[3], ":", as.numeric(x[4]) - as.numeric(x[3]),
                         "}", sep = ""), intern = T)


    }
  }
}

ids_seq <- ids_seq%>%mutate(ids_short = paste("sra_enterics-serratia_", row_number(), sep = ""))

for(i in 1:nrow(ids_seq)){
  for(j in 2:(ncol(ids_seq) -  2)){
    if(ids_seq[i, j] != ""){
      mat <- matrix(ncol = 1, nrow = 2)
      mat[1,1] <- paste(">", colnames(ids_seq)[j], sep = "")
      mat[2,1] <- ids_seq[i,j]
      write.table(mat, file = paste("~/phd/RNASeq/sRNAs/", ids_seq[i, ncol(ids_seq)], ".fasta", sep = ""),
                  row.names = F, col.names = F, quote = F, sep = "\t", append = T)

    }
  }
}



#
##phylogenetic distance vs prop conserved ####

dat <- read.table("~/phd/RNASeq/enterics-serratia_merged.gff", header = T, sep = "\t",
                  comment.char = "", quote = "", as.is = T )

colnames(dat)

dat[dat == "1-1"] <- "1"
dat[dat == "1-0"] <- "0-1"




dat_1 <- dat%>%select(escherichia)
dat_2 <- dat%>%select(klebsiella,
                      GCA_000438825.1.GCA_000747565.1,
                      GCA_001874505.1.GCA_002303275.1,
                      salmonella)
dat_3 <- dat%>%select(escherichia.salmonella, enterobacter.klebsiella)
dat_4 <- dat%>%select(enterics)
dat_5 <- dat%>%select(enterics.serratia)


dat_1 <- data.frame(set_val = as.vector(as.matrix(dat_1[,])))
dat_2 <- data.frame(set_val = as.vector(as.matrix(dat_2[,])))
dat_3 <- data.frame(set_val = as.vector(as.matrix(dat_3[,])))
dat_4 <- data.frame(set_val = as.vector(as.matrix(dat_4[,])))
dat_5 <- data.frame(set_val = as.vector(as.matrix(dat_5[,])))

# dat_1 <- dat_1%>%group_by(set_val)%>%summarise(count = n())%>%mutate(group = 1)
# dat_2 <- dat_2%>%group_by(set_val)%>%summarise(count = n())%>%mutate(group = 2)
# dat_3 <- dat_3%>%group_by(set_val)%>%summarise(count = n())%>%mutate(group = 3)
# dat_4 <- dat_4%>%group_by(set_val)%>%summarise(count = n())%>%mutate(group = 4)
# dat_5 <- dat_5%>%group_by(set_val)%>%summarise(count = n())%>%mutate(group = 5)

dat_1 <- dat_1%>%mutate(group = 1)
dat_2 <- dat_2%>%mutate(group = 2)
dat_3 <- dat_3%>%mutate(group = 3)
dat_4 <- dat_4%>%mutate(group = 4)
dat_5 <- dat_5%>%mutate(group = 5)


dat_all <- dat_1%>%
  bind_rows(dat_2)%>%
  bind_rows(dat_3)%>%
  bind_rows(dat_4)%>%
  bind_rows(dat_5)%>%
  mutate(is_1 = ifelse(set_val == "1", 1, 0))%>%
  mutate(not_0 = ifelse(set_val != "0" && set != "1", 1, 0))


ggplot() +
  geom_bar(data = dat_all%>%filter(set_val == "1"), aes(x = group, y = (..count.. / nrow(dat_all%>%filter(set_val == "1")))))+
  geom_bar(data = dat_all%>%filter(set_val != "0"), aes(x = group, y = (..count.. / nrow(dat_all%>%filter(set_val != "1")))))


ggplot() +
  geom_bar(data = dat_all%>%filter(set_val != "0"), aes(x = group, y = ..count.., group = set_val, fill = set_val), position = 'dodge')

#
## taxonomic levels vs prop conserved####

##read data
dat <- read.table("~/phd/RNASeq/enterics-serratia_merged.gff", header = T, sep = "\t",
                  comment.char = "", quote = "", as.is = T )
##make sets consistent
dat[dat == "1-1"] <- "1"
dat[dat == "1-0"] <- "0-1"

##
dat_files <- unlist(strsplit(dat$file_id[1], "-"))
head(dat_files)

dat_ids <- dat$id

ids_lookup <- data.frame(id = "", row = "0", genome = "")

for(i in 1:length(dat_ids)){
  uid <- dat_ids[i]
  uid <- unlist(strsplit(uid, "-"))
  for(j in 1:length(uid)){
    x <- paste(unlist(strsplit(uid[j], "_"))[1:2], collapse = "_")
    y <- unlist(strsplit(uid[j], "_"))[3]
    if(y != "0"){
      tmp <- data.frame(id = dat_ids[i], row = y, genome = x)
      ids_lookup <- ids_lookup%>%bind_rows(tmp)
    }

  }
}


datSmall <- dat[,c(8, 10, 12, 14, 15, 17:ncol(dat))]
ids_lookup <- ids_lookup%>%filter(row != "0")%>%left_join(datSmall, by = "id")


##classify each sRNA by the highest level of conservation found

main_col <- 7
genera_col <- c(14, 19, 20, 22)
species_col <- c(13:15, 18:22)
any_col <- c(7:ncol(ids_lookup))

assignConservationLevel <- function(ids_lookup, main_col = 7, genera_col, species_col, any_col = c(7:ncol(ids_lookup))){
  ids_lookup <- ids_lookup%>%mutate(type = "")
  for(i in 1:nrow(ids_lookup)){
    if("1" %in% ids_lookup[i, main_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "main_conserved"
    }else if("0-1" %in% ids_lookup[i, main_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "main_0-1"
    }else if("1" %in% ids_lookup[i, genera_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "genera_conserved"
    }else if("0-1" %in% ids_lookup[i, genera_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "genera_0-1"
    }else if("1" %in% ids_lookup[i, species_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "species_conserved"
    }else if("0-1" %in% ids_lookup[i, species_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "species_0-1"
    }else if("1" %in% ids_lookup[i, any_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "any_conserved"
    }else if("0-1" %in% ids_lookup[i, any_col]){
      ids_lookup[i, ncol(ids_lookup)] <- "any_0-1"
    }

  }
  return(ids_lookup)
}

ids_lookup <- assignConservationLevel(ids_lookup = ids_lookup, main_col = main_col, genera_col = genera_col,
                                      species_col = species_col, any_col = any_col)
id_count <- ids_lookup%>%select(id, type)%>%unique()%>%group_by(type)%>%summarise(count = n())
id_prop <- ids_lookup%>%group_by(type, new_feature)%>%
  summarise(count = n())%>%group_by(type)%>%mutate(prop = count /sum(count))%>%
  arrange(prop)%>%arrange(new_feature)%>%ungroup()

id_prop$type <- factor(id_prop$type, levels = unique(id_prop$type))

ggplot() +
  geom_col(data = id_prop, aes(x = type, y = prop, fill = new_feature))


#
#Venn Diagrams ####
dat <- read.table("~/phd/RNASeq/enterics-serratia_merged.gff", header = T, sep = "\t",
                  comment.char = "", quote = "", as.is = T )

colnames(dat)

dat[dat == "1-1"] <- "1"
dat[dat == "1-0"] <- "0-1"

plotKnownvsConserved(dat, columns = 15, not_zero = T)
plotKnownvsConserved(dat, columns = c(17:(ncol(dat) - 1)))
plotKnownvsConserved(dat, columns = c(23, 27:29, 31), not_zero = T)




#
##nwk format ordering ####

nwk <- readLines("~/Downloads/enterics_test.nwk")
nwkDat <- data.frame(lines = nwk, start = "", accession = "", group = "",stringsAsFactors = F)



groupNumber <- 0
groupLetter <- letters[1]
numberList <- c()
for(i in 1:nrow(nwkDat)){
  nwkDat[i,2] <- substr(as.character(nwkDat[i,1]), 1, 1)
  if(nwkDat[i, 2] == "("){
    groupNumber <- groupNumber + 1
  }

  if(nwkDat[i,2] == "'"){
    x <- substr(as.character(nwkDat[i,1]), 2, (nchar(nwkDat[i,1])))
    x <- unlist(strsplit(x, "'"))[1]
    if(groupNumber %in% numberList ==F){
      numberList <- c(numberList, groupNumber)
    }

    nwkDat[i,3] <- x
  }

}






#
##alifold results ####
filename <- readLines("~/phd/RNASeq/srna_all/new/main_conserved/known/alifold/main_conserved_known_alifold.txt")

alifoldDataFrameSetup <- function(filename, group, level, type, conserved = T){
lines <- readLines(filename)
alifold <- data.frame(lines = lines)
alifoldNames <- alifold%>%filter(substr(lines, 1, 1) == "#")%>%dplyr::rename(file = lines)
alifoldDat <- alifold%>%
  filter(lines != "")%>%
  filter(substr(lines, 1, 1) != "#")%>%
  filter(substr(lines, 1, 2) != " -")%>%
  filter(substr(lines, 1, 3) != "  F")

alifoldDat <- alifoldDat%>%
  mutate(lines = str_squish(lines))%>%
  separate(col = lines, into = c("From", "To", "Strand", "Native.MFE", "Mean.MFE", "STDV", "Z"),
                                    sep = " ", remove = T)
#alifoldDatFwd <- alifoldDat%>%filter(Strand == "+")%>%bind_cols(alifoldNames)
#alifoldDatRev <- alifoldDat%>%filter(Strand == "-")%>%bind_cols(alifoldNames)

#alifoldDat <- alifoldDatFwd%>%bind_rows(alifoldDatRev)%>%arrange(file)%>%
#  mutate(group = group)

alifoldDat <- alifoldDat%>%
  mutate(group = group,
         level = level,
         type = type,
         conserved = conserved)
return(alifoldDat)

}


mck <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/main_conserved/known/alifold/main_conserved_known_alifold.txt",
                             group = "mck",
                             level = "m",
                             type = "k")
gck <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/genera_conserved/known/alifold/genera_conserved_known_alifold.txt",
                             group = "gck",
                             level = "g",
                             type = "k")
gcn <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/genera_conserved/new/alifold/genera_conserved_new_alifold.txt",
                             group = "gcn",
                             level = "g",
                             type = "n")
gok <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/genera_0-1/known/alifold/genera_0-1_known_alifold.txt",
                             group = "gok",
                             level = "g",
                             type = "k",
                             conserved = F)
gon <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/genera_0-1/new/alifold/genera_0-1_new_alifold.txt",
                             group = "gon",
                             level = "g",
                             type = "n",
                             conserved = F)
sck <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/species_conserved/known/alifold/species_conserved_known_alifold.txt",
                             group = "sck",
                             level = "s",
                             type = "k")
scn <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/species_conserved/new/alifold/species_conserved_new_alifold.txt",
                             group = "scn",
                             level = "s",
                             type = "n")
sok <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/species_0-1/known/alifold/species_0-1_known_alifold.txt",
                             group = "sok",
                             level = "s",
                             type = "k",
                             conserved = F)
son <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/species_0-1/new/alifold/species_0-1_new_alifold.txt",
                             group = "son",
                             level = "s",
                             type = "n",
                             conserved = F)
aok <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/any_0-1/known/alifold/any_0-1_known_alifold.txt",
                             group = "aok",
                             level = "a",
                             type = "k",
                             conserved = F)
aon <- alifoldDataFrameSetup(filename = "~/phd/RNASeq/srna_all/new/any_0-1/new/alifold/any_0-1_new_alifold.txt",
                             group = "aon",
                             level = "a",
                             type = "n",
                             conserved = F)

alifoldDat <- mck%>%bind_rows(gck, gcn, gok, gon, sck, scn, sok, son, aok, aon)


ggplot(alifoldDat%>%filter(type == "k")) +
  geom_count(aes(x=as.numeric(Native.MFE), group = level, color = level), stat = "sum")



+
  stat_function(fun=dnorm, args = list(mean=0, sd=1), color="black")



+
  xlim(-10, 5)




ggplot(alifoldDat%>%filter(conserved == T)) +
  geom_count(aes(x=as.numeric(Z), group = group, color = type)) +
  stat_function(fun=dnorm, args = list(mean=0, sd=1), color="black") +
  xlim(-10, 5)


ggplot(alifoldDat) +
  geom_density(aes(x=as.numeric(Z), group = group, color = type)) +
  stat_function(fun=dnorm, args = list(mean=0, sd=1), color="black") +
  xlim(-10, 5)



ggplot(alifoldDat) +
  geom_density(aes(x=as.numeric(Native.MFE), group = group, color = type)) +
  xlim(-50, 1)


ggplot(alifoldDat) +
  geom_density(aes(x=as.numeric(Mean.MFE), group = group, color = group)) +
  xlim(-50, 1)

plot(as.numeric(alifoldDat$Native.MFE) ~ as.numeric(alifoldDat$Mean.MFE))



#
##other ####




library(limma)




hsb2 <- read.csv("https://stats.idre.ucla.edu/wp-content/uploads/2016/02/hsb2-3.csv")
attach(hsb2)
hw <- (write >= 60)
hm <- (math >= 60)
hr <- (read >= 60)
c3 <- cbind(hw, hm, hr)

a <- vennCounts(c3)
vennDiagram(a)



