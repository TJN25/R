library(tidyverse)
library(reshape2)
library(UpSetR)
library(ggpubr)

dat <- read.table("/Volumes/2TTJN/Priming_Required_files/cas_arrays_genomes_refseq_83.txt", 
                  sep = "\t", 
                  quote = "", 
                  comment.char = "", 
                  fill = T, header = T)


dat <- dat %>% filter(domainName != "", !is.na(subtype.present))

dat <- dat %>% select(assembly_accession, refseq_category, taxid, species_taxid, organism_name, assembly_level,
                      domainName, protein.count, cas.proteins, cas1, cas2, basic.genes, subtype.list, subtype.count, single.system, contig.count, proteins.present, arrays.present, subtype.present, array.and_subtype)


genes <- read.table("/Volumes/2TTJN/Priming_Required_files/cas_genes_refseq_83.txt", 
                    sep = "\t", 
                    quote = "", 
                    comment.char = "", 
                    fill = T, header = T)


arrays <- read.table("/Volumes/2TTJN/Priming_Required_files/all_CRISPRDetect_refseq_83.TEST.fna", 
                    sep = " ", 
                    quote = "", 
                    comment.char = "", 
                    fill = T, header = F)


genes <- genes %>% separate(col = filepath, into = c(NA, "t1"), sep = "#", remove = F, extra = "merge")
genes <- genes %>% separate(col = t1, into = c("assembly_accession", NA), sep = "\\/", remove = F, extra = "drop")


arrays <- arrays %>% separate(col = V1, into = c("contig", "t2"), sep = "-", extra = "merge", remove = F)

arrays <- arrays %>% separate(col = t2, into = c("genome_name", "array_coords", "t3"), sep = "\\|", extra = "merge", remove = T)

arrays <- arrays %>% separate(col = t3, into = c("array_number", "spacer_number", "t4"), sep = "_", extra = "merge", remove = T)

arrays <- arrays %>% separate(col = t4, into = c("repeat_seq", "array_data", "spacer_seq"), sep = "\\[", extra = "merge", remove = T)



arrays <- data.frame(lapply(arrays, function(x) {
  gsub(">", "", x)
}))

arrays_high <- arrays %>% filter(grepl(pattern = "LOW", x = spacer_seq) == F)


arrays_per_contig <- arrays_high %>% select(contig, array_number) %>% unique() %>% group_by(contig) %>% 
  summarise(array_count = n())

spacers_per_contig <- arrays_high  %>% group_by(contig) %>% 
  summarise(spacer_count = n())

spacers_per_array <- arrays_high %>% select(contig, array_number, spacer_number) %>% group_by(contig, array_number) %>% 
  summarise(spacer_per_array_count = n())

spacers_per_array_min <- spacers_per_array %>% group_by(contig) %>% 
  summarise(min_spacers_per_array = min(spacer_per_array_count))
spacers_per_array_mean <- spacers_per_array %>% group_by(contig) %>% 
  summarise(mean_spacers_per_array = mean(spacer_per_array_count))
spacers_per_array_max <- spacers_per_array %>% group_by(contig) %>% 
  summarise(max_spacers_per_array = max(spacer_per_array_count))


arrays_summary <- arrays_per_contig %>% 
  full_join(spacers_per_contig, by = "contig") %>% 
  full_join(spacers_per_array_min, by = "contig") %>% 
  full_join(spacers_per_array_mean, by = "contig") %>% 
  full_join(spacers_per_array_max, by = "contig")






genes_small <- genes %>% select(contig, Gene_name, assembly_accession) %>% mutate(count = 1)

contig_genome_lookup <- genes %>% select(contig, assembly_accession) %>% unique()

mat <- reshape2::acast(data = genes_small , formula = contig ~ Gene_name)

mat <- as.data.frame(mat)

mat$contig <- rownames(mat)

mat <- mat %>% left_join(contig_genome_lookup, by = "contig")


mat <- mat %>% full_join(arrays_summary, by = "contig")


mat$adaptation_complex[mat$cas01 > 0 & mat$cas2 > 0] <- T
matStored <- mat
mat <- matStored
mat$adaptation_complex[is.na(mat$adaptation_complex)] <- F

mat$class_i[mat$cas7 > 0 | mat$cas7b > 0 | mat$cas7f > 0 | mat$csc2gr7 > 0 | mat$csf2gr7 > 0  | mat$csm3gr7 > 0  | mat$csm5gr7 > 0 | mat$cmr1gr7 > 0  | mat$cmr6gr7 > 0 | mat$cmr4gr7 > 0] <- T
mat$class_i[is.na(mat$class_i)] <- F


mat$class_ii[mat$cas9 > 0 | mat$cpf1 > 0 | mat$c2c1 > 0 | mat$c2c2 > 0] <- T
mat$class_ii[is.na(mat$class_ii)] <- F

mat$cas8_generic[mat$cas8a1 > 0 | mat$cas8a2 > 0 | mat$cas8a3 > 0  | mat$cas8a4 > 0  | mat$cas8a5 > 0 | mat$cas8a6 > 0 | mat$cas8a7 > 0 | mat$cas8a8 > 0 | mat$cas8b1 > 0 | mat$cas8b10 > 0 | mat$cas8b2 > 0 | mat$cas8b3 > 0 | mat$cas8b4 > 0 | mat$cas8b5 > 0 | mat$cas8b6 > 0  | mat$cas8b8 > 0 | mat$cas8b9 > 0 | mat$cas8c > 0 | mat$cas8e > 0 | mat$cas8f > 0 | mat$cas8u1 > 0 | mat$cas8u2 > 0 | mat$Cse1 > 0] <- T
mat$cas8_generic[is.na(mat$cas8_generic)] <- F

mat$cas3_generic[mat$cas3 > 0 | mat$cas3f > 0 | mat$cas3HD > 0] <- T
mat$cas3_generic[is.na(mat$cas3_generic)] <- F

mat$cas5_generic[mat$cas5 > 0 | mat$cas5a > 0 | mat$cas5f > 0 | mat$cas5u > 0 | mat$csc1gr5 > 0 | mat$csf3gr5 > 0 | mat$csm4gr5 > 0 | mat$csx10gr5 > 0 | mat$cmr3gr5 > 0] <- T
mat$cas5_generic[is.na(mat$cas5_generic)] <- F
mat$cas5_type_i[mat$cas5 > 0 | mat$cas5a > 0 | mat$cas5f > 0 | mat$cas5u > 0 | mat$csc1gr5 > 0] <- T
mat$cas5_type_i[is.na(mat$cas5_type_i)] <- F
mat$cas5_type_iii[mat$csm4gr5 > 0 | mat$csx10gr5 > 0 | mat$cmr3gr5 > 0] <- T
mat$cas5_type_iii[is.na(mat$cas5_type_iii)] <- F
mat$cas5_type_iv[mat$csf3gr5 > 0] <- T
mat$cas5_type_iv[is.na(mat$cas5_type_iv)] <- F



mat$cas6_generic[mat$cas6 > 0 | mat$cas6e > 0 | mat$cas6f > 0 | mat$csf5gr6 > 0] <- T
mat$cas6_generic[is.na(mat$cas6_generic)] <- F
mat$cas6_type_iv[mat$csf5gr6 > 0] <- T
mat$cas6_type_iv[is.na(mat$cas6_type_iv)] <- F

mat$cas7_generic[mat$cas7 > 0 | mat$cas7b > 0 | mat$cas7f > 0 | mat$csc2gr7 > 0 | mat$csf2gr7 > 0  | mat$csm3gr7 > 0  | mat$csm5gr7 > 0 | mat$cmr1gr7 > 0  | mat$cmr6gr7 > 0 | mat$cmr4gr7 > 0] <- T
mat$cas7_generic[is.na(mat$cas7_generic)] <- F
mat$cas7_type_i[mat$cas7 > 0 | mat$cas7b > 0 | mat$cas7f > 0 | mat$csc2gr7 > 0] <- T
mat$cas7_type_i[is.na(mat$cas7_type_i)] <- F
mat$cas7_type_iii[mat$csm3gr7 > 0  | mat$csm5gr7 > 0 | mat$cmr1gr7 > 0  | mat$cmr6gr7 > 0 | mat$cmr4gr7 > 0] <- T
mat$cas7_type_iii[is.na(mat$cas7_type_iii)] <- F
mat$cas7_type_iv[mat$csf2gr7 > 0] <- T
mat$cas7_type_iv[is.na(mat$cas7_type_iv)] <- F

mat$cas11_generic[mat$cas11b > 0 | mat$cas11d > 0 | mat$csa5gr11 > 0 | mat$Cse2 > 0 | mat$cse2gr11 > 0 | mat$csm2gr11 > 0 | mat$cmr5gr11 > 0] <- T
mat$cas11_generic[is.na(mat$cas11_generic)] <- F

mat$cas8a[mat$cas8a1 > 0 | mat$cas8a2 > 0 | mat$cas8a3 > 0  | mat$cas8a4 > 0  | mat$cas8a5 > 0 | mat$cas8a6 > 0 | mat$cas8a7 > 0 | mat$cas8a8 > 0 ] <- T
mat$cas8b[mat$cas8b1 > 0 | mat$cas8b10 > 0 | mat$cas8b2 > 0 | mat$cas8b3 > 0 | mat$cas8b4 > 0 | mat$cas8b5 > 0 | mat$cas8b6 > 0  | mat$cas8b8 > 0 | mat$cas8b9 > 0 | mat$cas8c > 0 | mat$cas8e > 0 | mat$cas8f > 0 | mat$cas8u1 > 0 | mat$cas8u2 > 0 | mat$Cse1 > 0] <- T
mat$cas8c[mat$cas8c > 0 ] <- T
mat$cas8e[mat$cas8e > 0 | mat$Cse1 > 0] <- T
mat$cas8f[mat$cas8f > 0] <- T
mat$cas8g[mat$cas8u1 > 0 | mat$cas8u2 > 0] <- T



mat$type_i_a_specific[mat$csa5gr11 > 0 | mat$cas8a1 > 0 |  mat$cas8a2 > 0 | mat$cas8a3 > 0 | mat$cas8a4 > 0  | mat$cas8a5 > 0 | mat$cas8a6 > 0 | mat$cas8a7 > 0 | mat$cas8a8 > 0] <- T
mat$type_i_a_specific[is.na(mat$type_i_a_specific)] <- F
mat$type_i_b_specific[mat$cas8b1 > 0 | mat$cas8b10 > 0 | mat$cas8b2 > 0 | mat$cas8b3 > 0 | mat$cas8b4 > 0 | mat$cas8b5 > 0 | mat$cas8b6 > 0  | mat$cas8b8 > 0 | mat$cas8b9 > 0] <- T
mat$type_i_b_specific[is.na(mat$type_i_b_specific)] <- F
mat$type_i_c_specific[mat$cas8c > 0] <- T
mat$type_i_c_specific[is.na(mat$type_i_c_specific)] <- F
mat$type_i_g_specific[mat$cas8u1 > 0 | mat$cas8u2 > 0] <- T
mat$type_i_g_specific[is.na(mat$type_i_g_specific)] <- F
mat$type_i_d_specific[mat$cas10d > 0 | mat$csc2gr7 > 0 | mat$csc1gr5 > 0] <- T
mat$type_i_d_specific[is.na(mat$type_i_d_specific)] <- F
mat$type_i_e_specific[mat$cas8e > 0 | mat$Cse1 > 0 | mat$cse2gr11 > 0 | mat$Cse2 > 0 | mat$cas6e > 0] <- T
mat$type_i_e_specific[is.na(mat$type_i_e_specific)] <- F
mat$type_i_f_specific[mat$cas3f > 0 | mat$cas5f > 0 | mat$cas6f > 0 | mat$cas7f > 0 | mat$cas8f > 0] <- T
mat$type_i_f_specific[is.na(mat$type_i_f_specific)] <- F
mat$type_iv_specific[mat$csf5gr6 > 0 | mat$csf1gr8 > 0 | mat$csf2gr7 > 0 | mat$csf3gr5 > 0 | mat$csf4gr11 > 0] <- T
mat$type_iv_specific[is.na(mat$type_iv_specific)] <- F
mat$type_iii_adf_specific[mat$csm2gr11 > 0 | mat$csm3gr7 > 0 | mat$csm4gr5 > 0 | mat$csm5gr7 > 0 | mat$csx10gr5 > 0 | mat$csx1 > 0 | mat$csx15 > 0 | mat$csx16 > 0 | mat$csx18 > 0 | mat$csx19 > 0 | mat$csx20 > 0 | mat$csx21 > 0 | mat$csx22 > 0 | mat$csx23 > 0 | mat$csx24 > 0 | mat$csx25 > 0 | mat$csx26 > 0 | mat$csx3 > 0] <- T
mat$type_iii_adf_specific[is.na(mat$type_iii_adf_specific)] <- F
mat$type_iii_bc_specific[mat$cmr1gr7 > 0 | mat$cmr3gr5 > 0 | mat$cmr4gr7 > 0 | mat$cmr5gr11 > 0 | mat$cmr6gr7 > 0 | mat$cmr7 > 0 | mat$cmr8gr7 > 0] <- T
mat$type_iii_bc_specific[is.na(mat$type_iii_bc_specific)] <- F

mat$type_ii_a_specific[mat$cas9 > 0 & mat$csn2 > 0] <- T
mat$type_ii_a_specific[is.na(mat$type_ii_a_specific)] <- F
mat$type_ii_b_specific[mat$cas9 > 0 & mat$cas4 > 0] <- T
mat$type_ii_b_specific[is.na(mat$type_ii_b_specific)] <- F
mat$type_ii_c_specific[mat$cas9 > 0 & mat$csn2 == 0 & mat$cas4 == 0] <- T
mat$type_ii_c_specific[is.na(mat$type_ii_c_specific)] <- F


mat$type_i_a[mat$adaptation_complex & mat$cas4 > 0 & mat$cas3_generic & mat$cas8_generic & mat$cas5_type_i & mat$cas7_type_i & mat$cas11_generic & mat$cas6_generic & mat$type_i_a_specific] <- "I-A"
mat$type_i_a[mat$type_i_a != "I-A"] <- ""
mat$type_i_b[mat$adaptation_complex & mat$cas4 > 0 & mat$cas3_generic & mat$cas8_generic & mat$cas5_type_i & mat$cas7_type_i & mat$cas6_generic & mat$type_i_b_specific] <- "I-B"
mat$type_i_b[mat$type_i_b != "I-B"] <- ""
mat$type_i_c[mat$adaptation_complex & mat$cas4 > 0 & mat$cas3_generic & mat$cas8_generic & mat$cas7_type_i & mat$type_i_c_specific] <- "I-C"
mat$type_i_c[mat$type_i_c != "I-C"] <- ""
mat$type_i_d[mat$adaptation_complex & mat$cas4 > 0 & mat$cas3_generic & mat$cas10d > 0 & mat$cas5_type_i & mat$cas7_type_i & mat$cas6_generic & mat$type_i_d_specific] <- "I-D"
mat$type_i_d[mat$type_i_d != "I-D"] <- ""
mat$type_i_e[mat$adaptation_complex & mat$cas3_generic & mat$type_i_e_specific] <- "I-E"
mat$type_i_e[mat$type_i_e != "I-E"] <- ""
mat$type_i_f[mat$adaptation_complex & mat$cas3_generic & mat$type_i_f_specific] <- "I-F"
mat$type_i_f[mat$type_i_f != "I-F"] <- ""
mat$type_i_g[mat$adaptation_complex & mat$cas4 > 0 & mat$cas3_generic & mat$cas8_generic & mat$cas5_type_i & mat$cas7_type_i & mat$cas6_generic & mat$type_i_g_specific] <- "I-G"
mat$type_i_g[mat$type_i_g != "I-G"] <- ""

mat$type_iii_a[mat$cas6_generic & mat$cas10 > 0 & mat$cas11_generic & mat$cas7_type_iii & mat$cas5_type_iii & mat$type_iii_adf_specific & mat$csx19 == 0] <- "III-A"
mat$type_iii_a[mat$type_iii_a != "III-A"] <- ""
mat$type_iii_b[mat$cas6_generic & mat$cas10 > 0 & mat$cas11_generic & mat$cas7_type_iii & mat$cas5_type_iii & mat$type_iii_bc_specific] <- "III-C"
mat$type_iii_b[mat$type_iii_b != "III-B"] <- ""
mat$type_iii_c[mat$cas10 > 0 & mat$cas11_generic & mat$cas7_type_iii & mat$cas5_type_iii & mat$type_iii_bc_specific] <- "III-C"
mat$type_iii_c[mat$type_iii_c != "III-C"] <- ""
mat$type_iii_d[mat$cas10 > 0 & mat$cas11_generic & mat$cas7_type_iii & mat$cas5_type_iii & mat$type_iii_adf_specific & mat$csx19 > 0] <- "III-D"
mat$type_iii_d[mat$type_iii_d != "III-D"] <- ""
mat$type_iii_f[mat$cas10 > 0 & mat$cas11_generic & mat$cas7_type_iii & mat$cas5_type_iii & mat$type_iii_adf_specific & mat$csx19 == 0] <- "III-F"
mat$type_iii_f[mat$type_iii_f != "III-F"] <- ""

mat$type_iv_a[mat$DinG > 0 & mat$cas6_generic & mat$csf1gr8 > 0 & mat$cas7_type_iv & mat$cas5_type_iv] <- "IV"
mat$type_iv_a[mat$type_iv_a != "IV"] <- ""
mat$type_iv_bc[mat$cas7_type_iv & mat$cas5_type_iv & mat$cas11_generic] <- "IV"
mat$type_iv_bc[mat$type_iv_bc != "IV"] <- ""

mat$type_ii_a[mat$cas9 > 0 & mat$adaptation_complex & mat$csn2 > 0] <- "II-A"
mat$type_ii_a[mat$type_ii_a != "II-A"] <- ""
mat$type_ii_b[mat$cas9 > 0 & mat$adaptation_complex & mat$cas4 > 0] <- "II-B"
mat$type_ii_b[mat$type_ii_b != "II-B"] <- ""
mat$type_ii_c[mat$cas9 > 0 & mat$adaptation_complex & mat$cas4 == 0 & mat$csn2 == 0] <- "II-C"
mat$type_ii_c[mat$type_ii_c != "II-C"] <- ""

mat$type_i[mat$cas3_generic & mat$cas7_type_i] <- "I"
mat$type_ii[mat$cas9 > 0] <- "II"
mat$type_iii[mat$cas10 > 0] <- "III"
mat$type_iv[mat$type_iv_specific] <- "IV"
mat$array_count[is.na(mat$array_count)] <- 0
mat2 <- mat
mat2[is.na(mat2)] <- ""

tmp <- mat2 %>% mutate(types = paste(type_i, type_ii, type_iii, type_iv, sep = ","),
                      subtypes = paste(type_i_a, type_i_b, type_i_c, type_i_d, type_i_e, type_i_f, type_i_g,
                                       type_ii_a, type_ii_b, type_ii_c, 
                                       type_iii_a, type_iii_b, type_iii_c, type_iii_d, type_iii_f,
                                       type_iv_a, type_iv_bc, sep = ", ")) %>% 
  select(types, subtypes, contig)


tmp <- data.frame(lapply(tmp, function(x) {
  gsub(", , ", ", ", x)
}))


ggplot() +
  geom_bar(data = tmp %>% filter(types != ","), aes(x = types))

colnames(mat2)

names <- list()

names$specific <- c("type_i_a_specific", "type_i_b_specific", "type_i_c_specific", "type_i_d_specific", "type_i_e_specific", "type_i_f_specific", "type_i_f_specific",
           "type_iii_adf_specific", "type_iii_bc_specific",
           "type_iv_specific")

names$subtypes <- c("type_i_a", "type_i_b", "type_i_c", "type_i_d", "type_i_e", "type_i_f", "type_i_g",
           "type_ii_a", "type_ii_b", "type_ii_c",
           "type_iii_a", "type_iii_b", "type_iii_c", "type_iii_d", "type_iii_f",
           "type_iv_a", "type_iv_bc")
names$types <- c("type_i", "type_ii", "type_iii", "type_iv")


for(value in names$specific){
  colNum <- match(value, table = colnames(mat2))
#  print(value)
#  print(colNum)
  
  tmpDat <- mat2[mat2[,colNum] ==T, colNum]
  print(paste(value, ": ", length(tmpDat)), sep = "")
}


for(value in names$subtypes){
  colNum <- match(value, table = colnames(mat2))
  #  print(value)
  #  print(colNum)
  
  tmpDat <- mat2[mat2[,colNum] != "", colNum]
  print(paste(value, ": ", length(tmpDat)), sep = "")
}

for(value in names$types){
  colNum <- match(value, table = colnames(mat2))
  #  print(value)
  #  print(colNum)
  
  tmpDat <- mat2[mat2[,colNum] != "", colNum]
  print(paste(value, ": ", length(tmpDat)), sep = "")
}



mat3 <- mat %>% filter(array_count > 0)

mat3[is.na(mat3)] <- ""

for(value in names$specific){
  colNum <- match(value, table = colnames(mat3))
  #  print(value)
  #  print(colNum)
  
  tmpDat <- mat3[mat3[,colNum] ==T, colNum]
  print(paste(value, ": ", length(tmpDat)), sep = "")
}


for(value in names$subtypes){
  colNum <- match(value, table = colnames(mat3))
  #  print(value)
  #  print(colNum)
  
  tmpDat <- mat3[mat3[,colNum] != "", colNum]
  print(paste(value, ": ", length(tmpDat)), sep = "")
}

for(value in names$types){
  colNum <- match(value, table = colnames(mat3))
  #  print(value)
  #  print(colNum)
  
  tmpDat <- mat3[mat3[,colNum] != "", colNum]
  print(paste(value, ": ", length(tmpDat)), sep = "")
}


upsetDatTypes <- mat2 %>% select(array_count, adaptation_complex, type_i, type_ii, type_iii, type_iv)
upsetDatSubtypes <- mat2 %>% select(array_count, type_i_a, type_i_b, type_i_c, type_i_d, type_i_e, type_i_f, type_i_g,
                                    type_ii_a, type_ii_b, type_ii_c, 
                                    type_iii_a, type_iii_b, type_iii_c, type_iii_d, type_iii_f,
                                    type_iv_a, type_iv_bc)
upsetDatSubtypeSignatures <- mat2 %>% select(array_count, adaptation_complex, cas01, type_i_a_specific, type_i_b_specific, type_i_c_specific, type_i_d_specific, type_i_e_specific, type_i_f_specific, type_i_g_specific,
                                             type_ii_a_specific, type_ii_b_specific, type_ii_c_specific, 
                                             type_iii_adf_specific, type_iii_bc_specific,
                                    type_iv_specific)

upsetDatSubtypes[upsetDatSubtypes == T] <- 1
upsetDatSubtypes[upsetDatSubtypes == F] <- 0
upsetDatSubtypes[upsetDatSubtypes == ""] <- 0
upsetDatSubtypes[upsetDatSubtypes > 0] <- 1
upsetDatSubtypes <- upsetDatSubtypes %>% mutate_all(as.numeric)

upsetDatSubtypeSignatures[upsetDatSubtypeSignatures == T] <- 1
upsetDatSubtypeSignatures[upsetDatSubtypeSignatures == F] <- 0
upsetDatSubtypeSignatures[upsetDatSubtypeSignatures == ""] <- 0
upsetDatSubtypeSignatures[upsetDatSubtypeSignatures > 0] <- 1
upsetDatSubtypeSignatures <- upsetDatSubtypeSignatures %>% mutate_all(as.numeric)


upsetDatTypes[upsetDatTypes == T] <- 1
upsetDatTypes[upsetDatTypes == F] <- 0
upsetDatTypes[upsetDatTypes == ""] <- 0
upsetDatTypes[upsetDatTypes > 0] <- 1
upsetDatTypes <- upsetDatTypes %>% mutate_all(as.numeric)

length(upsetDatTypes$type_i[upsetDatTypes$type_i ==1])
length(upsetDatSubtypeSignatures$type_i_e_specific[upsetDatSubtypeSignatures$type_i_e_specific ==1])

UpSetR::upset(upsetDatTypes, sets = colnames(upsetDatTypes), mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = T)
UpSetR::upset(upsetDatSubtypes, sets = colnames(upsetDatSubtypes), mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = T, nintersects = 100)
UpSetR::upset(upsetDatSubtypeSignatures, sets = colnames(upsetDatSubtypeSignatures), mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = T, nintersects = 100)

svg(filename="~/phd/RNASeq/figures/upsetr_crispr_types.svg",
    width=15,
    height=10,
    pointsize=12)
UpSetR::upset(upsetDatTypes, sets = colnames(upsetDatTypes), mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = T)

dev.off()

svg(filename="~/phd/RNASeq/figures/upsetr_crispr_subtypes.svg",
    width=15,
    height=10,
    pointsize=12)
UpSetR::upset(upsetDatSubtypes, sets = colnames(upsetDatSubtypes), mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = T, nintersects = 100)

dev.off()

svg(filename="~/phd/RNASeq/figures/upsetr_crispr_signatures.svg",
    width=15,
    height=10,
    pointsize=12)
UpSetR::upset(upsetDatSubtypeSignatures, sets = colnames(upsetDatSubtypeSignatures), mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = T, nintersects = 100)

dev.off()


arrayLengthsTypes <- mat2 %>% select(array_count, spacer_count, mean_spacers_per_array, 
                                     adaptation_complex, type_i, type_ii, type_iii, type_iv) 

arrayLengthsTypes[arrayLengthsTypes == T] <- 1
arrayLengthsTypes[arrayLengthsTypes == F] <- 0
arrayLengthsTypes[arrayLengthsTypes == ""] <- 0


type_i_dat <- arrayLengthsTypes %>% filter(type_i == "I", type_ii != "II", type_iii != "III") %>% 
  select(array_count, spacer_count, mean_spacers_per_array, adaptation_complex) %>% mutate_all(as.numeric) %>% 
  mutate(group = "I")

type_ii_dat <- arrayLengthsTypes %>% filter(type_ii == "II", type_iii != "III", type_i != "I") %>% 
  select(array_count, spacer_count, mean_spacers_per_array, adaptation_complex) %>% mutate_all(as.numeric) %>% 
  mutate(group = "II")

type_iii_dat <- arrayLengthsTypes %>% filter(type_iii == "III", type_i != "I", type_ii != "ii") %>% 
  select(array_count, spacer_count, mean_spacers_per_array, adaptation_complex) %>% mutate_all(as.numeric) %>% 
  mutate(group = "III")

type_iv_dat <- arrayLengthsTypes %>% filter(type_iv == "IV") %>% 
  select(array_count, spacer_count, mean_spacers_per_array, adaptation_complex) %>% mutate_all(as.numeric) %>% 
  mutate(group = "IV")


no_type_dat <- arrayLengthsTypes %>% filter(type_i != "I", type_ii != "II", type_iii != "III", type_iv != "IV") %>% 
  select(array_count, spacer_count, mean_spacers_per_array, adaptation_complex) %>% mutate_all(as.numeric) %>% 
  mutate(group = "No system")

arraysDat <- type_i_dat %>% bind_rows(type_ii_dat, type_iii_dat, type_iv_dat, no_type_dat) %>% filter(array_count > 0)

types.p <- ggplot() +
  geom_boxplot(data = arraysDat %>% mutate(group = ifelse(group == "IV", ifelse(adaptation_complex == 0, "IV", "NA"),group)) %>% filter(group != "NA") , aes(x = group, y = spacer_count, fill = group), outlier.shape = NA, show.legend = F)
types.p

adapt.p <- ggplot() +
  geom_boxplot(data = arraysDat, aes(x = as.character(adaptation_complex), y = spacer_count, fill = as.character(adaptation_complex)), outlier.shape = NA, show.legend = F)
adapt.p



all.p <- ggarrange(types.p, adapt.p + rremove("x.text"),  
                   labels = LETTERS[1:2],
                   ncol = 2, nrow = 1)
all.p

ggsave(filename = "~/phd/RNASeq/figures/spacers_boxplot.svg", plot = all.p, width = 450, height = 300, units = "mm")


cumulativeCounts <- function(dists, smooth = T){
  
  groups <- unique(dists$group)
  for(i in groups){
    dat <- dists %>% filter(group == i)
    dat <- dat %>% mutate(count = 1) %>% 
      arrange(-max_dist) %>% group_by(group) %>% 
      mutate(cumulativeCount = cumsum(count)) %>% ungroup() %>% 
      group_by(group, max_dist) %>% summarise(cumulative_prop = max(cumulativeCount)/ nrow(dat))
    
    if(smooth){
      dat <- as.data.frame(spline(x = dat$max_dist,y =  dat$cumulative_prop))
    }
    dat <- dat %>% ungroup() %>% mutate(group = i)
    if(exists('combinedDat')){
      combinedDat <- combinedDat %>% bind_rows(dat)
    }else{
      combinedDat <- dat 
    }
  }
  return(combinedDat)  
  
}

meanSpacersCounts <-cumulativeCounts(dists = arraysDat %>% filter(group != "No system", group != "IV") %>% dplyr::rename(max_dist = mean_spacers_per_array))
totslSpacersCounts <-cumulativeCounts(dists = arraysDat %>% filter(group != "No system", group != "IV") %>% dplyr::rename(max_dist = spacer_count))




p <- ggplot() +
  geom_line(data = meanSpacersCounts, aes(x= x, y = y, group = group, colour = group)) +
  scale_x_continuous(trans = "log10") +
  theme_classic() 
p

p <- ggplot() +
  geom_line(data = totslSpacersCounts, aes(x= x, y = y, group = group, colour = group)) + 
  scale_x_continuous(trans = "log10") +
  theme_classic() 
p
ggsave(filename = "~/phd/RNASeq/figures/crispr_array_counts.svg", plot = p, width = 450, height = 307, units = "mm")

ggplot() +
  geom_freqpoly(data = arraysDat, aes(x = spacer_count, y = ..density.., group = group, color = group), binwidth = 5)







ggplot() +
  geom_freqpoly(data = arraysDat, aes(x = mean_spacers_per_array, y = ..density.., group = group, color = group), binwidth = 5)


p <-ggplot() +
  geom_freqpoly(data = arraysDat %>% filter(group != "No system", group != "IV"), aes(x = spacer_count, y = ..count.., group = group, color = group), binwidth = 10) +
  scale_y_continuous(trans = "log10")
p
ggsave(filename = "~/phd/RNASeq/figures/crispr_array_counts.svg", plot = p, width = 450, height = 307, units = "mm")


ggplot() +
  geom_freqpoly(data = arraysDat %>% filter(group != "No system", group != "IV"), aes(x = mean_spacers_per_array, y = ..density.., group = group, color = group), binwidth = 10) +
  scale_y_continuous(trans = "log10")



ggplot() +
  geom_freqpoly(data = arraysDat, aes(x = spacer_count, y = ..density.., group = as.character(adaptation_complex), color = as.character(adaptation_complex)), binwidth = 5)


ggplot() +
  geom_freqpoly(data = arraysDat, aes(x = mean_spacers_per_array, y = ..density.., group = as.character(adaptation_complex), color = as.character(adaptation_complex)), binwidth = 5)








blastAphage <- read.table("~/phd/priming/phageq.archaea.blastn")
blastBphage <- read.table("~/phd/priming/phageq.bacteria.blastn", fill = T)

colnames(blastAphage) <- c("query.acc.", "subject.acc.", "percent.identity", "alignment.length", 
                           "mismatches", "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score")
colnames(blastBphage) <- c("query.acc.", "subject.acc.", "percent.identity", "alignment.length", 
                           "mismatches", "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score")

swipeAphage <- read.table("~/phd/priming/phageq2.archaea.swipe._-1_1")
swipeBphage <- read.table("~/phd/priming/phageq2.bacteria.swipe._-1_1", fill = T)
swipeALookup <- read.table("~/phd/priming/phageq2.archaea.swipe._-1_1.lookup", fill = T)
swipeBLookup <- read.table("~/phd/priming/phageq2.bacteria.swipe._-1_.lookup", fill = T)

colnames(swipeAphage) <- c("query.acc.", "subject.name", "percent.identity", 
                           "alignment.length", "mismatches", "gap.opens", 
                           "q.start", "q.end", "s.start", "s.end", "bit.score")
colnames(swipeBphage) <- c("query.acc.", "subject.name", "percent.identity", 
                           "alignment.length", "mismatches", "gap.opens", 
                           "q.start", "q.end", "s.start", "s.end", "bit.score")

colnames(swipeALookup) <- c("subject.name", "subject.acc.")
colnames(swipeBLookup) <- c("subject.name", "subject.acc.")

swipeALookup <- swipeALookup %>% filter(subject.acc. != "")
swipeBLookup <- swipeBLookup %>% filter(subject.acc. != "")

swipeALookup <- data.frame(lapply(swipeALookup, function(x) {
  gsub(">", "", x)
}))

swipeBLookup <- data.frame(lapply(swipeBLookup, function(x) {
  gsub(">", "", x)
}))

swipeAphage <- swipeAphage %>% full_join(swipeALookup, by = "subject.name")
swipeBphage <- swipeBphage %>% full_join(swipeBLookup, by = "subject.name")


swipePhage <- swipeAphage %>% bind_rows(swipeBphage)

accessions <- swipePhage %>% select(subject.acc.)

blastPhage <- blastAphage %>% bind_rows(blastBphage) %>% full_join(accessions)

swipePhage[is.na(swipePhage)] <- 0
blastPhage[is.na(blastPhage)] <- 0

swipeCounts <- swipePhage %>% group_by(subject.acc.) %>% summarise(count = n())
blastCounts <- blastPhage %>% group_by(subject.acc.) %>% summarise(count = n())
swipeScores <- swipePhage %>% group_by(subject.acc.) %>% summarise(score = max(bit.score))
blastScores  <- blastPhage %>% group_by(subject.acc.) %>% summarise(score = max(bit.score))
swipeLengths <- swipePhage %>% group_by(subject.acc.) %>% summarise(alignment.length = median(alignment.length))
blastLengths  <- blastPhage %>% group_by(subject.acc.) %>% summarise(alignment.length = median(alignment.length))
swipeMismatches <- swipePhage %>% group_by(subject.acc.) %>% summarise(mismatches = median(mismatches))
blastMismatches  <- blastPhage %>% group_by(subject.acc.) %>% summarise(mismatches = median(mismatches))

swipeData <- swipeCounts %>% 
  full_join(swipeScores, by = "subject.acc.") %>% 
  full_join(swipeLengths, by = "subject.acc.") %>% 
  full_join(swipeMismatches, by = "subject.acc.") %>% 
  mutate(group = "swipe")

blastData <- blastCounts %>% 
  full_join(blastScores, by = "subject.acc.") %>% 
  full_join(blastLengths, by = "subject.acc.") %>% 
  full_join(blastMismatches, by = "subject.acc.") %>% 
  mutate(group = "blast")

dat <- swipeData %>% bind_rows(blastData)

mergeDat <- swipeData %>% full_join(blastData, by = "subject.acc.")


align.p <- ggplot() +
  geom_freqpoly(data = dat, aes(x = alignment.length, group = group, color = group), binwidth = 1, show.legend = F) +
  xlim(min = 0, max = NA)

mismatch.p <- ggplot() +
  geom_freqpoly(data = dat, aes(x = mismatches, group = group, color = group), binwidth = 1, show.legend = F) +
  xlim(min = 0, max = NA)

score.p <- ggplot() +
  geom_freqpoly(data = dat, aes(x = score, group = group, color = group), binwidth = 3, show.legend = F) +
  xlim(min = 0, max = NA)



all.p <- ggarrange(align.p, mismatch.p, score.p + rremove("x.text"),  
                   labels = LETTERS[1:3],
                   ncol = 1, nrow = 3)
all.p

ggsave(filename = "~/phd/RNASeq/figures/swipe_vs_blast.svg", plot = all.p, width = 300, height = 450, units = "mm")



length(mergeDat$score.x[is.na(mergeDat$score.x)])
length(mergeDat$score.x[!is.na(mergeDat$score.x)])

length(mergeDat$score.y[is.na(mergeDat$score.y)])
length(mergeDat$score.y[!is.na(mergeDat$score.y)])


cas_genes <- read.table("/Volumes/2TTJN/Priming_Required_files/cas_genes_refseq_83.txt", sep = "\t", comment.char = "", quote = "", 
                        header = T)

tmp <- cas_genes %>% select(Gene_name, CRISPR.Cas_system_associated, x) %>% unique()



