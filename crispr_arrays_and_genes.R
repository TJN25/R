library(tidyverse)
library(reshape2)
library(UpSetR)
library(ggpubr)
# library(multidplyr)
library(parallel)
library(lme4)
library(readr)

# dat <- read.table("/Volumes/2TTJN/Priming_Required_files/cas_arrays_genomes_refseq_83.txt", 
                  # sep = "\t", 
                  # quote = "", 
                  # comment.char = "", 
                  # fill = T, header = T)



 # cas_arrays_genomes_refseq_83 <- dat
 # save(cas_arrays_genomes_refseq_83,file = "~/bin/large_files/cas_arrays_genomes_refseq_83.Rda")
load("~/bin/large_files/cas_arrays_genomes_refseq_83.Rda") 
dat <- cas_arrays_genomes_refseq_83

final_data <- read_csv("~/phd/CRISPR/R-script-and-files/final_hits_with_domain_wiki_added.csv")
head(final_data)

table(final_data$domain_2, final_data$Subtype)


dat <- dat %>% filter(domainName != "", !is.na(subtype.present))

dat <- dat %>% select(assembly_accession, refseq_category, taxid, species_taxid, organism_name, assembly_level,
                      domainName, protein.count, cas.proteins, cas1, cas2, basic.genes, subtype.list, subtype.count, single.system, contig.count, proteins.present, arrays.present, subtype.present, array.and_subtype)

head(cas_arrays_genomes_refseq_83)

genomes_taxa <- cas_arrays_genomes_refseq_83 %>% 
  select(domainName, gbrs_paired_asm, organism_name, assembly_accession, taxid) %>% 
  unique()

final_genomes <- final_data %>% select(host.target.pair, hits.count, host.acc., phages.targeted, 
                                       hosts.targeting, Subtype) %>% unique() %>% 
  separate(col = host.acc., into = c("assembly_accession", NA), sep = "\\$", remove = F)

final_genomes <- final_genomes %>% left_join(genomes_taxa) %>% arrange(organism_name) %>% 
  separate(col = organism_name, into = c("genus", NA), sep = " ", remove = F)

final_genomes <- data.frame(mclapply(final_genomes, function(x) {
  gsub("\\[", "", x)
}))
final_genomes <- data.frame(mclapply(final_genomes, function(x) {
  gsub("\\]", "", x)
}))

write.csv(final_genomes, file = "~/phd/CRISPR/R-script-and-files/final_hits_with_domain.csv", quote = T,row.names = F)



genes <- read.table("~/phd/CRISPR/Priming_Required_files/cas_genes_refseq_83.txt", 
                    sep = "\t", 
                    quote = "", 
                    comment.char = "", 
                    fill = T, header = T)


# arrays <- read.table("/Volumes/2TTJN/Priming_Required_files/all_CRISPRDetect_refseq_83.TEST.fna", 
#                     sep = " ", 
#                     quote = "", 
#                     comment.char = "", 
#                     fill = T, header = F)
# save(arrays,file = "~/bin/large_files/arrays.Rda")
load("~/bin/large_files/arrays.Rda")

genes <- genes %>% separate(col = filepath, into = c(NA, "t1"), sep = "#", remove = F, extra = "merge")
genes <- genes %>% separate(col = t1, into = c("assembly_accession", NA), sep = "\\/", remove = F, extra = "drop")



arrays_high <- arrays %>% filter(grepl(pattern = "LOW", x = V1) == F)

arrays_high <- arrays_high %>% separate(col = V1, into = c("contig", "t2"), sep = "-", extra = "merge", remove = F)

arrays_high <- arrays_high %>% separate(col = t2, into = c("genome_name", "array_coords", "t3"), sep = "\\|", extra = "merge", remove = T)

arrays_high <- arrays_high %>% separate(col = t3, into = c("array_number", "spacer_number", "t4"), sep = "_", extra = "merge", remove = T)

arrays_high <- arrays_high %>% separate(col = t4, into = c("repeat_seq", "array_data", "spacer_seq"), sep = "\\[", extra = "merge", remove = T)

arrays_high <- data.frame(mclapply(arrays_high, function(x) {
  gsub(">", "_", x)
}))


save(arrays_high,file = "~/bin/large_files/arrays_high.Rda")






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


swipe_counts2 <- swipePhage  %>% group_by(query.acc.) %>% summarise(counts = n())
archaea_crispr <- swipeAphage %>% filter(bit.score >= 23) 
archaea_crispr <- swipePhage %>% separate(col = subject.acc., into = c(NA, "t1", NA), sep = "\\[", remove = F)
archaea_crispr <- archaea_crispr %>% separate(col = t1, into = c(NA, NA, NA, "taxonomy"), sep = "_", remove = F)
archaea_crispr <- archaea_crispr %>% select(query.acc., taxonomy, bit.score) %>% separate(col = taxonomy, into = c("do", "ki", "ph", "cl", "or", "fa", "ge", "sp"), sep = ";", remove = F)

archaea_crispr <- archaea_crispr %>% filter(!is.na(ph))



targets_kept <- archaea_crispr %>% filter(bit.score >= 23) %>% group_by(query.acc.) %>% summarise(min_bit_score = min(bit.score))

archaea_summary <- archaea_crispr %>% group_by(query.acc.) %>% summarise(domains = paste(unique(do), collapse = ","),
                                                                  phylums = paste(unique(ph), collapse = ","),
                                                                  classes = paste(unique(cl), collapse = ","),
                                                                  orders = paste(unique(or), collapse = ","),
                                                                  families = paste(unique(fa), collapse = ","))

archaea_counts <- archaea_crispr %>% group_by(query.acc.) %>% summarise(domains = length(unique(do)),
                                                                 phylums = length(unique(ph)),
                                                                 classes = length(unique(cl)),
                                                                 orders = length(unique(or)),
                                                                 families = length(unique(fa))
)

archaea_freq <- archaea_crispr %>% group_by(query.acc.) %>% summarise(counts = n())

archaea <- archaea_summary %>% full_join(archaea_counts, by = "query.acc.") %>% 
  full_join(archaea_freq, by = "query.acc.") %>% 
  full_join(targets_kept, by = "query.acc.")

archaea_table_high <- as.data.frame(table(archaea$domains.x[archaea$min_bit_score >=23], archaea$domains.y[archaea$min_bit_score >=23]))
archaea_table_high <- as.data.frame(acast(data = archaea_table_high, formula = Var1 ~ Var2))


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
mergeDat <- mergeDat %>% 
  mutate(align.x.bool = ifelse(alignment.length.x > 0, 1, 0),
         align.y.bool = ifelse(alignment.length.y > 0, 1, 0)) %>% 
  mutate(align.sum = align.x.bool + align.y.bool) %>% 
  mutate(align.diff = alignment.length.x - alignment.length.y)


table(mergeDat$align.sum)



align.p <- ggplot() +
  geom_freqpoly(data = dat, aes(x = alignment.length, group = group, color = group), binwidth = 1, show.legend = F) +
  xlim(min = 0, max = NA)

mismatch.p <- ggplot() +
  geom_freqpoly(data = dat, aes(x = mismatches, group = group, color = group), binwidth = 1, show.legend = F) +
  xlim(min = 0, max = NA)

score.p <- ggplot() +
  geom_freqpoly(data = dat, aes(x = score, group = group, color = group), binwidth = 3, show.legend = F) +
  xlim(min = 0, max = NA)
len.p <- ggplot() +
  geom_histogram(data = mergeDat %>% filter(align.sum == 2), aes(x = align.diff), binwidth = 1)



all.p <- ggarrange(align.p, mismatch.p, score.p, len.p + rremove("x.text"),  
                   labels = LETTERS[1:4],
                   ncol = 1, nrow = 4)
all.p

ggsave(filename = "~/phd/thesis_figures/SVG/swipe_vs_blast_2.svg", plot = all.p, width = 300, height = 600, units = "mm")


mergeDat$subject.acc.[1]
spacers <- read.table("~/phd/CRISPR/Priming_Required_files/spacer_information_for_spacer_id_KEEP_THIS.txt", 
                    sep = "\t", 
                    quote = "", 
                    comment.char = "", 
                    fill = T, header = T)

spacers <- spacers %>% filter(spacer_id != ">NA")

spacers <- data.frame(lapply(spacers, function(x) {
  gsub(">", "", x)
}))

spacers <- spacers %>% mutate(spacer.length = nchar(as.character(V2)))

spacers <- spacers %>% filter(spacer.length < 60, spacer.length > 15)

spacer.small <- spacers %>% select(V1, spacer.length) %>% 
  dplyr::rename(subject.acc. = V1)
spacer.small <- spacer.small %>% separate(subject.acc., into = c("Contig", "position", "other", "other.2"), sep = "\\|", remove = T)
spacer.small <- spacer.small %>% separate(other, into = c("array.num", "spacer.num", NA), sep = "_",remove = F)
spacer.small <- spacer.small %>% separate(Contig, into = c("contig.small", NA), sep = "-",remove = F)
spacer.small <- spacer.small %>% mutate(spacer.id.2 = paste(contig.small, array.num, spacer.num, sep = "-"))

final_data <- read.csv("~/phd/CRISPR/Priming_Required_files/FinalData.csv",as.is=T)


mergeDat <- mergeDat %>% separate(subject.acc., into = c("Domain", "Subtype", "Contig", "position", "other"), sep = "\\|", remove = F)
mergeDat$other[1]
mergeDat <- mergeDat %>% separate(other, into = c("other.2", NA), sep = "\\[",remove = F)
mergeDat <- mergeDat %>% separate(other.2, into = c("array.num", "spacer.num"), sep = "_",remove = F)
mergeDat <- mergeDat %>% separate(Contig, into = c("contig.small", NA), sep = "-",remove = F)
mergeDat <- mergeDat %>% mutate(spacer.id.2 = paste(contig.small, array.num, spacer.num, sep = "-"))


spacer.small.2 <- spacer.small %>% select(spacer.id.2, spacer.length)
spacer.small$subject.acc.[1]


mergeDat.2 <- mergeDat %>% full_join(spacer.small.2, by = "spacer.id.2")

tmp$Contig[1]


tmp <- mergeDat %>% filter(align.sum == 2, !is.na(spacer.length))

tmp <- tmp %>% mutate(alignment.length.diff = spacer.length - alignment.length.x)



final_data <- final_data %>% separate(array.id, into = c("array.id.2", NA), sep = "$", extra = 'drop', remove = F)

final_data <- final_data %>% mutate(spacer.id.2 = paste(array.id.2, spacer.number, sep = "-"))


tmp <- tmp %>% full_join(final_data, by = "spacer.id.2")

tmp %>% filter(alignment.length.diff >= 0, alignment.length.diff < 3) %>%  nrow()


tmp2 <- tmp %>% filter(!is.na(PPS.score), !is.na(Subtype.x))


ggplot() +
  geom_histogram(data = tmp2, aes(x = alignment.length.diff), binwidth = 1)

ggplot() +
  geom_point(data = tmp2 %>% filter(hits.count > 3, spacer.order.number == 2), aes(x = score.x, y =alignment.length.diff))
ggplot() +
  geom_point(data = tmp2 %>% filter(hits.count > 3, spacer.order.number == 2), aes(x = score.x, y =alignment.length.diff))
ggplot() +
  geom_histogram(data = final_data %>% filter(hits.count > 3, spacer.order.number == 2), aes(x = bit.score.num), binwidth = 1)



img_crispr_input <- read.table("~/phd/CRISPR/Priming_Required_files/refseq_83_swipe_setup_IMG_ALL.txt", 
                         sep = "\t", header = T, comment.char = "", quote = "")

img_crispr <- img_crispr_input %>% select(target.acc., spacer.info, bit.score) %>% unique()

img_crispr <- img_crispr %>% separate(col = spacer.info, into = c(NA, "t1", NA), sep = "\\[", remove = F)
img_crispr <- img_crispr %>% separate(col = t1, into = c(NA, NA, NA, "taxonomy"), sep = "_", remove = F)
img_crispr <- img_crispr %>% select(target.acc., taxonomy, bit.score) %>% separate(col = taxonomy, into = c("do", "ki", "ph", "cl", "or", "fa", "ge", "sp"), sep = ";", remove = F)

img_crispr <- img_crispr %>% filter(!is.na(ph))

img_crispr[img_crispr == "20001Bacteria"] <- "Bacteria"


targets_kept <- img_crispr %>% filter(bit.score >= 23) %>% group_by(target.acc.) %>% summarise(min_bit_score = min(bit.score))

img_summary <- img_crispr %>% group_by(target.acc.) %>% summarise(domains = paste(unique(do), collapse = ","),
                                                                  phylums = paste(unique(ph), collapse = ","),
                                                                  classes = paste(unique(cl), collapse = ","),
                                                                  orders = paste(unique(or), collapse = ","),
                                                                  families = paste(unique(fa), collapse = ","))

img_counts <- img_crispr %>% group_by(target.acc.) %>% summarise(domains = length(unique(do)),
                                                                  phylums = length(unique(ph)),
                                                                  classes = length(unique(cl)),
                                                                  orders = length(unique(or)),
                                                                  families = length(unique(fa))
                                                                 )

img_freq <- img_crispr %>% group_by(target.acc.) %>% summarise(counts = n())

img <- img_summary %>% full_join(img_counts, by = "target.acc.") %>% 
  full_join(img_freq, by = "target.acc.") %>% 
  full_join(targets_kept, by = "target.acc.")

img$min_bit_score[is.na(img$min_bit_score)] <- 0

ggplot() +
  geom_histogram(data = img_counts, aes(x = families), binwidth = 1)


img_table <- as.data.frame(table(img$domains.x, img$families.y))
img_table <- as.data.frame(acast(data = img_table, formula = Var1 ~ Var2))

img_table_high <- as.data.frame(table(img$domains.x[img$min_bit_score >=23], img$families.y[img$min_bit_score >=23]))
img_table_high <- as.data.frame(acast(data = img_table_high, formula = Var1 ~ Var2))

cas_arrays_genomes_refseq_83 <- cas_arrays_genomes_refseq_83 %>% mutate(proteins_and_arrays = as.logical(proteins.present) + as.logical(arrays.present))


cassystems <- cas_arrays_genomes_refseq_83 %>% select(domainName, subtype.list) %>% 
  filter(!is.na(subtype.list),
         subtype.list != "",
         subtype.list != "FALSE")


cassystems <- data.frame(lapply(cassystems, function(x) {
  gsub("-A", "", x)
}))
cassystems <- data.frame(lapply(cassystems, function(x) {
  gsub("-B", "", x)
}))
cassystems <- data.frame(lapply(cassystems, function(x) {
  gsub("-C", "", x)
}))
cassystems <- data.frame(lapply(cassystems, function(x) {
  gsub("-D", "", x)
}))
cassystems <- data.frame(lapply(cassystems, function(x) {
  gsub("-E", "", x)
}))
cassystems <- data.frame(lapply(cassystems, function(x) {
  gsub("-F", "", x)
}))

cassystems <- cassystems %>% mutate(types = paste(strsplit(as.character(subtype.list), split = "_")[[1]], collapse = ", "))

for(i in 1:nrow(cassystems)){
  cassystems$types[i] <- paste(sort(unique(strsplit(as.character(cassystems$subtype.list[i]), split = "_")[[1]])), collapse = ",")
}
table(cas_arrays_genomes_refseq_83$domainName, cas_arrays_genomes_refseq_83$proteins_and_arrays)
typesDat <- as.data.frame(table(cassystems$domainName, cassystems$types))

typesDat <- as.data.frame(acast(data = typesDat, formula = Var1 ~ Var2))

typesDat$total <- typesDat %>% rowSums()

length(mergeDat$score.x[is.na(mergeDat$score.x)])
length(mergeDat$score.x[!is.na(mergeDat$score.x)])

length(mergeDat$score.y[is.na(mergeDat$score.y)])
length(mergeDat$score.y[!is.na(mergeDat$score.y)])


cas_genes <- read.table("/Volumes/2TTJN/Priming_Required_files/cas_genes_refseq_83.txt", sep = "\t", comment.char = "", quote = "", 
                        header = T)

tmp <- cas_genes %>% select(Gene_name, CRISPR.Cas_system_associated, x) %>% unique()




quad_dat <- read_csv("~/phd/CRISPR/R-script-and-files/Output/Graph_QuadDataAllClust.csv")

quad_dat <- quad_dat %>% mutate(group = ifelse(Var1 == "t_3" | Var1 == "n_3", "3'", 
                                               ifelse(Var1 == "t_5" | Var1 == "n_5", "5'", NA)))

direction_dat <- quad_dat %>% group_by(group) %>% 
  summarise_if(is.numeric, sum, na.rm = TRUE) %>% 
  filter(!is.na(group)) %>% 
  gather("subtype", "hits.count", -group) %>% 
  filter(grepl(pattern = "freq", x = subtype))


ggplot() +
  geom_bar(data = direction_dat, aes(x = subtype, y = hits.count, group = group, fill = group), stat = "identity", position = "dodge")

subtype_val <- "freq_II-C"
binom.test(direction_dat$hits.count[direction_dat$subtype == subtype_val],p=0.5)


binom.test(c(43, 56),p=0.5)



p-value = 0.002857
p-value = 0.3426
p-value = 0.9362
p-value = 7.572e-09
p-value = 0.651
p-value = 0.04191