#!/usr/bin/env Rscript

##no idea which packages are needed
library(tidyverse)
library(devtools)
library(lubridate)
library(genoPlotR)
library(drake)
library(ape)
library(Biostrings)
library(ROSE)
library(reshape2)

suppressMessages(library(tidyverse))
suppressMessages(library(comparativeSRA))
genera_list <- list.files("~/phd/RNASeq/genera/")

genera <- "Acinetobacter"

for(genera in genera_list){
  print(genera)
  if(genera == "Escherichia"){
    next
  }
  
  
  accession_list <- list.files(paste("~/phd/RNASeq/genera/",  genera, "/", sep = ""), pattern = ".data$")
  accessionsDat <- data.frame(accession_list = accession_list)
  
  accessionsDat <- accessionsDat %>% mutate_all(as.character)
  accession_folder <-  accessionsDat$accession_list[1]
  
  
  for(accession_folder in accessionsDat$accession_list){
    accession <- unlist(strsplit(accession_folder, "\\."))[c(1,2)]
    accession <- paste(accession, collapse = ".")
    files <- list.files(paste("~/phd/RNASeq/genera/",  genera, "/", accession ,".data/plot_files/", sep = ""), pattern = ".plot$")
    
    filesDat <- data.frame(files = files)
    
    filesDat <- filesDat %>% filter(grepl(pattern = "ncRNA", x = files)==F,
                                    grepl(pattern = "fwd", x = files)==F,
                                    grepl(pattern = "rev", x = files)==F,
                                    grepl(pattern = accession, x = files)==F) %>% mutate_all(as.character)
    
    
    if(nrow(filesDat) == 0){
      next
    }
    print(accession)
    
    gffDat <- read.table(paste("~/phd/RNASeq/genera/",  genera, "/", accession ,".data/gff_files/", accession, ".gff", sep = ""), sep = "\t", fill = T, comment.char = "#", quote = "")
    colnames(gffDat) <- c("sequence", "source", "feature", "start", "end", "score", "strand", "phase", "Atrribute")
    
    regions <- gffDat %>% filter(feature == "region", start == "1") %>% select(sequence, start, end)
    
    current_pos <- 0
    for(i in 1:nrow(regions)){
      # print(regions$start[i])
      regions$end[i] <- regions$end[i] + current_pos
      regions$start[i] <- regions$start[i] + current_pos
      current_pos <- regions$end[i]
      
    }
    
    file_name <-  filesDat$files[1]
    
    for(file_name in filesDat$files){
      print(file_name)
      sra_id <- unlist(strsplit(file_name, "\\."))[1]
      
      plotFile <- read.table(paste("~/phd/RNASeq/genera/",  genera, "/", accession ,".data/plot_files/", file_name, sep = ""))
      
      
      
      
      
      
      cdsDat <- gffDat %>% filter(feature == "CDS")
      
      for(i in 1:nrow(regions)){
        contig <- regions$sequence[i]
        increase_val <- regions$start[i]
        cdsDat <- cdsDat %>% mutate(start = ifelse(sequence == contig, start + increase_val, start),
                                    end = ifelse(sequence == contig, end + increase_val, end))
        
      }
      
      
      
      positions <- cdsDat %>% select(sequence,start, end) %>% mutate(start = start - 149) %>% mutate(start = ifelse(start < 1, 1, start))  %>% mutate(end = end + 149) %>% mutate(end = ifelse(end > max(regions$end), max(regions$end), end))
      
      
      
      
      
      
      positions$new <- do.call(paste, c(positions, sep=":")) 
      positions$new <- Map(":", positions$start, positions$end)
      
      
      values <- unlist(positions$new[1:1000])
      
      
      positions <- data.frame(rowNum = values, type="CDS")
      
      starts <- data.frame(rowNum = cdsDat$start, type2 ="start")
      stops <- data.frame(rowNum = cdsDat$end, type2 ="stop")
      starts <- starts %>% mutate(rowNum = rowNum - 149) %>% mutate(rowNum = ifelse(rowNum < 1, 1, rowNum))
      stops <- stops %>% mutate(rowNum = rowNum + 149) %>% mutate(rowNum = ifelse(rowNum > max(regions$end), max(regions$end), rowNum))
      
      starts_and_stops <- starts %>% bind_rows(stops)
      
      positions <- positions %>% 
        left_join(starts_and_stops, by = "rowNum") %>% 
        mutate(type2 = ifelse(is.na(type2), "remove", type2))
      
      plotDat <- plotFile %>% mutate(value = ifelse(V1 > V2, V1, V2)) %>% select(value) %>% mutate(rowNum = row_number())
      
      plotDat <- plotDat %>% left_join(positions, by = "rowNum") %>% mutate(type2 = ifelse(is.na(type2), "intergenic", type2)) %>% 
        filter(type2 != "remove")
      
      
      plotDat <- plotDat %>% unique()
      
      
      
      plotDat <- plotDat %>% mutate(type = ifelse(is.na(type), "intergenic", "CDS"))
      
      total <- sum(plotDat$value)/1000000
      
      plotDat$value <- plotDat$value/total
      
      # save(plotDat, file = "~/phd/RNASeq/tmp/plotDat.Rda")
      # 
      # load("~/phd/RNASeq/tmp/plotDat.Rda")
      
      checkDat <- plotDat %>% filter(value > 15/total & type == "intergenic" )
      
      
      # checkDat <- data.frame(value = "0.04", rowNum = c(1,2,3,4,5,15,16,17,18,21,22,23,50,51,52,56,61,62,63,64,100,101,102,103),
      #                        type = "intergenic")
      # checkDat <- checkDat %>% mutate(type = ifelse(rowNum > 15 & rowNum < 30, "CDS", "intergenic"))
      
      plotDat <- checkDat
      
      callsDat <- data.frame(start = rep(as.character("0"), 20000), stop = rep(as.character("0"), 20000), stringsAsFactors = F)
      start <- 0
      stop <- 0
      current_feature <- F
      cds <- F
      id_pos <- 0
      
      
      i <- 1
      for(i in 1:nrow(plotDat)){
        printRemaining(i = i, length = nrow(plotDat))
        if(i == nrow(plotDat)){
          if(current_feature){
            if((plotDat$rowNum[i] - 25) < stop){
              stop <- plotDat$rowNum[i]
              if(plotDat$type[i] == "CDS"){
                cds <- T
              }
            }
            if(cds){
              current_feature <- F
              cds <- F
              next
            }
            id_pos <- id_pos + 1
            # print(paste(start,stop, id_pos))
            callsDat$start[id_pos] <- start
            callsDat$stop[id_pos] <- stop
            current_feature <- F
            cds <- F
          }
        }
        
        if(current_feature == F){
          start <- plotDat$rowNum[i]
          stop <- plotDat$rowNum[i]
          if(plotDat$type[i] == "CDS"){
            cds <- T
          }
          current_feature <- T
          next
        }
        
        if(plotDat$type2[i] == "stop"){
          stop <- plotDat$rowNum[i]
          cds <- T
          next
        }
        
        if((plotDat$rowNum[i] - 25) < stop){
          stop <- plotDat$rowNum[i]
          if(plotDat$type[i] == "CDS"){
            cds <- T
          }
          next
        }else{
          if(cds){
            current_feature <- F
            cds <- F
            start <- plotDat$rowNum[i]
            stop <- plotDat$rowNum[i]
            next
          }
          id_pos <- id_pos + 1
          # print(paste(start,stop, id_pos))
          callsDat$start[id_pos] <- start
          callsDat$stop[id_pos] <- stop
          current_feature <- F
          cds <- F 
        }
      }  
      
      
      
      
      
      
      combinedCalls <- callsDat %>% filter(start != 0) %>% mutate_all(as.numeric) %>% filter((stop - start) > 50)
      
      
      
      
      gffMain <- readLines(paste("~/phd/RNASeq/genera/",  genera, "/", accession ,".data/gff_files/", accession, ".gff", sep = ""))
      gffMain <- data.frame(text = gffMain)
      genomeInfo <- as.character(gffMain[8,1])
      genomeBuild <- as.character(gffMain[4,1])
      genomeSpecies <- as.character(gffMain[9,1])
      accession_contig <- strsplit(genomeInfo, " ")[[1]][2]
      
      
      
      gffFwd <- combinedCalls%>%mutate(strand = "+",
                                       source = "sraAlignedncRNAExpression",
                                       seqname = accession_contig,
                                       median.val = 0,
                                       feature = "ncRNA",
                                       frame = ".",
                                       attribute = paste("ID=rna_fwd_", row_number(), sep = ""))%>%
        select(seqname, source, feature, start, stop, median.val, strand, frame, attribute)
      
      
      
      
      fileConn<-file(paste("~/phd/RNASeq/genera/",  genera, "/", accession ,".data/gff_files/", sra_id, "_snra_calls.gff", sep = ""))
      writeLines(c("##gff-version 3",
                   "#!gff-spec-version 1.21",
                   "#!processor R script (local) with manual add of top section",
                   genomeBuild,
                   paste("#!genome-build-accession NCBI_Assembly:", "GCA_000017745.1", sep = ""),
                   paste("#!annotation-date ", Sys.Date(), sep = ""),
                   "#!annotation-source snraCalls (local version)",
                   genomeInfo,
                   genomeSpecies), fileConn)
      close(fileConn)
      
      print(paste("Writing to ~/phd/RNASeq/genera/",  genera, "/", accession ,".data/gff_files/", sra_id, "_snra_calls.gff", sep = ""))
      write.table(x = gffFwd, file = paste("~/phd/RNASeq/genera/",  genera, "/", accession ,".data/gff_files/", sra_id, "_snra_calls.gff", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t", append = T)
      
      
      # plotDat <- plotDat %>% mutate(value = ifelse(is.na(type), value, 0))
      # 
      # 
      # outDat <- plotDat %>% select(value)
      # 
      # 
      # write.table(outDat, paste("~/phd/RNASeq/genera/",  genera, "/", accession ,".data/plot_files/", file_name, "_ncRNA_no_buffer.plot", sep = ""), quote = F, row.names = F, col.names = F)
    }
    
    
  }
}