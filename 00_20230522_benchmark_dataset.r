library(readxl)
library(utils)

is_unix=.Platform$OS.type == "unix"
if(is_unix) {
  source('/home/shengq2/program/scDemultiplex_analysis/common.r')
} else {
  source('C:/Users/sheng/Programs/scDemultiplex_analysis/common.r')
}

#https://www.biorxiv.org/content/biorxiv/early/2023/01/16/2022.12.20.521313.full.pdf
#https://github.com/Oshlack/hashtag-demux-paper


if(!dir.exists("hashtag-demux-paper")){
  system("git clone https://github.com/Oshlack/hashtag-demux-paper")
}

if(1){
  data_dir = paste0(root_dir, "hashtag-demux-paper/data/BAL_data/")

  #For each of the batches need a list of the HTOs and the associated genetic donors.
  donor_list = list(
    "batch1" = list("donor_A" = "Human-HTO-3", "donor_B" = "Human-HTO-1", "donor_C" = "Human-HTO-4", "donor_D" = "Human-HTO-2", "donor_E" = "Human-HTO-8", "donor_F" = "Human-HTO-6", "donor_G" = "Human-HTO-5", "donor_H" = "Human-HTO-7", "doublet" = "Doublet", "unassigned" = "Negative"),
    "batch2" = list("donor_A" = "Human-HTO-7", "donor_B" = "Human-HTO-13", "donor_C" = "Human-HTO-15", "donor_D" = "Human-HTO-12", "donor_E" = "Human-HTO-6", "donor_F" = "Human-HTO-9", "donor_G" = "Human-HTO-14", "donor_H" = "Human-HTO-10", "doublet" = "Doublet", "unassigned" = "Negative"),
    "batch3" = list("donor_A" = "Human-HTO-6", "donor_B" = "Human-HTO-10", "donor_C" = "Human-HTO-14", "donor_D" = "Human-HTO-13", "donor_E" = "Human-HTO-15", "donor_F" = "Human-HTO-7", "donor_G" = "Human-HTO-12", "donor_H" = "Human-HTO-9", "doublet" = "Doublet", "unassigned" = "Negative")
  )

  batch="batch1"
  capture="c1"
  for(batch in names(donor_list)){
    for(capture in c("c1", "c2")){
      cur_sample = paste0(batch, "_", capture)

      sample_folder=paste0(root_dir, cur_sample)
      if(!dir.exists(sample_folder)){
        dir.create(sample_folder)
      }
      setwd(sample_folder)
      
      obj_file = paste0(cur_sample, ".obj.rds")
      if(file.exists(obj_file)){
        next
      }
      
      print(cur_sample)
      
      counts <- read.csv(paste0(data_dir, cur_sample, "_hto_counts.csv"), check.names = FALSE, row.names = 1)
      rownames(counts)<-gsub(" ","_", rownames(counts))
      donors <- read.csv(paste0(data_dir, cur_sample, "_donors.csv"), row.names = 1)
      donors$genetic_donor <- gsub(" ","_", donors$genetic_donor)

      save_to_matrix(counts=as.sparse(counts), target_folder="data")

      rds_file=paste0(cur_sample, ".counts.rds")
      saveRDS(counts, rds_file)

      obj <- scDemultiplex:::read_hto(rds_file)
      
      stopifnot(all(donors$Barcode == colnames(obj)))

      obj$ground_truth=donors$genetic_donor

      obj<-hto_umap(obj)
      obj<-RunTSNE(obj, dims = 1:nrow(obj), perplexity = 100, check_duplicates = FALSE, verbose = FALSE)
      saveRDS(obj, obj_file)  
    }
  }
}

if(1){
  #no doublets and negatives
  #https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-019-0433-8/MediaObjects/41592_2019_433_MOESM3_ESM.xlsx
  library(readxl)
  library(utils)

  cur_sample="barnyard"
  sample_folder=paste0(root_dir, cur_sample)
  if(!dir.exists(sample_folder)){
    dir.create(sample_folder)
  }
  setwd(sample_folder)

  obj_file = paste0(cur_sample, ".obj.rds")
  if(!file.exists(obj_file)){
    print(cur_sample)

    hto_file = '41592_2019_433_MOESM3_ESM.xlsx'
    if(!file.exists(hto_file)){
      download.file('https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-019-0433-8/MediaObjects/41592_2019_433_MOESM3_ESM.xlsx', hto_file)
    }

    hs_data <- read_excel(hto_file, sheet="POC_Nuc-hs_BarMatrix_MetaData", skip=1)
    hs_data = hs_data[hs_data$nUMI_Bar > 0,,drop=F]

    mm_data <- read_excel(hto_file, sheet="POC_Nuc-mm_BarMatrix_MetaData", skip=1)
    mm_data = mm_data[mm_data$nUMI_Bar > 0,,drop=F]

    hs_data$CellID = paste0("hs_", hs_data$CellID)
    mm_data$CellID = paste0("mm_", mm_data$CellID)

    all_data = rbind(hs_data, mm_data)
    counts = data.frame(all_data[,paste0('Bar',c(1:12))])
    rownames(counts) = all_data$CellID
    counts = t(counts)

    save_to_matrix(counts=as.sparse(counts), target_folder="data")

    rds_file=paste0(cur_sample, ".counts.rds")
    saveRDS(counts, rds_file)

    obj <- scDemultiplex:::read_hto(rds_file)

    stopifnot(all(all_data$CellID == colnames(obj)))
    obj$ground_truth=all_data$CellType

    obj<-hto_umap(obj)
    obj<-RunTSNE(obj, dims = 1:nrow(obj), perplexity = 100, check_duplicates = FALSE, verbose = FALSE)
    saveRDS(obj, obj_file)  
  }
}

if(1){
  library(utils)

  cur_sample="pbmc8"
  sample_folder=paste0(root_dir, cur_sample)
  if(!dir.exists(sample_folder)){
    dir.create(sample_folder)
  }
  setwd(sample_folder)

  obj_file = paste0(cur_sample, ".obj.rds")
  if(!file.exists(obj_file)){
    print(cur_sample)

    hto_file = 'stoeckius_pbmc.RData'
    if(!file.exists(hto_file)){
      download.file('https://github.com/Gartner-Lab/deMULTIplex2/raw/main/data/stoeckius_pbmc.RData', hto_file)
    }

    load(hto_file)
    counts = stoeckius_pbmc
    counts = t(counts)

    #we got the ground_truth from Dr. Zhu, Qin, author of paper deMULTIplex2
    donors = read.table('/home/shengq2/program/scDemultiplex_analysis/pbmc_donor_ids.tsv', sep="\t", header=T, row.names=1, stringsAsFactors = F)
    rownames(donors)<-gsub("-.","", rownames(donors))

    stopifnot(all(colnames(counts) %in% rownames(donors)))

    donors = donors[colnames(counts),]

    ids = table(donors$donor_id)
    ids = ids[ids < 20]

    donors$donor_id[donors$donor_id %in% names(ids)] = "unassigned"

    save_to_matrix(counts=as.sparse(counts), target_folder="data")

    rds_file=paste0(cur_sample, ".counts.rds")
    saveRDS(counts, rds_file)

    obj <- scDemultiplex:::read_hto(rds_file)

    obj$ground_truth=donors$donor_id

    obj<-hto_umap(obj)
    obj<-RunTSNE(obj, dims = 1:nrow(obj), perplexity = 100, check_duplicates = FALSE, verbose = FALSE)
    saveRDS(obj, obj_file)  
  }
}
