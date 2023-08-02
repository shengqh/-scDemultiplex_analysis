library(readxl)
library(utils)

is_unix=.Platform$OS.type == "unix"
if(is_unix) {
  source('/home/shengq2/program/scDemultiplex_analysis/common.r')
} else {
  source('C:/Users/sheng/Programs/scDemultiplex_analysis/common.r')
}

load_install("scDemultiplex", c('shengqh/cutoff', 'shengqh/scDemultiplex'))

#https://www.biorxiv.org/content/biorxiv/early/2023/01/16/2022.12.20.521313.full.pdf
#https://github.com/Oshlack/hashtag-demux-paper


if(!dir.exists("hashtag-demux-paper")){
  system("git clone https://github.com/Oshlack/hashtag-demux-paper")
}

if(1){
  data_dir = paste0(root_dir, "hashtag-demux-paper/data/BAL_data/")

  bal_donor_map = list(
    "BAL A" = "BAL-01",
    "BAL B" = "BAL-02",
    "BAL C" = "BAL-03",
    "BAL D" = "BAL-04",
    "BAL E" = "BAL-05",
    "BAL F" = "BAL-06",
    "BAL G" = "BAL-07",
    "BAL H" = "BAL-08",
    "BAL I" = "BAL-09",
    "BAL J" = "BAL-10",
    "BAL K" = "BAL-11",
    "BAL L" = "BAL-12",
    "BAL M" = "BAL-13",
    "BAL N" = "BAL-14",
    "BAL O" = "BAL-15",
    "BAL P" = "BAL-16",
    "BAL Q" = "BAL-17",
    "BAL R" = "BAL-18",
    "BAL S" = "BAL-19",
    "BAL T" = "BAL-20",
    "BAL U" = "BAL-21",
    "BAL V" = "BAL-22",
    "BAL W" = "BAL-23",
    "BAL X" = "BAL-24",
    "Doublet" = "Doublet",
    "Negative" = "Negative"
  )

  batch="batch1"
  capture="c1"
  for(batch in c("batch1", "batch2", "batch3")){
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
      rownames(counts)<-gsub(" ","-", rownames(counts))
      stopifnot(all(rownames(counts) %in% bal_donor_map))

      donors <- read.csv(paste0(data_dir, cur_sample, "_donors.csv"), row.names = 1)
      stopifnot(all(donors$genetic_donor %in% names(bal_donor_map)))
      donors$genetic_donor <- unlist(bal_donor_map[donors$genetic_donor])
      stopifnot(all(donors$genetic_donor %in% bal_donor_map))

      save_to_matrix(counts=as.sparse(counts), target_folder="data")

      rds_file=paste0(cur_sample, ".counts.rds")
      saveRDS(counts, rds_file)

      if(use_rlog){
        obj <- CreateSeuratObject(counts = counts, assay="HTO")
        # Normalize HTO data, here we use centered log-ratio (CLR) transformation
        obj <- NormalizeData(obj, assay = "HTO", normalization.method = "LogNormalize")
        DefaultAssay(object = obj) <- "HTO"
      }else{
        obj <- scDemultiplex:::read_hto(rds_file)
      }
      stopifnot(all(obj$nCount_HTO > 0))
      stopifnot(all(donors$Barcode == colnames(obj)))

      obj$ground_truth=donors$genetic_donor

      obj<-hto_umap(obj)
      #obj<-RunTSNE(obj, dims = 1:nrow(obj), perplexity = 100, check_duplicates = FALSE, verbose = FALSE)
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
    stopifnot(all(obj$nCount_HTO > 0))
    stopifnot(all(all_data$CellID == colnames(obj)))
    obj$ground_truth=all_data$CellType

    obj<-hto_umap(obj)
    #obj<-RunTSNE(obj, dims = 1:nrow(obj), perplexity = 100, check_duplicates = FALSE, verbose = FALSE)
    saveRDS(obj, obj_file)  
  }
}

if(0){#pbmc from original paper
  library(readxl)
  library(utils)

  cur_sample="pbmc8"
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

  hto_file = 'GSM2895283_Hashtag-HTO-count.csv.gz'
  if(!file.exists(hto_file)){
    download.file('https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2895nnn/GSM2895283/suppl/GSM2895283_Hashtag-HTO-count.csv.gz', hto_file)
  }
  exp_file = 'GSM2895282_Hashtag-RNA.umi.txt.gz'
  if(!file.exists(exp_file)){
    download.file('https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2895nnn/GSM2895282/suppl/GSM2895282_Hashtag-RNA.umi.txt.gz', exp_file)
  }

  if(!dir.exists("scSplit_paper_data")){
    system("git clone https://github.com/jon-xu/scSplit_paper_data/")
  }
  
  n='A'
  df = NULL
  for(n in c('A', 'B', 'C','D', 'E', 'F', 'G', 'H', 'doublets')){
    barcodes = readLines(paste0("scSplit_paper_data/Table 3/Hashtag/", n))
    cur_df = data.frame("cell"=barcodes, "ground_truth"=ifelse(n == "doublets", "Doublet", paste0("Batch", n)))
    df = rbind(df, cur_df)
  }

  hto_data <- data.frame(fread(hto_file), row.names=1)
  counts <- hto_data[c(1:8),]
  rownames(counts) <- gsub("-.*", "", rownames(counts))

  stopifnot(all(df$cell %in% colnames(counts)))
  counts = counts[,df$cell]

  save_to_matrix(counts=as.sparse(counts), target_folder="data")

  rds_file=paste0(cur_sample, ".counts.rds")
  saveRDS(counts, rds_file)

  obj <- scDemultiplex:::read_hto(rds_file)

  stopifnot(all(df$cell == colnames(obj)))
  obj$ground_truth=df$ground_truth

  obj<-hto_umap(obj)
  obj<-RunTSNE(obj, dims = 1:nrow(obj), perplexity = 100, check_duplicates = FALSE, verbose = FALSE)
  saveRDS(obj, obj_file)  
}


if(1){#pbmc with new annotation
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
    unlink(hto_file)

    counts = stoeckius_pbmc
    counts = t(counts)
    rownames(counts)<-gsub("_","-", rownames(counts))

    #we got the ground_truth from Dr. Zhu, Qin, author of paper deMULTIplex2
    donors = read.table('/home/shengq2/program/scDemultiplex_analysis/pbmc_donor_ids.tsv', sep="\t", header=T, row.names=1, stringsAsFactors = F)
    rownames(donors)<-gsub("-.","", rownames(donors))

    donor_map = list(
      "donor6" = "HTO-A",
      "donor0" = "HTO-B",
      "donor8" = "HTO-C",
      "donor10" = "HTO-D",
      "donor5" = "HTO-E",
      "donor11" = "HTO-F",
      "donor2" = "HTO-G",
      "donor7" = "HTO-H",
      "doublet" = "Doublet",
      "unassigned" = "Negative"
    )

    stopifnot(all(colnames(counts) %in% rownames(donors)))

    donors = donors[colnames(counts),]

    ids = table(donors$donor_id)
    ids = ids[ids < 20]

    donors$donor_id[donors$donor_id %in% names(ids)] = "unassigned"
    donors$hto = unlist(donor_map[donors$donor_id])

    save_to_matrix(counts=as.sparse(counts), target_folder="data")

    rds_file=paste0(cur_sample, ".counts.rds")
    saveRDS(counts, rds_file)

    obj <- scDemultiplex:::read_hto(rds_file)
    stopifnot(all(obj$nCount_HTO > 0))

    obj$ground_truth=donors$hto

    obj<-hto_umap(obj)
    #obj<-RunTSNE(obj, dims = 1:nrow(obj), perplexity = 100, check_duplicates = FALSE, verbose = FALSE)
    saveRDS(obj, obj_file)  
  }
}



