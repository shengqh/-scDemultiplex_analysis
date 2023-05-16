source('/home/shengq2/program/scDemultiplex_analysis/common.r')

setwd("/panfs/accrepfs.vampire/nobackup/h_cqs/collaboration/20230301_scrna_hto/hashtag-demux-paper")

library(here)

batch1_c1_counts <- read.csv(here("data", "batch1_c1_hto_counts.csv"), check.names = FALSE, row.names = 1)
batch1_c2_counts <- read.csv(here("data", "batch1_c2_hto_counts.csv"), check.names = FALSE, row.names = 1)

batch2_c1_counts <- read.csv(here("data", "batch2_c1_hto_counts.csv"), check.names = FALSE, row.names = 1)
batch2_c2_counts <- read.csv(here("data", "batch2_c2_hto_counts.csv"), check.names = FALSE, row.names = 1)

batch3_c1_counts <- read.csv(here("data", "batch3_c1_hto_counts.csv"), check.names = FALSE, row.names = 1)
batch3_c2_counts <- read.csv(here("data", "batch3_c2_hto_counts.csv"), check.names = FALSE, row.names = 1)

batch1_c1_donors <- read.csv(here("data", "batch1_c1_donors.csv"), row.names = 1)
batch1_c2_donors <- read.csv(here("data", "batch1_c2_donors.csv"), row.names = 1)

batch2_c1_donors <- read.csv(here("data", "batch2_c1_donors.csv"), row.names = 1)
batch2_c2_donors <- read.csv(here("data", "batch2_c2_donors.csv"), row.names = 1)

batch3_c1_donors <- read.csv(here("data", "batch3_c1_donors.csv"), row.names = 1)
batch3_c2_donors <- read.csv(here("data", "batch3_c2_donors.csv"), row.names = 1)

seu1_c1 <- CreateSeuratObject(counts = batch1_c1_counts, assay = "HTO")
seu1_c2 <- CreateSeuratObject(counts = batch1_c2_counts, assay = "HTO")

seu2_c1 <- CreateSeuratObject(counts = batch2_c1_counts, assay = "HTO")
seu2_c2 <- CreateSeuratObject(counts = batch2_c2_counts, assay = "HTO")

seu3_c1 <- CreateSeuratObject(counts = batch3_c1_counts, assay = "HTO")
seu3_c2 <- CreateSeuratObject(counts = batch3_c2_counts, assay = "HTO")

seu1_c1$Barcode <- colnames(seu1_c1)
seu1_c1$capture <- "capture_1"
seu1_c2$Barcode <- colnames(seu1_c2)
seu1_c2$capture <- "capture_2"

seu2_c1$Barcode <- colnames(seu2_c1)
seu2_c1$capture <- "capture_1"
seu2_c2$Barcode <- colnames(seu2_c2)
seu2_c2$capture <- "capture_2"

seu3_c1$Barcode <- colnames(seu3_c1)
seu3_c1$capture <- "capture_1"
seu3_c2$Barcode <- colnames(seu3_c2)
seu3_c2$capture <- "capture_2"

seu1_c1$genetic_donor <- batch1_c1_donors
seu1_c2$genetic_donor <- batch1_c2_donors

seu2_c1$genetic_donor <- batch2_c1_donors
seu2_c2$genetic_donor <- batch2_c2_donors

seu3_c1$genetic_donor <- batch3_c1_donors
seu3_c2$genetic_donor <- batch3_c2_donors

seu1 <- merge(seu1_c1, seu1_c2)
seu2 <- merge(seu2_c1, seu2_c2)
seu3 <- merge(seu3_c1, seu3_c2)

#For each of the batches need a list of the HTOs and the associated genetic donors.
donor_list_batch1 <- list("donor_A" = "Human-HTO-3", "donor_B" = "Human-HTO-1", "donor_C" = "Human-HTO-4", "donor_D" = "Human-HTO-2", "donor_E" = "Human-HTO-8", "donor_F" = "Human-HTO-6", "donor_G" = "Human-HTO-5", "donor_H" = "Human-HTO-7", "doublet" = "Doublet", "unassigned" = "Negative")
donor_list_batch2 <- list("donor_A" = "Human-HTO-7", "donor_B" = "Human-HTO-13", "donor_C" = "Human-HTO-15", "donor_D" = "Human-HTO-12", "donor_E" = "Human-HTO-6", "donor_F" = "Human-HTO-9", "donor_G" = "Human-HTO-14", "donor_H" = "Human-HTO-10",
                          "doublet" = "Doublet", "unassigned" = "Negative")
donor_list_batch3 <- list("donor_A" = "Human-HTO-6", "donor_B" = "Human-HTO-10", "donor_C" = "Human-HTO-14", "donor_D" = "Human-HTO-13", "donor_E" = "Human-HTO-15", "donor_F" = "Human-HTO-7", "donor_G" = "Human-HTO-12", "donor_H" = "Human-HTO-9", "doublet" = "Doublet", "unassigned" = "Negative")

seu1$genetic_HTO=unlist(donor_list_batch1[seu1$genetic_donor])
seu2$genetic_HTO=unlist(donor_list_batch2[seu2$genetic_donor])
seu3$genetic_HTO=unlist(donor_list_batch3[seu3$genetic_donor])

sample_objs=list("batch1"=seu1, "batch2" = seu2, "batch3" = seu3)
for(cur_sample in names(sample_objs)){  
  sample_folder=paste0(root_dir, cur_sample)
  if(!dir.exists(sample_folder)){
    dir.create(sample_folder)
  }
  setwd(sample_folder)

  obj = sample_objs[[cur_sample]]
  counts <- GetAssayData(obj, assay="HTO", slot = "counts")
  
  save_to_matrix(counts=counts, target_folder="data")
  
  rds_file=paste0(cur_sample, ".counts.rds")
  saveRDS(counts, rds_file)
  
  newobj <- scDemultiplex:::read_hto(rds_file)
  stopifnot(all(colnames(newobj) == colnames(obj)))
  newobj$genetic_HTO = obj$genetic_HTO
  newobj<-hto_umap(newobj)
  saveRDS(newobj, paste0(cur_sample, ".obj.rds"))  
}
