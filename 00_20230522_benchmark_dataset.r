is_unix=.Platform$OS.type == "unix"
if(is_unix) {
  source('/home/shengq2/program/scDemultiplex_analysis/common.r')
} else {
  source('C:/Users/sheng/Programs/scDemultiplex_analysis/common.r')
}

#https://www.biorxiv.org/content/biorxiv/early/2023/01/16/2022.12.20.521313.full.pdf
#https://github.com/Oshlack/hashtag-demux-paper

data_dir = paste0(root_dir, "/hashtag-demux-paper/data/")

#For each of the batches need a list of the HTOs and the associated genetic donors.
donor_list = list(
  "batch1" = list("donor_A" = "Human-HTO-3", "donor_B" = "Human-HTO-1", "donor_C" = "Human-HTO-4", "donor_D" = "Human-HTO-2", "donor_E" = "Human-HTO-8", "donor_F" = "Human-HTO-6", "donor_G" = "Human-HTO-5", "donor_H" = "Human-HTO-7", "doublet" = "Doublet", "unassigned" = "Negative"),
  "batch2" = list("donor_A" = "Human-HTO-7", "donor_B" = "Human-HTO-13", "donor_C" = "Human-HTO-15", "donor_D" = "Human-HTO-12", "donor_E" = "Human-HTO-6", "donor_F" = "Human-HTO-9", "donor_G" = "Human-HTO-14", "donor_H" = "Human-HTO-10", "doublet" = "Doublet", "unassigned" = "Negative"),
  "batch3" = list("donor_A" = "Human-HTO-6", "donor_B" = "Human-HTO-10", "donor_C" = "Human-HTO-14", "donor_D" = "Human-HTO-13", "donor_E" = "Human-HTO-15", "donor_F" = "Human-HTO-7", "donor_G" = "Human-HTO-12", "donor_H" = "Human-HTO-9", "doublet" = "Doublet", "unassigned" = "Negative")
)

batch="batch2"
capture="c2"
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
    donors <- read.csv(paste0(data_dir, cur_sample, "_donors.csv"), row.names = 1)

    save_to_matrix(counts=as.sparse(counts), target_folder="data")

    rds_file=paste0(cur_sample, ".counts.rds")
    saveRDS(counts, rds_file)

    obj <- scDemultiplex:::read_hto(rds_file)
    
    stopifnot(all(rownames(donors) == colnames(obj)))

    donor_map = donor_list[[batch]]
    obj$genetic_HTO=unlist(donor_map[donors$V1])

    obj<-hto_umap(obj)
    obj<-RunTSNE(obj, dims = 1:nrow(obj), perplexity = 100, check_duplicates = FALSE, verbose = FALSE)
    saveRDS(obj, obj_file)  
  }
}
