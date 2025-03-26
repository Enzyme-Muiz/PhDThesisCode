library(CytofRUV)
library(CATALYST)
library(flowCore)
library(ggplot2)
library(readxl)
library(ruv)
library(purrr)
library(FlowSOM)
library(SummarizedExperiment)
library(ConsensusClusterPlus)
library(SingleCellExperiment)
library(shiny)
library(shinyjs)
library(shinydashboard)
library(writexl)
library(ComplexHeatmap)
library(shinycssloaders)
library(readxl)
library(stringr)
source("make_fcs_file.R")
# source("normalize_function.R")
loading_data <- function(clusters_nb = 20, seed = 1234, rep_samples, k_value, copy_results_here = "fcs_untransformed") {
  try({
    output_dir <- "CytofRUV_output"
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    2 + 2
    wd_data <- file.path(getwd(), output_dir)
    print(wd_data)
    metadata_filename <- "Metadata.xlsx"
    panel_filename <- "Panel.xlsx"
    seed <- 1234
    clusters_nb <- clusters_nb
    print("abc2")
    data <- load_data(wd_data, metadata_filename, panel_filename, cofactor = 5)
    print("abc1")
    print(data$daf)
    print(getwd())
    xx <- read_excel("Panel.xlsx")
    print(xx$antigen)
    # xx$fcs_colname
    print(data$lineage_markers)
    data$daf <- cluster_data(data$daf, seed,
      markers_to_use = data$lineage_markers,
      clusters_nb
    )

    daf <- data$daf
    md <- data$md
    seed <- seed
    n_subset_marker_specific <- 1000
    print("abc")
    TSNE_subset <- 1000
    print("Running TSNE")
    daf <- CATALYST::runDR(daf, dr = "TSNE", cells = TSNE_subset)


    UMAP_subset <- 1000
    print("Running UMAP")
    daf <- runDR(daf, "UMAP", cells = UMAP_subset)

    n_subset <- 1000
    sub_daf <- daf[, sample(ncol(daf), n_subset)]

    panel <- data$panel




    print("-1")
    dir_name_norm_data <- "CytofRUV_Norm_data_HC2_all_cl_20"
    print("0")
    print(data$daf)
    clusters <- paste0("meta", as.character(clusters_nb))
    raw_data <- data.frame(sample = data$daf$sample_id, cluster = cluster_ids(data$daf, clusters), t(SummarizedExperiment::assay(data$daf, "exprs")))
    print("1")
    colnames(raw_data) <- gsub("^X", "", colnames(raw_data))
    print("2")
    rep_samples <- list(rep_samples)
    print("3")
    cluster_list_rep_samples <- list(seq(1, as.numeric(clusters_nb)))
    print(cluster_list_rep_samples)
    print("4")
    k_value <- k_value
    seed <- 1234
    normalise_data(data = data, raw_data = raw_data, rep_samples = rep_samples, norm_clusters = cluster_list_rep_samples, k = k_value, num_clusters = clusters_nb, wd_data = wd_data, dir_norm_data = dir_name_norm_data)
    print(getwd())

    setwd("..")
    unarcsinh_fcs_files(
      "CytofRUV_output//CytofRUV_Norm_data_HC2_all_cl_20//",
      "CytofRUV_output//CytofRUV_Norm_data_HC2_all_cl_20//"
    )
    path <- "CytofRUV_output//CytofRUV_Norm_data_HC2_all_cl_20//"

    list_in_fcs_untransformed <- list.files(path, full.names = TRUE)
    list2 <- list_in_fcs_untransformed[str_detect(list_in_fcs_untransformed, ".fcs")]
    file.rename(list2, paste0(path, str_match(basename(list2), "(Batch_.+)")[, 1]))
    list_in_fcs_untransformed <- list.files(path, full.names = TRUE)
    list2 <- list_in_fcs_untransformed[str_detect(list_in_fcs_untransformed, ".fcs")]
    file.copy(list2, copy_results_here, overwrite = TRUE)
  })
}



# CytofRUV::launch_Shiny(daf)



# dir_name_norm_data="CytofRUV_Norm_data_HC2_all_cl_20"
# raw_data <- data.frame(sample = data$daf$sample_id, cluster=cluster_ids(data$daf,"meta20"), t(SummarizedExperiment::assay(data$daf,"exprs")))
# colnames(raw_data) <- gsub("^X", "",  colnames(raw_data))
# rep_samples=list(c("HC2_B1","HC2_B2"))
# cluster_list_rep_samples <- list(seq(1,20))
# k_value <- 5
# seed=1234

# normalise_data(data=data,raw_data=raw_data,rep_samples=rep_samples, norm_clusters=cluster_list_rep_samples, k=k_value, num_clusters=clusters_nb,wd_data=wd_data,dir_norm_data=dir_name_norm_data)
normalize_by_cytofruv <- function(rep_samples) {
  try({
    print("-1")
    dir_name_norm_data <- "CytofRUV_Norm_data_HC2_all_cl_20"
    print("0")
    print(data)
    raw_data <- data.frame(sample = data$daf$sample_id, cluster = cluster_ids(data$daf, "meta20"), t(SummarizedExperiment::assay(data$daf, "exprs")))
    print("1")
    colnames(raw_data) <- gsub("^X", "", colnames(raw_data))
    print("2")
    rep_samples <- list(rep_samples)
    print("3")
    cluster_list_rep_samples <- list(seq(1, 20))
    print("4")
    k_value <- 5
    seed <- 1234
    # normalise_data(data=data,raw_data=raw_data,rep_samples=rep_samples, norm_clusters=cluster_list_rep_samples, k=k_value, num_clusters=clusters_nb,wd_data=wd_data,dir_norm_data=dir_name_norm_data)
    print(getwd())
  })
}













# library(CytofRUV)
# library(CATALYST)
# library(flowCore)
# library(ggplot2)
# library(readxl)
# library(ruv)
# library(purrr)
# library(FlowSOM)
# library(SummarizedExperiment)
# library(ConsensusClusterPlus)
# library(SingleCellExperiment)
# library(shiny)
# library(shinyjs)
# library(shinydashboard)
# library(writexl)
# library(ComplexHeatmap)
# library(shinycssloaders)
# library(readxl)
# library(stringr)


# setwd("abc")
# rep_samples<- c()
# for (i in c(1:6)){
#     anchor_vali<- read_excel("Metadata.xlsx")$sample_id[str_detect(read_excel("Metadata.xlsx")$sample_id, "anchor|vali")]
#     results<- anchor_vali[str_detect(anchor_vali, as.character(i))]
#     if (!identical(results, character(0)) ) {
#         #print(results)
#         rep_samples<- append(rep_samples, results[1])
#     }
# }
# #rep_s<- rep_samples

# #rep_samples <- append(rep_s[1], paste0("vali", str_match(rep_s[1], "\\d")[ ,1]))

# setwd("..")




# #setwd("C://Users/oaona/OneDrive/Documents")
# #rep_samples=list(c("HC2_B1","HC2_B2"))
# rep_samples = list(rep_samples)
# output_dir= "abc"
# if (!dir.exists(output_dir)){
#     dir.create(output_dir)
# }
# wd_data=file.path(getwd(),output_dir)
# print(wd_data)
# metadata_filename="Metadata.xlsx"
# panel_filename="Panel.xlsx"
# seed=1234
# clusters_nb=2
# data <- load_data(wd_data,metadata_filename,panel_filename, cofactor = 5)
# print(data$lineage_markers)
# data$daf<- cluster_data(data$daf,seed,
#     markers_to_use=data$lineage_markers,
#     clusters_nb)

# daf<- data$daf
# md<- data$md
# seed=seed
# n_subset_marker_specific <- 1000
# TSNE_subset <- 1000
# print("Running TSNE")
# daf <- CATALYST::runDR(daf, dr = "TSNE", cells = TSNE_subset)
# UMAP_subset <- 1000
# print("Running UMAP")
# daf <- runDR(daf, "UMAP", cells = UMAP_subset)
# n_subset <- 1000
# sub_daf <- daf[, sample(ncol(daf), n_subset)]
# panel<- data$panel
# dir_name_norm_data="CytofRUV_Norm_data_HC2_all_cl_20"
# clusters = paste0("meta", as.character(clusters_nb))
# raw_data <- data.frame(sample = data$daf$sample_id, cluster=cluster_ids(data$daf,clusters), t(SummarizedExperiment::assay(data$daf,"exprs")))
# colnames(raw_data) <- gsub("^X", "",  colnames(raw_data))
# cluster_list_rep_samples <- list(seq(1,as.numeric(clusters_nb)))
# k_value <- 5
# seed=1234
# normalise_data(data=data,raw_data=raw_data,rep_samples=rep_samples, norm_clusters=cluster_list_rep_samples, k=k_value, num_clusters=clusters_nb,wd_data=wd_data,dir_norm_data=dir_name_norm_data)
# setwd("..")













































































# output_dir="CytofRUV_output2"
# if (!dir.exists(output_dir)){
# dir.create(output_dir)
# }
# wd_data=file.path(getwd(),output_dir)
# write.FCS(x=CytofRUV::A1,filename =file.path(wd_data,"A1.fcs"))
# write.FCS(x=CytofRUV::A2,filename = file.path(wd_data,"A2.fcs"))
# write.FCS(x=CytofRUV::Run3_A1,filename = file.path(wd_data,"Run3_A1.fcs"))
# write.FCS(x=CytofRUV::Run3_A2,filename = file.path(wd_data,"Run3_A2.fcs"))
# write_xlsx(x=CytofRUV::md,path = file.path(wd_data,"Metadata.xlsx"))
# write_xlsx(x=CytofRUV::panel,path = file.path(wd_data,"Panel.xlsx"))









# #rep_samples=list(c("HC2_B1","HC2_B2"))

# output_dir="CytofRUV_output3"
# if (!dir.exists(output_dir)){
# dir.create(output_dir)
# }
# wd_data=file.path(getwd(),output_dir)
# metadata_filename="Metadata.xlsx"
# panel_filename="Panel.xlsx"
# seed=1234
# clusters_nb=20
# data=load_data(wd_data,metadata_filename,panel_filename)
# data$daf=cluster_data(data$daf,seed,markers_to_use=data$lineage_markers,clusters_nb)
# daf=data$daf
# md=data$md
# seed=1234
# n_subset_marker_specific <- 1000
# set.seed(seed)
# TSNE_subset <- 1000
# daf <- CATALYST::runDR(daf, dr = "TSNE", cells = TSNE_subset)
# UMAP_subset <- 2000
# daf <- runDR(daf, "UMAP", cells = UMAP_subset)
# n_subset <- 1000
# sub_daf <- daf[, sample(ncol(daf), n_subset)]
# panel=data$panel
# dir_name_norm_data="CytofRUV_Norm_data_HC2_all_cl_20"
# raw_data <- data.frame(sample = data$daf$sample_id, cluster=cluster_ids(data$daf,"meta20"), t(SummarizedExperiment::assay(data$daf,"exprs")))
# colnames(raw_data) <- gsub("^X", "",  colnames(raw_data))
# rep_samples=list(c("HC2_B1","HC2_B2"))
# cluster_list_rep_samples <- list(seq(1,20))
# k_value <- 5
# seed=1234
# normalise_data(data=data,raw_data=raw_data,rep_samples=rep_samples, norm_clusters=cluster_list_rep_samples, k=k_value, num_clusters=clusters_nb,wd_data=wd_data,dir_norm_data=dir_name_norm_data)
