source("MyPhDCodes/UsedSourceCodes/mass_cytometry_functions.R")


pqr<- cluster_proportion_from_concat(
    "concattransformed.fcs",
    FALSE,
    concat_result$cell_lengths,
    20,
    concat_result$columns_length,
    concat_result$file_names
)

write.csv(pqr$proportion, "clusters_proportions.csv")