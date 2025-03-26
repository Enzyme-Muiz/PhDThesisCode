source("MyPhDCodes/UsedSourceCodes/mass_cytometry_functions.R")
options(repr.plot.width=12, repr.plot.height=7)

needed_antibody1<- c('CD196','CD19', 'CD5',  'CD38', 'IgD', 'CD11c', 'CD43', 'CD69',
                    'CD21', 'CXCR5', 'CD62L', 'CD27', 'CD22', 'CXCR3', 'CD23', 'CD24', 'CCR7', 
                    'CD20', 'IgM', 'HLA-DR', 'CD49d', 'CXCR4')


result <- colmedian_input_marker_of_fcs("data/all_cytonorm_normalized/",
                                    c("_1_", "_2_", "_3_", "_5_","_6_"),
                                    "yes",
                                    needed_antibody1, 
                                    needed_antibody1
)

write.csv(result, "median_markers_expressions.csv")