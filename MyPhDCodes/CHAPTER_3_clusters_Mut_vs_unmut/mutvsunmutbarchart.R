source("MyPhDCodes/UsedSourceCodes/mass_cytometry_functions.R")
options(repr.plot.width=12, repr.plot.height=7)

needed_antibody1<- c('CD196','CD19', 'CD5',  'CD38', 'IgD', 'CD11c', 'CD43', 'CD69',
                    'CD21', 'CXCR5', 'CD62L', 'CD27', 'CD22', 'CXCR3', 'CD23', 'CD24', 'CCR7', 
                    'CD20', 'IgM', 'HLA-DR', 'CD49d', 'CXCR4')


image1 <- boxplot_groupby_input_marker_of_fcs("data/all_cytonorm_normalized/",
                                    c("_6_","_1_", "_3_", "_2_"),
                                    "yes",
                                    needed_antibody1, 
                                    FALSE,
                                    TRUE,
                                    "Mutation Status",
                                    "The median intensity expression of cells by markers \n using 51 CLL samples"
)

print(image1)