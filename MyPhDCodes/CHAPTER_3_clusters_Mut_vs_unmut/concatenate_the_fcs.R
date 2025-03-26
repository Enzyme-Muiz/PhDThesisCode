source("MyPhDCodes/UsedSourceCodes/mass_cytometry_functions.R")


needed_antibody1<- c('CD196','CD19', 'CD5',  'CD38', 'IgD', 'CD11c', 'CD43', 'CD69',
                    'CD21', 'CXCR5', 'CD62L', 'CD27', 'CD22', 'CXCR3', 'CD23', 'CD24', 'CCR7', 
                    'CD20', 'IgM', 'HLA-DR', 'CD49d', 'CXCR4')

number_of_cell = 2000
concat_result <- generate_concat(c("_1_", "_2_", "_3_", "_5_","_6_"), 
                "data/all_cytonorm_normalized/",
                number_of_cell,
                column_needed = needed_antibody1,
                needed_antibody1,
                "yes",
                "concattransformed"
               )