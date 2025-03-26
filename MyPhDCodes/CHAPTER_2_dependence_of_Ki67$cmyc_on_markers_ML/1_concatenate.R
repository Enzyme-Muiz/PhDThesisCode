source("MyPhDCodes/UsedSourceCodes/mass_cytometry_functions.R")


# Vector with INS identifiers
ins_ids <- c("INS14", "INS17", "INS5", "INS19", "INS18", "INS20", "INS16", "INS7", 
             "INS24","INS3", "INS26", "INS8", "INS25", "INS11", "INS13", 
             "INS4", "INS6", "INS27", "INS1", "CD19", "CD5", "IgD", "CD21", "CD27", "CD22",
             "Ki67", "CD38", "IgM", "CXCR4")

# Vector with gene names
gene_names <- c("KLF10", "MYC", "CD83", "NR4A2", "PPP1R15A", "JUND", "FOSB", "RGS1", 
                "KLF6", "TCL1A", "DDIT3", "YPEL5", "HSPA5", "TXNIP", "RGS2", 
                "EGR1", "CD69", "FOS", "Actin", "CD19", "CD5", "IgD", "CD21", "CD27", "CD22",
                "Ki67", "CD38", "IgM", "CXCR4")

# # Vector with INS identifiers
# ins_ids <- c( "CD19", "CD5", "IgD", "CD21", "CD27", "CD22",
#              "Ki67", "CD38", "IgM", "CXCR4")

# # Vector with gene names
# gene_names <- c("CD19", "CD5", "IgD", "CD21", "CD27", "CD22",
#              "Ki67", "CD38", "IgM", "CXCR4")

number_of_cell = 4000
concat_result <- generate_concat(NULL, 
                "data/20 CLL/",
                number_of_cell,
                column_needed = ins_ids,
                gene_names,
                "yes",
                "concattransformedPlayR", 
                "fcs"
               )