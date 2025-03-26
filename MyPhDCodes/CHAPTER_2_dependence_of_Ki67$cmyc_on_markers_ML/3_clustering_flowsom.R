setwd("C:/Users/rajim/Desktop/PHDcodes")
source("MyPhDCodes/UsedSourceCodes/mass_cytometry_functions.R")
library(sweep)

set.seed(45)
fSOM <- FlowSOM('data/other_data_for_20_cll/concattransformedPlayR.fcs',
compensate = FALSE, transform = FALSE,
scale = FALSE,
colsToUse = c(1:29), xdim = 10, ydim = 10,
nClus = 10)

results<- FlowSOM::GetMetaclusterMFIs(fsom = fSOM)
plot3<- pheatmap::pheatmap(t(results),  
                           border_color = "black", 
                           #scale="row",
                           fontsize=10, 
                           #fontface= "bold",
                           
                   #annotation_row = rownames(average),
                   #treeheight_row = 0,
                   treeheight_col = 0,
                   #cutree_rows = 2,
                    #show_colnames = F,
                    #show_rownames = F,
                   #kmeans_k= 4,
                    #color = c("blue", "yellow","red"),
                    #breaks = c(-4, -0.5, 1, 4), 
                    main = "Cluster Analysis of Cells sampled from all 18 Samples"
                    #main = "Cluster Analysis of Cells"
                  )
#print(plot3)
#View(results)
percentage = FlowSOM::GetPercentages(fsom= fSOM, level = "metaclusters")
print(percentage)
barplot(percentage, col = "skyblue", 
        main = "proportions of the 10 clusters",
        border = "black", ylim = c(0, max(percentage) + 0.5),
        )
box(lwd= 1)
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 2) 

df_minmax <- sweep(t(results), 1, apply(t(results), 1, min), FUN = "-")
df_minmax <- sweep(df_minmax, 1, apply(df_minmax, 1, max), FUN = "/")
View(df_minmax)
plot4<- pheatmap::pheatmap(df_minmax,  
                           border_color = "black", 
                           #scale="row",
                           fontsize=10, 
                           #fontface= "bold",
                           
                           #annotation_row = rownames(average),
                           #treeheight_row = 0,
                           treeheight_col = 0,
                           #cutree_rows = 2,
                           #show_colnames = F,
                           #show_rownames = F,
                           #kmeans_k= 4,
                           #color = c("blue", "yellow","red"),
                           #breaks = c(-4, -0.5, 1, 4), 
                           main = "Cluster Analysis of Cells sampled from all 18 Samples"
                           #main = "Cluster Analysis of Cells"
)#View(percentage)
#playR = read.csv("data/other_data_for_20_cll/concattransformedPlayR.csv")
gene_names <- c("KLF10", "MYC", "CD83", "NR4A2", "PPP1R15A", "JUND", "FOSB", "RGS1", 
                "KLF6", "TCL1A", "DDIT3", "YPEL5", "HSPA5", "TXNIP", "RGS2", 
                "EGR1", "CD69", "FOS", "Actin", "CD19", "CD5", "IgD", "CD21", "CD27", "CD22",
                "Ki67", "CD38", "IgM", "CXCR4")

playR <- fcs_to_dataframe('data/other_data_for_20_cll/concattransformedPlayR.fcs', gene_names, transform_or_not = FALSE)
xyz = FlowSOM::GetMetaclusters(fsom = fSOM)
print(table(xyz))



playR["label"] = FlowSOM::GetMetaclusters(fsom = fSOM)
write.csv(playR, "data/other_data_for_20_cll/concattransformedPlayRWithLabel.csv")
