setwd("C:/Users/rajim/Desktop/PHDcodes")
source("MyPhDCodes/UsedSourceCodes/mass_cytometry_functions.R")
library(randomcoloR)
options(
repr.plot.width=15,
repr.plot.height=7

)

playR = read.csv("data/other_data_for_20_cll/concattransformedPlayRWithLabel.csv", row.names = NULL)
#print(head(playR))
playR <- playR[sample(nrow(playR), 2000), ] 
rownames(playR) <- NULL
label = playR$label
print(table(label))
colors1 = randomcoloR::distinctColorPalette(length(unique(label)))
#colors1 = c("green", "cyan", "purple", "yellow", "red", "blue", "orange", "black", "grey", "skyblue")
colors = colors1[label]
print(colors)
print(table(colors))
playRdata <- playR[, !names(playR) %in% c("label", "X.1", "X")]
set.seed(45)
gg<- umap::umap(playRdata, n_neighbours= 20, min_dist= 0.25)

df1<- gg$layout
colnames(df1)<- c("umap1", "umap2")
clusters<- label
df1<- cbind(df1, clusters)
first<- aggregate(umap1 ~ label, data = df1, mean)
second<- aggregate(umap2 ~ label, data = df1, mean)
coordinate<- cbind(first, second$umap2)
colnames(coordinate)<- c("clusters", "umap1", "umap2")


plot(gg$layout, main= paste0("UMAP of sampled CLL cells from all samples\n with ", " FLOWSOM clusters"),
     xlab=substitute(paste(bold('umap_1'))), 
     ylab=substitute(paste(bold('umap_2'))),
     col=colors, 
     pch=19,
     cex= 1.5,
     cex.main=2.0
) 


head(gg$layout)
legend("topright", legend= c(1:length(unique(label))), col=colors1, lwd=10:15, cex=0.8)
box(lwd= 5)
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 2) 
for (i in c(1:dim(coordinate)[1])){
  text(x=coordinate[i, 2], y=coordinate[i, 3], coordinate[i, 1], cex=2, font=3)}
