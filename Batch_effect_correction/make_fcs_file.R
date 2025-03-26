for (i in c("MASS","ggplot2", "CytoNorm", "reshape2", "flowCore", "cluster", "stringr",
           "Biobase", "flowViz", "Rtsne", "reticulate", "R.utils", "emdist", "FlowSOM","umap", "dplyr",
            "matrixStats", "combinat", "kernlab", "gridExtra", "writexl"
           )){
    
     suppressWarnings(suppressMessages(library(i, character.only = TRUE)))
     }   



#### function that makes fcs file from dataframe
make_fcs_from_df<- function(df, storage_folder, filename){
    meta <- data.frame(name=colnames(df), desc=colnames(df))
    meta$range <- apply(apply(df,2,range),2,diff)
    meta$minRange <- apply(df,2,min)
    meta$maxRange <- apply(df,2,max)
    ff <- new("flowFrame", exprs= data.matrix(df), parameters=Biobase::AnnotatedDataFrame(meta))
    if (!dir.exists(storage_folder)){
    dir.create(storage_folder)
    } else {}
    #flowCore::write.FCS(ff, file.path(storage_folder, paste0(filename, ".fcs")))
    flowCore::write.FCS(ff, file.path(storage_folder, filename))
    

    
}





###function that changes untransform all fcs files in the given path

unarcsinh_fcs_files<- function(path= ".", storage_folder){
    if (length(list.files(path, pattern = "\\.fcs", full.names = TRUE))>= 1){
                for (fcs_files in list.files(path, pattern = "\\.fcs", full.names = TRUE)){
                    fcs_file<- flowCore::read.FCS(fcs_files)
                    fcs_expression_data<- flowCore::exprs(fcs_file)   
                    Time<- fcs_expression_data[ , "Time"] ##Time column need not arcsinh transformed
                    fcs_expression_data <- subset(fcs_expression_data, select = -c(Time))
                    fcs_expression_data<- cytofTransform.reverse(fcs_expression_data)
                    fcs_expression_data<- cbind(Time, fcs_expression_data)
                    make_fcs_from_df(fcs_expression_data, storage_folder, basename(fcs_files))
                    #print(basename(fcs_files))
                }
        }
    else{}
}