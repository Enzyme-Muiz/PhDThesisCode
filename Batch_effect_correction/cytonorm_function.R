PlotOverviewCV2 <- function(fsom, cv_res, max_cv = 2.5, show_cv = 1.5){
    cvs <- cv_res$cvs
    pctgs <- cv_res$pctgs
    nMetaClus <- length(levels(fsom$metaclustering))
    nClus <- fsom$FlowSOM$map$nNodes

    cluster_values <- as.numeric(names(cvs))[-length(cvs)]
    width <- max(cluster_values)
    chosen <- which(cluster_values == nMetaClus)
    cv_matrix <- do.call(rbind,
                         lapply(cvs[as.character(cluster_values)],
                                function (x) {
                                    c(x,
                                      rep(NA,
                                          (width - length(x))))
                                }))
    cv_matrix <- rbind(cv_matrix,
                       matrix(c(cvs[[as.character(nClus)]],
                                rep(NA, (width - (nClus %% width)))),
                              ncol = width,
                              byrow = TRUE))
    rownames(cv_matrix)[length(cluster_values) + 1] <- "Original\nclustering"
    colnames(cv_matrix) <- NULL

    disp <- apply(cv_matrix, 2, function(x) as.character(round(x, 2)))
    disp[cv_matrix < show_cv] <- ""
    disp[is.na(cv_matrix)] <- ""
    #print(cv_matrix)
    p_cvs <- pheatmap::pheatmap(cv_matrix,
                                cluster_cols = FALSE,
                                cluster_rows = FALSE,
                                na_col = "white",
                                border_color = "white",
                                gaps_row = c(chosen-1, chosen,
                                             length(cluster_values)),
                                breaks = seq(0, max_cv, length.out = 100),
                                display_numbers = disp)
  
   return(p_cvs) }




# generate_cytonorm_model(file_list222,
#                         input$integer1_cyto_cell_number,
#                         input$integer1_cyto,
#                         input$quantiles,
#                         input$integer2_cyto)


generate_cytonorm_model<- function(files_anchor, number_of_cells, number_of_clusters, quantiles, goal){
      model <- CytoNorm::CytoNorm.train(files = files_anchor, #file_list222, 
      outputDir=  ".", 
      labels = as.numeric(gsub("\\D", "",  files_anchor
                            #file_list222
                            )), 
      channels, 
      transformList = NULL, 
      FlowSOM.params = list(nCells = number_of_cells,
        #input$integer1_cyto_cell_number, 
        xdim = 5, ydim = 5, 
        nClus = number_of_clusters,
        #input$integer1_cyto,
         scale = FALSE), 
      normMethod.train = CytoNorm::QuantileNorm.train,
      normParams = list(nQ = as.numeric(
        #input$quantiles
        ), 
      goal = goal, 
      #input$integer2_cyto), 
      seed = 1, 
      verbose = TRUE))
      return(model)
   }