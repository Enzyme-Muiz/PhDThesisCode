for (i in c(
  "MASS", "ggplot2", "CytoNorm", "reshape2", "flowCore", "cluster", "stringr",
  "Biobase", "flowViz", "Rtsne", "reticulate", "R.utils", "emdist", "FlowSOM", "umap", "dplyr",
  "matrixStats", "tidyr", "docstring", "scales", "RColorBrewer",
  "plotly", "ggfortify", "ggpubr"
)) {
  suppressWarnings(suppressMessages(library(i, character.only = TRUE)))
}

# given an fcs file
# given transform of yes or not
# given columns_number to use
# return subset of the dataframe
remove_time_and_subset_by_given_column_numbers <- function(BB1, transform, columns_needed, rename_columns_needed) {
  internalbb <- flowCore::exprs(BB1)
  if (transform == "yes") {
    translist <- flowCore::transformList(colnames(internalbb), cytofTransform) #
    BB1 <- transform(BB1, translist) #
  } else {}
  if (!is.null(columns_needed)) {
    for (i in attributes(colnames(internalbb)))
    {
      jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
    }
    choice <- c()
    for (i in jk) {
      jkl <- flowCore::description(BB1)[[i]]
      if (!is.null(jkl)) {
        choice <- append(choice, jkl)
      } else {
        choice <- append(choice, "nothing_here")
      }
    }
    print(choice)
    positions <- c()
    for (i in columns_needed) {
      new_mark <- paste0(i, "$")
      position <- grep(new_mark, choice)
      positions <- append(positions, position)
    }
    print(positions)

    BBT1 <- flowCore::exprs(BB1)[, positions]
    # print(colnames(BBT1))
  } else {
    BBT1 <- flowCore::exprs(BB1)
  }
 
  if ("Time" %in% colnames(BBT1)) {
    BBT1 <- subset(BBT1, select = -c(Time))
  }

  if (!is.null(rename_columns_needed)){
   colnames(BBT1) <- rename_columns_needed
    }
  return(BBT1)
}

boxplot_groupby_input_marker_of_fcs <- function(pathx, substringy, transformTrueorFalse, columns_needed, rename_columns_needed, legend_use, legend_title, main_title) {
  if (!is.null(substringy)) {
    jfcs_files <- list.files(path = pathx, full.names = TRUE)
    sample_numbers <- c()
    for (j in substringy) {
      sample_numbers1 <- str_which(basename(jfcs_files), j)
      sample_numbers <- append(sample_numbers, sample_numbers1)
    }
    anchor_path <- basename(jfcs_files)[sample_numbers]
    list_file <- paste0(pathx, anchor_path)
  } else {
    list_file <- list.files(path = pathx, , full.names = TRUE)
  }
  random_names <- paste0("a", c(1:length(list_file)))
  for (i in c(1:length(list_file))) {
    BB1 <- flowCore::read.FCS(list_file[i])
    text_to_evaluate <- paste0(
      random_names[i], " =",
      "remove_time_and_subset_by_given_column_numbers(BB1, transformTrueorFalse, columns_needed, rename_columns_needed)"
    )
    eval(parse(text = text_to_evaluate))
  }
  result <- c()
  for (i in random_names) {
    text_to_evaluate <- paste0(
      "expression_data", " = as.matrix(",
      i, ")"
    )
    eval(parse(text = text_to_evaluate))
    expression_data <- as.matrix(expression_data)
    res <- colMedians(expression_data)
    result <- rbind(result, res)
  }
  text_to_evaluate <- paste0(
    "column_names", " = colnames(",
    random_names[1], ")"
  )
  eval(parse(text = text_to_evaluate))
  colnames(result) <- column_names

  if (legend_use == TRUE) {
    n <- length(list_file)
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
    colourss <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    # colourss<- c("red", "green", "brown", "black", "pink")
    legend <- c()
    for (i in c(1:length(list_file))) {
      legend_value <- readline(prompt = paste0("Enter legend value for ", list_file[i], ":"))
      legend <- append(legend, as.character(legend_value))
    }
    rownames(result) <- legend
    check <- melt(result)
    colnames(check) <- c("mutation_status", "Var2", "median_intensity")

    p <- ggplot(check, aes(x = Var2, y = median_intensity, fill = mutation_status)) +
      geom_boxplot() +
      ggtitle(main_title) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
      ) +
      scale_fill_manual(
        breaks = unique(legend), #### editable
        values = colourss[1:length(unique(legend))]
      ) +
      scale_fill_discrete(name = legend_title) +
      stat_compare_means(aes(label = ..p.signif..), label.y = 4, hide.ns = TRUE, size = 10)
  } else {
    legend <- c(1:length(list_file))

    rownames(result) <- legend
    check <- melt(result)
    colnames(check) <- c("mutation_status", "Var2", "median_intensity")

    p <- ggplot(check, aes(x = Var2, y = median_intensity)) +
      geom_boxplot(fill = "#1d70b8", color = "black") +
      ggtitle(main_title) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16, face = "bold"),
        # panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
      )
  }
  return(p)
}

colmedian_input_marker_of_fcs <- function(pathx, substringy, transformTrueorFalse, columns_needed, rename_columns_needed){
  if (!is.null(substringy)) {
    jfcs_files <- list.files(path = pathx, full.names = TRUE)
    sample_numbers <- c()
    for (j in substringy) {
      sample_numbers1 <- str_which(basename(jfcs_files), j)
      sample_numbers <- append(sample_numbers, sample_numbers1)
    }
    anchor_path <- basename(jfcs_files)[sample_numbers]
    list_file <- paste0(pathx, anchor_path)
  } else {
    list_file <- list.files(path = pathx, , full.names = TRUE)
  }
  random_names <- paste0("a", c(1:length(list_file)))
  for (i in c(1:length(list_file))) {
    BB1 <- flowCore::read.FCS(list_file[i])
    text_to_evaluate <- paste0(
      random_names[i], " =",
      "remove_time_and_subset_by_given_column_numbers(BB1, transformTrueorFalse, columns_needed, rename_columns_needed)"
    )
    eval(parse(text = text_to_evaluate))
  }
  result <- c()
  for (i in random_names) {
    text_to_evaluate <- paste0(
      "expression_data", " = as.matrix(",
      i, ")"
    )
    eval(parse(text = text_to_evaluate))
    expression_data <- as.matrix(expression_data)
    res <- colMedians(expression_data)
    result <- rbind(result, res)
  }
  text_to_evaluate <- paste0(
    "column_names", " = colnames(",
    random_names[1], ")"
  )
  eval(parse(text = text_to_evaluate))
  colnames(result) <- column_names
  legend <- basename(list_file)

  rownames(result) <- legend

  return(result)
}

pca_plot_from_proportion <- function(df, use_exist = TRUE, legend_title) {
  abc <- df
  if (exists("legendabcde") == TRUE & use_exist == TRUE) {
    legendabcde <<- legendabcde
  } else {
    legendabcde <<- c()
    for (i in rownames(abc)) {
      legend_value <- readline(prompt = paste0("Enter legend value for ", i, ":"))
      while (legend_value == "") {
        legend_value <- readline(prompt = paste0("Enter legend value for ", i, ":"))
      }
      legendabcde <<- append(legendabcde, legend_value)
    }
  }
  legend_title1 <- legend_title
  text_to_evaluate <- paste0(legend_title, " <- ", "legendabcde")
  eval(parse(text = text_to_evaluate))
  # legend_title<- legendabcde
  # print(get(legend_title))
  abc <- cbind(get(legend_title), abc)
  abc <- data.frame(abc)
  colnames(abc)[1] <- legend_title1

  df <- abc[-c(1)]
  df <- data.frame(apply(df, 2, function(x) as.numeric(as.character(x))))
  pca_res <- prcomp(df, scale. = TRUE)
  p <- autoplot(pca_res, data = abc, colour = legend_title1, size = 6)
  p +
    ggplot2::theme(axis.text = ggplot2::element_text(face = "bold", colour = "black", size = 12)) +
    ggplot2::theme(axis.title = ggplot2::element_text(face = "bold", colour = "black", size = 20)) +
    ggplot2::theme(
      legend.title = element_text(face = "bold", colour = "black", size = 12), legend.text = element_text(face = "bold", colour = "black", size = 12),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

fcs_to_dataframe <- function(path_to_fcs_file, list_of_antibody, transform_or_not = TRUE) {
  #' given antibody
  #' given an fcs file
  #' given transform or not
  #' output an R dataframe
  dataframe_result <- exprs(read.FCS(path_to_fcs_file))
  if (transform_or_not == TRUE) {
    translist <- flowCore::transformList(colnames(dataframe_result), cytofTransform)
    ok <- transform(read.FCS(path_to_fcs_file), translist)
  } else {
    ok <- read.FCS(path_to_fcs_file)
  }
  for (i in attributes(colnames(dataframe_result)))
  {
    jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
  }
  choice <- c()
  for (i in jk) {
    jkl <- flowCore::description(ok)[[i]]
    choice <- append(choice, jkl)
    choicelist2 <- c()
    j <- 0
    for (i in choice) {
      j <- j + 1
      choicelist2[[i]] <- j
    }
  }

  columns_number <- match(list_of_antibody, choice)
  dataframe_to_return <- dataframe_result[, columns_number]
  colnames(dataframe_to_return) <- list_of_antibody
  dataframe_to_return <- as.data.frame(dataframe_to_return)
  return(dataframe_to_return)
}

output_path <- function(sample_number, path_to_fcs_files, batch_needed_number = "3rd") {
  # ' given the 4 digit sample number needed
  # ' given the path of the fcs files
  # ' given the batch needed
  # ' output the path of the fcs file
  jfcs_files <- list.files(path = path_to_fcs_files, full.names = TRUE)
  batch_needed <- basename(jfcs_files)[str_starts(basename(jfcs_files), batch_needed_number)]
  sample_numbers <- str_extract(batch_needed, "[0-9]{4}")
  anchor_path <- batch_needed[which(sample_numbers %in% sample_number)]
  full_path_anchor <- paste0(path_to_fcs_files, anchor_path)
  return(full_path_anchor)
}

histogram_of_columns <- function(dataframe, title, xlab, bin) {
  #' given dataframe
  #' given title name
  #' given xlab value
  #' give bin size
  #' output histogram of all columns
  p <- ggplot(gather(dataframe), aes(value)) +
    geom_histogram(bins = bin) +
    facet_wrap(~key, scales = "free_x") +
    xlab(xlab) +
    theme(text = element_text(size = 20)) +
    ggtitle(title)
  return(p)
}

mmd_heatmap <- function(x, y = NULL, columns_number, transform = "yes") {
 if (!is.null(y)) {
    jfcs_files <- list.files(path = x, full.names = TRUE)
    sample_numbers <- c()
    for (j in y) {
      sample_numbers1 <- str_which(basename(jfcs_files), j)
      sample_numbers <- append(sample_numbers, sample_numbers1)
    }
    anchor_path <- basename(jfcs_files)[sample_numbers]
    list_file <- paste0(x, anchor_path)
  } else {
    list_file <- list.files(path = x, , full.names = TRUE)
  }
  # print(list_file)
  number_of_file <- length(list_file)
  distance_matrix <- matrix(, nrow = number_of_file, ncol = number_of_file)
  rownames(distance_matrix) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(list_file))
  colnames(distance_matrix) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(list_file))

  if (number_of_file == 2) {
    BB1 <- flowCore::read.FCS(list_file[1])
    BBT1 <- remove_time_and_subset_by_given_column_numbers(BB1, transform, columns_number)
    BB2 <- flowCore::read.FCS(list_file[2])
    BBT2 <- remove_time_and_subset_by_given_column_numbers(BB2, transform, columns_number)

    sampling_number <- min(dim(BBT1)[1], dim(BBT2)[1], 2000)
    BBT1 <- BBT1[sample(nrow(BBT1), sampling_number), ]
    BBT1 <- data.matrix(BBT1)
    BBT2 <- BBT2[sample(nrow(BBT2), sampling_number), ]
    BBT2 <- data.matrix(BBT2)
    # print(colnames(BBT2))
    # print(typeof(BBT2))
    mmdo <- invisible(kernlab::kmmd(BBT1, BBT2))
    distance_matrix[1, 2] <- kernlab::mmdstats(mmdo)[1]
    distance_matrix[2, 1] <- kernlab::mmdstats(mmdo)[1]
  } else {
    combination <- combinat::combn(number_of_file, 2)
    for (i in c(1:dim(combination)[2]))
    {
      matrix_position <- combination[, i]
      BB1 <- flowCore::read.FCS(list_file[matrix_position[1]])
      BBT1 <- remove_time_and_subset_by_given_column_numbers(BB1, transform, columns_number)

      BB2 <- flowCore::read.FCS(list_file[matrix_position[2]])
      BBT2 <- remove_time_and_subset_by_given_column_numbers(BB2, transform, columns_number)

      sampling_number <- min(dim(BBT1)[1], dim(BBT2)[1], 2000)
      BBT1 <- BBT1[sample(nrow(BBT1), sampling_number), ]
      BBT1 <- data.matrix(BBT1)
      BBT2 <- BBT2[sample(nrow(BBT2), sampling_number), ]
      BBT2 <- data.matrix(BBT2)
      mmdo <- invisible(kernlab::kmmd(BBT1, BBT2))
      distance_matrix[matrix_position[1], matrix_position[2]] <- kernlab::mmdstats(mmdo)[1]
      distance_matrix[matrix_position[2], matrix_position[1]] <- kernlab::mmdstats(mmdo)[1]
    }
  }
  distance_matrix[is.na(distance_matrix)] <- 0
  return_list <- list("dist" = distance_matrix, "files" = list_file)
  return(return_list)
}

produce_heatmap <- function(distance_matrix, list_file, legend_used = FALSE, legend_name = "Legend", use_exist = TRUE) {
  if (legend_used == TRUE) {
    if (exists("legendabcde") == TRUE & use_exist == TRUE) {
      legendabcde <<- legendabcde
    } else {
      legendabcde <<- c()
      for (i in c(1:length(list_file))) {
        legend_value <- readline(prompt = paste0("Enter legend value for ", list_file[i], ":"))
        while (legend_value == "") {
          legend_value <- readline(prompt = paste0("Enter legend value for ", list_file[i], ":"))
        }
        legendabcde <<- append(legendabcde, legend_value)
      }
    }
    my_sample_col <- data.frame(legendabcde)
    colnames(my_sample_col) <- legend_name
    row.names(my_sample_col) <- colnames(distance_matrix)
    pheatmap(distance_matrix, annotation_col = my_sample_col)
  } else {
    pheatmap(distance_matrix)
  }
}

generate_concat <- function(y, x, breadth, column_needed, rename_columns_needed, transform = "yes", name_to_save_concat = "concatuntransformed", fcs_or_csv = "fcs") {
  if (!is.null(y)) {
    jfcs_files <- list.files(path = x, full.names = TRUE)
    sample_numbers <- c()
    for (j in y) {
      sample_numbers1 <- str_which(basename(jfcs_files), j)
      sample_numbers <- append(sample_numbers, sample_numbers1)
    }
    anchor_path <- basename(jfcs_files)[sample_numbers]
    list_file <- paste0(x, anchor_path)
  } else {
    list_file <- list.files(path = x, , full.names = TRUE)
  }

  concatBBT <- c()
  cell_lengths <- c()
  original_indexes_of_sampled <- c()
  for (xyz in list_file)
  {

    BB1 <- read.FCS(xyz)
    BBT <- remove_time_and_subset_by_given_column_numbers(BB1, transform, column_needed, rename_columns_needed)

    if (!is.null(breadth)) {
      sample_indexes <- sample(nrow(BBT), breadth, replace = TRUE)
      BBT <- BBT[sample_indexes, ]
      original_indexes_of_sampled <- append(original_indexes_of_sampled, sample_indexes)
    } else {
      original_indexes_of_sampled <- append(original_indexes_of_sampled, c(1:dim(BBT)[1]))
    }
    cell_length <- dim(BBT)[1]
    concatBBT <- rbind(concatBBT, BBT)
    cell_lengths <- append(cell_lengths, cell_length)
  }

  meta <- data.frame(
    name = colnames(concatBBT),
    desc = paste(colnames(concatBBT))
  )
  if (fcs_or_csv == "fcs"){
      meta$range <- apply(apply(concatBBT, 2, range), 2, diff)
  meta$minRange <- apply(concatBBT, 2, min)
  meta$maxRange <- apply(concatBBT, 2, max)
  ff <- new("flowFrame", exprs = data.matrix(concatBBT), parameters = AnnotatedDataFrame(meta))
  write.FCS(ff, paste0(name_to_save_concat, ".fcs"))
  }
  else if (fcs_or_csv== "csv") {
     write.csv(concatBBT, paste0(name_to_save_concat, ".csv"))
  }


  return_result <- list(
    "cell_lengths" = cell_lengths,
    "columns_length" = dim(concatBBT)[2],
    "file_names" = basename(list_file),
    "original_indexes" = original_indexes_of_sampled
  )
  return(return_result)
}

cluster_proportion_from_concat <- function(name_of_concat_fcs_file, transform1, each_file_length, metaclusters, columns_length, row_names) {
  set.seed(45)
  fSOM <- FlowSOM(name_of_concat_fcs_file,
    compensate = FALSE, transform = transform1,
    toTransform = c(1:columns_length),
    scale = FALSE,
    colsToUse = c(1:columns_length),
    xdim = 10,
    ydim = 10,
    nClus = metaclusters
  )


  pp <- cumsum(each_file_length)
  qq <- pp
  qq <- R.utils::insert(qq, 1, 1)
  qq[1] <- 0
  pqr <- c()
  for (i in c(2:(length(each_file_length) + 1)))
  {
    kjk <- GetMetaclusters(fSOM)[(qq[i - 1] + 1):qq[i]]
    lengths <- length(kjk)
    kjk <- table(kjk)
    check <- sum(kjk)
    # ok_check<- append(ok_check, check)
    kjk <- as.numeric(kjk)
    new2 <- (kjk / lengths) * 100
    pqr <- rbind(pqr, new2)
  }
  rownames(pqr) <- row_names
  return_result <- list(
    "proportion" = pqr,
    "fsom" = fSOM
  )
  return(return_result)
}

barchart_from_proportion <- function(pqr, legend_title, main_title, use_exist, statistics) {
  if (exists("legendabcde") == TRUE & use_exist == TRUE) {
    legendabcde <<- legendabcde
  } else {
    legendabcde <<- c()
    for (i in rownames(pqr)) {
      legend_value <- readline(prompt = paste0("Enter legend value for ", i, ":"))
      while (legend_value == "") {
        legend_value <- readline(prompt = paste0("Enter legend value for ", i, ":"))
      }
      legendabcde <<- append(legendabcde, legend_value)
    }
  }
  col_names <- sprintf("cluster%s", seq(1:dim(pqr)[2]))
  colnames(pqr) <- col_names
  pqr <- cbind(legendabcde, pqr)
  pqr <- data.frame(pqr)

  res <- melt(pqr, id.vars = "legendabcde")
  res$legendabcde <- as.factor(res$legendabcde)
  res$value <- as.numeric(res$value)
  rownames(pqr) <- NULL

  if (statistics == FALSE) {
    p <- ggplot(res, aes(x = variable, y = value, fill = legendabcde)) +
      geom_boxplot(outlier.size = 3, lwd = 1) +
      ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 2)) +
      ggplot2::xlab("Clusters") +
      ggtitle(main_title) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black")) +
      ggplot2::theme(axis.text = ggplot2::element_text(face = "bold", colour = "black", size = 12)) +
      ggplot2::theme(axis.title = ggplot2::element_text(face = "bold", colour = "black", size = 20)) +
      ggplot2::theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
      ggplot2::ylab("Proportion") +
      scale_fill_discrete(name = legend_title) +
      ggplot2::theme(
        legend.title = element_text(face = "bold", colour = "black", size = 12), legend.text = element_text(face = "bold", colour = "black", size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  } else {
    colnames(res) <- c(legend_title, "variable", "value")
    p <- ggpubr::ggboxplot(res,
      x = "variable", y = "value",
      fill = legend_title, palette = "jco"
    )
    p <- p + ggpubr::stat_compare_means(aes(group = get(legend_title)), label = "p.signif", size = 15, hide.ns = TRUE, label.y = aggregate(value ~ variable, data = res, max)$value - 3) +
      ylab("Percentage (%)") +
      xlab(" ") +
      ggplot2::theme(
        legend.title = element_text(face = "bold", colour = "black", size = 12), legend.text = element_text(face = "bold", colour = "black", size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      ggtitle(main_title) +
      ggplot2::theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title.y = element_text(size = 18, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2)
      )
  }
  p
  # return(res)
}

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
