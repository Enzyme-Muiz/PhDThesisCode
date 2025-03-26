for (i in c(
  "MASS", "ggplot2", "CytoNorm", "reshape2", "flowCore", "cluster", "stringr",
  "Biobase", "flowViz", "Rtsne", "reticulate", "R.utils", "emdist", "FlowSOM", "umap", "dplyr",
  "matrixStats", "combinat", "kernlab", "gridExtra"
)) {
  suppressWarnings(suppressMessages(library(i, character.only = TRUE)))
}

## as.numeric(gsub("\\D", "", c(basename(file_list222),basename(vali_files_to_be_corrected))     ))

mmd_heatmap <- function(x, y = NULL, columns_number, transform = "yes") {
  if (!is.null(y)) {
    jfcs_files <- list.files(path = x, full.names = TRUE)
    sample_numbers <- str_which(basename(jfcs_files), y)
    anchor_path <- basename(jfcs_files)[sample_numbers]
    list_file <- paste0(x, anchor_path)
  } else {
    list_file <- list.files(path = x, , full.names = TRUE)
  }

  number_of_file <- length(list_file)
  distance_matrix <- matrix(, nrow = number_of_file, ncol = number_of_file)
  rownames(distance_matrix) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(list_file))
  colnames(distance_matrix) <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(list_file))

  if (number_of_file == 2) {
    BB1 <- flowCore::read.FCS(list_file[1])

    if (transform == "yes") {
      BBT1 <- flowCore::exprs(BB1)
      if ("Time" %in% colnames(BBT1)) {
        BBT1 <- subset(BBT1, select = -c(Time))
      }
      translist <- flowCore::transformList(colnames(BBT1), cytofTransform) #
      BB1 <- transform(BB1, translist) #
    } else {}

    BBT1 <- flowCore::exprs(BB1)[, columns_number]
    if ("Time" %in% colnames(BBT1)) {
      BBT1 <- subset(BBT1, select = -c(Time))
    }
    BB2 <- flowCore::read.FCS(list_file[2])
    if (transform == "yes") {
      BBT2 <- flowCore::exprs(BB2)
      if ("Time" %in% colnames(BBT2)) {
        BBT2 <- subset(BBT2, select = -c(Time))
      }
      translist <- flowCore::transformList(colnames(BBT2), cytofTransform) #
      BB2 <- transform(BB2, translist) #
    } else {}
    BBT2 <- flowCore::exprs(BB2)[, columns_number]
    if ("Time" %in% colnames(BBT2)) {
      BBT2 <- subset(BBT2, select = -c(Time))
    }
    sampling_number <- min(dim(BBT1)[1], dim(BBT2)[1], 2000)
    BBT1 <- BBT1[sample(nrow(BBT1), sampling_number), ]
    BBT1 <- data.matrix(BBT1)
    BBT2 <- BBT2[sample(nrow(BBT2), sampling_number), ]
    BBT2 <- data.matrix(BBT2)
    print(colnames(BBT2))
    print(typeof(BBT2))
    mmdo <- kernlab::kmmd(BBT1, BBT2)
    distance_matrix[1, 2] <- kernlab::mmdstats(mmdo)[1]
    distance_matrix[2, 1] <- kernlab::mmdstats(mmdo)[1]
  } else {
    combination <- combinat::combn(number_of_file, 2)
    for (i in c(1:dim(combination)[2]))
    {
      matrix_position <- combination[, i]
      BB1 <- flowCore::read.FCS(list_file[matrix_position[1]])
      BBT1 <- flowCore::exprs(BB1) # [ , columns_number]
      # BBT1 <- BBT1[,!names(BBT1) %in% c("Time", "time", "Times", "times")]
      # BBT1 <- subset(BBT1, select = -c(Time))
      if ("Time" %in% colnames(BBT1)) {
        BBT1 <- subset(BBT1, select = -c(Time))
      }
      if (transform == "yes") {
        translist <- flowCore::transformList(colnames(BBT1), cytofTransform) #
        BB1 <- transform(BB1, translist) #
      } else {}
      BBT1 <- flowCore::exprs(BB1)[, columns_number] #
      # BBT1 <- subset(BBT1, select = -c(Time))
      if ("Time" %in% colnames(BBT1)) {
        BBT1 <- subset(BBT1, select = -c(Time))
      }

      BB2 <- flowCore::read.FCS(list_file[matrix_position[2]])
      BBT2 <- flowCore::exprs(BB2) # [ , columns_number]
      # BBT2 <- BBT2[,!names(BBT2) %in% c("Time", "time", "Times", "times")]
      # BBT2 <- subset(BBT2, select = -c(Time))
      if ("Time" %in% colnames(BBT2)) {
        BBT2 <- subset(BBT2, select = -c(Time))
      }
      if (transform == "yes") {
        translist <- flowCore::transformList(colnames(BBT2), cytofTransform) #
        BB2 <- transform(BB2, translist) #
      } else {}
      BBT2 <- flowCore::exprs(BB2)[, columns_number] #
      # BBT2 <- BBT2[,!names(BBT2) %in% c("Time", "time", "Times", "times")]
      # BBT2 <- subset(BBT2, select = -c(Time))
      if ("Time" %in% colnames(BBT2)) {
        BBT2 <- subset(BBT2, select = -c(Time))
      }
      sampling_number <- min(dim(BBT1)[1], dim(BBT2)[1], 2000)
      BBT1 <- BBT1[sample(nrow(BBT1), sampling_number), ]
      BBT1 <- data.matrix(BBT1)
      BBT2 <- BBT2[sample(nrow(BBT2), sampling_number), ]
      BBT2 <- data.matrix(BBT2)
      mmdo <- kernlab::kmmd(BBT1, BBT2)
      distance_matrix[matrix_position[1], matrix_position[2]] <- kernlab::mmdstats(mmdo)[1]
      distance_matrix[matrix_position[2], matrix_position[1]] <- kernlab::mmdstats(mmdo)[1]
    }
  }
  distance_matrix[is.na(distance_matrix)] <- 0
  my_list <- list("color" = "red", "size" = distance_matrix, "shape" = "round")
  return(my_list)
}








mmd_heatmap2 <- function(x) {
  list_file <- list.files(path = x, , full.names = TRUE)
  # print(list_file)

  file_list_before <- stringr::str_starts(list_file, "fcs_untransformed/Batch_")
  # print(file_list_before)
  file_list_before <- list_file[file_list_before]
  file_list_before_1 <- stringr::str_ends(file_list_before, "anchorstim.fcs")
  file_list_before <- file_list_before[file_list_before_1]

  file_list_after <- stringr::str_starts(list_file, "fcs_untransformed/Norm_")
  file_list_after <- list_file[file_list_after]
  file_list_after_1 <- stringr::str_ends(file_list_after, "anchorstim.fcs")
  file_list_after <- file_list_after[file_list_after_1]

  # print(paste("file_list_before:",file_list_before))



  file_list_before2 <- stringr::str_starts(list_file, "fcs_untransformed/Batch_")
  # print(file_list_before)
  file_list_before2 <- list_file[file_list_before2]
  file_list_before_12 <- stringr::str_ends(file_list_before2, "vali.fcs")
  file_list_before2 <- file_list_before2[file_list_before_12]

  file_list_after2 <- stringr::str_starts(list_file, "fcs_untransformed/Norm_")
  file_list_after2 <- list_file[file_list_after2]
  file_list_after_12 <- stringr::str_ends(file_list_after2, "vali.fcs")
  file_list_after2 <- file_list_after2[file_list_after_12]




  table_to_display <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("batches", "mmd_distance_before", "mmd_distance_after")
  colnames(table_to_display) <- x

  number_of_file <- length(file_list_before)

  Batches_num2 <- c()
  anchor_mmd_distance_before2 <- c()
  for (i in c(1:(length(file_list_before) - 1)))
  {
    for (j in c((i + 1):length(file_list_before)))
    {
      BB1 <- flowCore::read.FCS(file_list_before[i])
      ii <- as.numeric(gsub("\\D", "", file_list_before[i]))

      BBT1 <- flowCore::exprs(BB1)
      # BBT1 <- subset(BBT1, select = -c(Time))
      if ("Time" %in% colnames(BBT1)) {
        BBT1 <- subset(BBT1, select = -c(Time))
      }


      BB2 <- flowCore::read.FCS(file_list_before[j])
      jj <- as.numeric(gsub("\\D", "", file_list_before[j]))
      ijij <- paste0(ii, jj)
      Batches_num2 <- append(Batches_num2, ijij)
      BBT2 <- flowCore::exprs(BB2)
      # BBT2 <- subset(BBT2, select = -c(Time))
      if ("Time" %in% colnames(BBT2)) {
        BBT2 <- subset(BBT2, select = -c(Time))
      }
      sampling_number <- min(dim(BBT1)[1], dim(BBT2)[1], 2000)
      BBT1 <- BBT1[sample(nrow(BBT1), sampling_number), ]
      BBT1 <- data.matrix(BBT1)
      BBT2 <- BBT2[sample(nrow(BBT2), sampling_number), ]
      BBT2 <- data.matrix(BBT2)
      mmdo_before <- kernlab::kmmd(BBT1, BBT2)

      anchor_mmd_distance_before2 <- append(anchor_mmd_distance_before2, kernlab::mmdstats(mmdo_before)[1])
    }
  }
  anchor_before <- cbind(Batches_num2, anchor_mmd_distance_before2)




  Batches_num2 <- c()
  vali_mmd_distance_before2 <- c()
  for (i in c(1:(length(file_list_before2) - 1)))
  {
    for (j in c((i + 1):length(file_list_before2)))
    {
      BB1 <- flowCore::read.FCS(file_list_before2[i])
      ii <- as.numeric(gsub("\\D", "", file_list_before2[i]))

      BBT1 <- flowCore::exprs(BB1)
      # BBT1 <- subset(BBT1, select = -c(Time))
      if ("Time" %in% colnames(BBT1)) {
        BBT1 <- subset(BBT1, select = -c(Time))
      }

      BB2 <- flowCore::read.FCS(file_list_before2[j])
      jj <- as.numeric(gsub("\\D", "", file_list_before2[j]))
      ijij <- paste0(ii, jj)
      Batches_num2 <- append(Batches_num2, ijij)
      BBT2 <- flowCore::exprs(BB2)
      # BBT2 <- subset(BBT2, select = -c(Time))
      if ("Time" %in% colnames(BBT2)) {
        BBT2 <- subset(BBT2, select = -c(Time))
      }
      sampling_number <- min(dim(BBT1)[1], dim(BBT2)[1], 2000)
      BBT1 <- BBT1[sample(nrow(BBT1), sampling_number), ]
      BBT1 <- data.matrix(BBT1)
      BBT2 <- BBT2[sample(nrow(BBT2), sampling_number), ]
      BBT2 <- data.matrix(BBT2)
      mmdo_before <- kernlab::kmmd(BBT1, BBT2)

      vali_mmd_distance_before2 <- append(vali_mmd_distance_before2, kernlab::mmdstats(mmdo_before)[1])
    }
  }
  vali_before <- cbind(Batches_num2, vali_mmd_distance_before2)








  Batches_num2 <- c()
  anchor_mmd_distance_after2 <- c()
  for (i in c(1:(length(file_list_after) - 1)))
  {
    for (j in c((i + 1):length(file_list_after)))
    {
      BB1 <- flowCore::read.FCS(file_list_after[i])
      ii <- as.numeric(gsub("\\D", "", file_list_after[i]))

      BBT1 <- flowCore::exprs(BB1)
      # BBT1 <- subset(BBT1, select = -c(Time))
      if ("Time" %in% colnames(BBT1)) {
        BBT1 <- subset(BBT1, select = -c(Time))
      }
      BB2 <- flowCore::read.FCS(file_list_after[j])
      jj <- as.numeric(gsub("\\D", "", file_list_after[j]))
      ijij <- paste0(ii, jj)
      Batches_num2 <- append(Batches_num2, ijij)
      BBT2 <- flowCore::exprs(BB2)
      # BBT2 <- subset(BBT2, select = -c(Time))
      if ("Time" %in% colnames(BBT2)) {
        BBT2 <- subset(BBT2, select = -c(Time))
      }

      sampling_number <- min(dim(BBT1)[1], dim(BBT2)[1], 2000)
      BBT1 <- BBT1[sample(nrow(BBT1), sampling_number), ]
      BBT1 <- data.matrix(BBT1)
      BBT2 <- BBT2[sample(nrow(BBT2), sampling_number), ]
      BBT2 <- data.matrix(BBT2)
      mmdo_before <- kernlab::kmmd(BBT1, BBT2)

      anchor_mmd_distance_after2 <- append(anchor_mmd_distance_after2, kernlab::mmdstats(mmdo_before)[1])
    }
  }
  anchor_after <- cbind(Batches_num2, anchor_mmd_distance_after2)

  Batches_num2 <- c()
  vali_mmd_distance_after2 <- c()
  for (i in c(1:(length(file_list_after2) - 1)))
  {
    for (j in c((i + 1):length(file_list_after2)))
    {
      BB1 <- flowCore::read.FCS(file_list_after2[i])
      ii <- as.numeric(gsub("\\D", "", file_list_after2[i]))

      BBT1 <- flowCore::exprs(BB1)
      # BBT1 <- subset(BBT1, select = -c(Time))
      if ("Time" %in% colnames(BBT1)) {
        BBT1 <- subset(BBT1, select = -c(Time))
      }

      BB2 <- flowCore::read.FCS(file_list_after2[j])
      jj <- as.numeric(gsub("\\D", "", file_list_after2[j]))
      ijij <- paste0(ii, jj)
      Batches_num2 <- append(Batches_num2, ijij)
      BBT2 <- flowCore::exprs(BB2)
      # BBT2 <- subset(BBT2, select = -c(Time))
      if ("Time" %in% colnames(BBT2)) {
        BBT2 <- subset(BBT2, select = -c(Time))
      }

      sampling_number <- min(dim(BBT1)[1], dim(BBT2)[1], 2000)
      BBT1 <- BBT1[sample(nrow(BBT1), sampling_number), ]
      BBT1 <- data.matrix(BBT1)
      BBT2 <- BBT2[sample(nrow(BBT2), sampling_number), ]
      BBT2 <- data.matrix(BBT2)
      mmdo_before <- kernlab::kmmd(BBT1, BBT2)

      vali_mmd_distance_after2 <- append(vali_mmd_distance_after2, kernlab::mmdstats(mmdo_before)[1])
    }
  }
  vali_after <- cbind(Batches_num2, vali_mmd_distance_after2)



  merged <- merge(anchor_before, vali_before, by = "Batches_num2", all = TRUE)

  merged <- merge(merged, anchor_after, by = "Batches_num2", all = TRUE)
  merged <- merge(merged, vali_after, by = "Batches_num2", all = TRUE)

  colnames(merged) <- c("batch", "anchor_b", "vali_b", "anchor_a", "vali_a")


  print(merged)




  return(merged)
}
