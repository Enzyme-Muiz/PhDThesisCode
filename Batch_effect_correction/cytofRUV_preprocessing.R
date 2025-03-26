for (i in c(
  "MASS", "ggplot2", "CytoNorm", "reshape2", "flowCore", "cluster", "stringr",
  "Biobase", "flowViz", "Rtsne", "reticulate", "R.utils", "emdist", "FlowSOM", "umap", "dplyr",
  "matrixStats", "combinat", "kernlab", "gridExtra", "writexl"
)) {
  suppressWarnings(suppressMessages(library(i, character.only = TRUE)))
}



Panel_making_function <- function(path = ".", path_to_save = ".") {
  if (length(list.files(path)) != 0 &
    TRUE %in% str_detect(basename(list.files(path, full.names = TRUE)), ".+\\.fcs")
  ) {
    print("abcd")
    hjk <- exprs(flowCore::read.FCS(list.files(path, full.names = TRUE)[str_detect(list.files(path, full.names = TRUE), ".+\\.fcs")][1]))
    translist <- flowCore::transformList(colnames(hjk), cytofTransform)
    ok <- transform(read.FCS(list.files(path, full.names = TRUE)[str_detect(list.files(path, full.names = TRUE), ".+\\.fcs")][1]), translist)
    print("abcde")

    for (i in attributes(colnames(hjk)))
    {
      jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
      print("abcdef")
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
    antigen1 <- str_detect(choice, "(\\d+\\D+)_(.+)")
    choice1 <- choice[antigen1]
    print("abc1")
    antigen <- str_match(choice1, "(\\d+\\D+)_(.+)")[, 3]
    index_of_duplicate <- (duplicated(antigen) | duplicated(antigen, fromLast = TRUE))
    antigen <- antigen[!index_of_duplicate]
    halfway <- str_detect(choice, "(\\d+\\D+)_(.+)")
    print("abc2")
    choice1 <- choice[halfway]
    halfway <- str_match(choice1, "(\\d+\\D+)_(.+)")[, 2]
    print("abc4")
    halfway <- halfway[!index_of_duplicate]
    print("abc3")
    #    fcs_colname<- paste0(str_match(halfway, "(\\d+)(\\D+)")[ ,3], str_match(halfway, "(\\d+)(\\D+)")[ ,2], "Di")
    names <- colnames(exprs(flowCore::read.FCS(list.files(path, full.names = TRUE)[str_detect(list.files(path, full.names = TRUE), ".+\\.fcs")][1])))

    fcs_colname <- names[!names %in% c("Time")]

    print(length(antigen))
    type_number <- length(antigen) %/% 2
    state_number <- length(antigen) - type_number
    panel <- data.frame(
      fcs_colname = fcs_colname,
      antigen = antigen,
      marker_class = c(rep(c("type"), times = type_number), rep(c("state"), times = state_number))
    )
    print("abc7")
    write_xlsx(panel, paste0(path_to_save, "/Panel.xlsx"), col_names = TRUE)
  }
  return(length(antigen))
}







Metadata_making_function <- function(path = ".", path_to_save = ".") {
  if (length(list.files(path)) != 0 &
    TRUE %in% str_detect(basename(list.files(path, full.names = TRUE)), ".+\\.fcs")
  ) {
    others <- list.files(path)[str_detect(list.files(path), "(\\d)_\\d+nonanchor1")]
    anchor <- list.files(path)[str_detect(list.files(path), "anchorstim")]
    vali <- list.files(path)[str_detect(list.files(path), "nonanchor_vali")]
    if (!identical(vali, character(0))) {
      vali_id <- paste0(str_match(vali, "_(\\d+)_nonanchor_(vali)")[, 3], str_match(vali, "_(\\d+)_nonanchor_(vali)")[, 2])
      file_vali <- paste0(vali)
    } else {
      vali_id <- c()
      file_vali <- c()
    }

    if (!identical(others, character(0))) {
      others_id <- tools::file_path_sans_ext(others)
      file_other <- paste0(others)
    } else {
      others_id <- c()
      file_other <- c()
    }
    if (!identical(anchor, character(0))) {
      anchor_id <- paste0(str_match(anchor, "_(\\d+)_(anchor)")[, 3], str_match(anchor, "_(\\d+)_(anchor)")[, 2])
      file_anchor <- paste0(anchor)
    } else {
      anchor_id <- c()
      file_anchor <- c()
    }


    file_name <- c(file_anchor, file_vali, file_other)

    sample_id <- c(anchor_id, vali_id, others_id)
    condition <- c(rep("HC", length(sample_id)))
    anchor_id1 <- str_match(anchor, "_(\\d+)_(anchor)")[, 3]
    vali_id1 <- str_match(vali, "_(\\d+)_nonanchor_(vali)")[, 3]
    patient_id <- c(anchor_id1, vali_id1, others_id)
    anchor_id2 <- str_match(anchor, "_(\\d+)_(anchor)")[, 2]
    vali_id2 <- str_match(vali, "_(\\d+)_nonanchor_(vali)")[, 2]
    batch_other <- str_match(others, "Batch_(\\d)_\\d+")[, 2]
    batch <- c(anchor_id2, vali_id2, batch_other)
    Metadata <- data.frame(
      file_name = file_name,
      sample_id = sample_id,
      condition = condition,
      patient_id = patient_id,
      batch = batch
    )
    write_xlsx(Metadata, paste0(path_to_save, "/Metadata.xlsx"), col_names = TRUE)
  }
}
