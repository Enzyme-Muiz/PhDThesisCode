library("shinyWidgets")
library(shinybusy)
library(dplyr)
library(shinycookie)
library(cookies)
library(tidyverse)
options(shiny.host = "127.0.0.1")
options(shiny.port = 8899)
source("cytofRUVUI.R")
source("BatchAdjustUI.R")
source("cytoNormUI.R")
source("vaeUI.R")
source("BatchVizUI.R")
source("cytonorm_function.R")
source("cytofRUV_preprocessing.R")
source("for_cytofRUV.R")
source("make_fcs_file.R")

ui <- fluidPage(
  title = "Mass_Cytometry_Batch_Effect_Normalization",
  tags$head(tags$style(
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/css/materialize.min.css")
  ),
  tags$link(rel = "icon", type = "image/png", href = "fav.png")),
  
  tags$head(tags$style(HTML("
                          .fade {opacity: 1;
                      transition: opacity .25s ease-in-out;
                      -moz-transition: opacity .25s ease-in-out;
                      }"))),
  add_busy_spinner(spin = "fading-circle"),
  tabsetPanel(
    BatchVizUI(),
    BatchAdjustUI(),
    cytofRUVUI(),
    cytoNormUI(),
    vaeUI()
  ),
  shinyWidgets::setBackgroundColor("lightgreen"),
  tags$script(HTML(
    "document.body.style.backgroundColor = 'skyblue';"
  ))
)


server <- function(input, output, session) {
  options(shiny.maxRequestSize = 800 * 1024^2)
  observe({
    file11 <- input$file1
    file.copy(file11$datapath, file.path("to_be_corrected", "Batch_1_anchorstim.fcs"), overwrite = TRUE)
    file90AAB <- input$file90AA
    file.copy(file90AAB$datapath, "to_be_corrected\\Batch_1_nonanchor_vali.fcs", overwrite = TRUE)
    file91AAB <- input$file91AA
    file.copy(file91AAB$datapath, "to_be_corrected\\Batch_2_nonanchor_vali.fcs", overwrite = TRUE)
    file92 <- input$file92
    file.copy(file92$datapath, "to_be_corrected\\Batch_3_nonanchor_vali.fcs", overwrite = TRUE)
    file93 <- input$file93
    file.copy(file93$datapath, "to_be_corrected\\Batch_4_nonanchor_vali.fcs", overwrite = TRUE)
    file94 <- input$file94
    file.copy(file94$datapath, "to_be_corrected\\Batch_5_nonanchor_vali.fcs", overwrite = TRUE)
    file95 <- input$file95
    file.copy(file95$datapath, "to_be_corrected\\Batch_6_nonanchor_vali.fcs", overwrite = TRUE)
    file900 <- input$file90
    for (ii in c(1:length(file900$name)))
    {
      file.copy(file900$datapath[ii], paste0("to_be_corrected\\Batch_1_1", ii, "nonanchor1.fcs"), overwrite = TRUE)
    }
    file900A <- input$file90A
    for (ii in c(1:length(file900A$name)))
    {
      file.copy(file900A$datapath[ii], paste0("to_be_corrected\\Batch_1_2", ii, "nonanchor1.fcs"), overwrite = TRUE)
    }
    file22 <- input$file2
    file.copy(file22$datapath, "to_be_corrected\\Batch_2_anchorstim.fcs", overwrite = TRUE)
    file910 <- input$file91
    for (ii in c(1:length(file910$name)))
    {
      file.copy(file910$datapath[ii], paste0("to_be_corrected\\Batch_2_1", ii, "nonanchor1.fcs"), overwrite = TRUE)
    }
    file910A <- input$file91A
    for (ii in c(1:length(file910A$name)))
    {
      file.copy(file910A$datapath[ii], paste0("to_be_corrected\\Batch_2_2", ii, "nonanchor1.fcs"), overwrite = TRUE)
    }
    file33 <- input$file3
    file.copy(file33$datapath, "to_be_corrected\\Batch_3_anchorstim.fcs", overwrite = TRUE)
    file920A <- input$file92A
    for (ii in c(1:length(file920A$name)))

    {
      file.copy(file920A$datapath[ii], paste0("to_be_corrected\\Batch_3_1", ii, "nonanchor1.fcs"), overwrite = TRUE)
    }
    file44 <- input$file4
    file.copy(file44$datapath, "to_be_corrected\\Batch_4_anchorstim.fcs", overwrite = TRUE)
    file930A <- input$file93A
    for (ii in c(1:length(file930A$name)))
    {
      file.copy(file930A$datapath[ii], paste0("to_be_corrected\\Batch_4_1", ii, "nonanchor1.fcs"), overwrite = TRUE)
    }
    file55 <- input$file5
    file.copy(file55$datapath, "to_be_corrected\\Batch_5_anchorstim.fcs", overwrite = TRUE)
    file940A <- input$file94A
    for (ii in c(1:length(file940A$name)))
    {
      file.copy(file940A$datapath[ii], paste0("to_be_corrected\\Batch_5_1", ii, "nonanchor1.fcs"), overwrite = TRUE)
    }
    file66 <- input$file6
    file.copy(file66$datapath, "to_be_corrected\\Batch_6_anchorstim.fcs", overwrite = TRUE)
    file950A <- input$file95A
    for (ii in c(1:length(file950A$name)))
    {
      file.copy(file950A$datapath[ii], paste0("to_be_corrected\\Batch_6_1", ii, "nonanchor1.fcs"), overwrite = TRUE)
    }
    onSessionEnded(function() {
      do.call(file.remove, list(list.files("to_be_corrected\\", pattern = "*.fcs", full.names = TRUE)))
    })
    onSessionEnded(function() {
      do.call(file.remove, list(list.files("CytofRUV_output\\", recursive = TRUE, full.names = TRUE)))
    })
    onSessionEnded(function() {
      do.call(file.remove, list(list.files("res_of_cytofRUV_before_confirmed", full.names = TRUE)))
    })
    onSessionEnded(function() {
      do.call(file.remove, list(list.files("fcs_untransformed\\", pattern = "*.fcs", full.names = TRUE)))
    })
    onSessionEnded(function() {
      do.call(file.remove, list(list.files("www\\concat\\", pattern = "*.fcs", full.names = TRUE)))
    })
    onSessionEnded(function() {
      do.call(file.remove, list(list.files("fcs_aligned\\", pattern = "*.fcs", full.names = TRUE)))
    })
    onSessionEnded(function() {
      do.call(file.remove, unlink("fcs_aligned\\DistributionPlots", recursive = TRUE))
    })

    choicelist2 <- append(c("SD"), paste(input$integer1, "p", sep = ""))
    updatePickerInput(session, inputId = "integer2", choices = choicelist2)
    ############### to update the goal in the
    if (length(list.files("fcs_untransformed\\")) != 0) {
      file_list <- list.files("fcs_untransformed\\")
      file_list2 <- stringr::str_subset(file_list, "Batch_")
      do.call(file.remove, list(list.files("fcs_aligned\\", full.names = TRUE)))
      unlink("fcs_aligned\\DistributionPlots", recursive = TRUE)
      file_list22 <- paste0("fcs_untransformed\\", file_list2)
      where_the_pattern_is <- stringr::str_which(file_list22, "anchorstim")
      file_list222 <- file_list22[where_the_pattern_is]
      goal_list <- as.numeric(gsub("\\D", "", file_list222))
      goal_list <- append(goal_list, "mean")

      updatePickerInput(session, inputId = "integer2_cyto", choices = goal_list)
    }
    if (length(list.files("to_be_corrected")) != 0) #
      {
        jj <- file.path("to_be_corrected", list.files("to_be_corrected")[1])
        # if (file.exists(file.path("to_be_corrected", "Batch_1_anchorstim.fcs")))
        # {jj<- file.path("to_be_corrected", "Batch_1_anchorstim.fcs")
        hjk <- flowCore::read.FCS(jj, truncate_max_range = FALSE)
        for (i in attributes(colnames(flowCore::exprs(hjk))))
        {
          jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
        }
        choice <- c()
        for (i in jk)
        {
          jkl <- flowCore::description(hjk)[[i]]
          choice <- append(choice, jkl)
        }
        choicelist <- c()
        j <- 0
        for (i in choice)
        {
          j <- j + 1
          choicelist[[i]] <- j
        }
        updatePickerInput(session, inputId = "variable", choices = choicelist)
        updateCookie(session, "val" = choicelist)
      } else {
      choice <- c()
      cc <- c("all", "ball")
      for (i in cc)
      {
        choice[[i]] <- 2
      }

      updatePickerInput(session, inputId = "variable", choices = choice)
      updateCookie(session, "val2" = choice)
    }
  })

  observe({
    if (file.exists("fcs_untransformed\\Batch_1_anchorstim.fcs")) {
      jm <- "fcs_untransformed\\Batch_1_anchorstim.fcs"

      hjk <- flowCore::read.FCS(jm, truncate_max_range = FALSE)
      for (i in attributes(colnames(flowCore::exprs(hjk))))
      {
        jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
      }
      choice <- c()
      for (i in jk)
      {
        jkl <- flowCore::description(hjk)[[i]]
        choice <- append(choice, jkl)
      }
      choicelist2 <- c()
      j <- 0
      for (i in choice)
      {
        j <- j + 1
        choicelist2[[i]] <- j
      }
      updatePickerInput(session, inputId = "variable7", choices = choicelist2)
    } else {
      choice <- c()
      cc <- c("all", "ball")
      for (i in cc)
      {
        choice[[i]] <- 2
      }

      updatePickerInput(session, inputId = "variable7", choices = choice)
    }
  })

  observeEvent(input$action66_cyto, {
    try({
      if (length(list.files("fcs_untransformed\\")) != 0) {
        file_list <- list.files("fcs_untransformed\\")
        file_list2 <- stringr::str_subset(file_list, "Batch_")
        do.call(file.remove, list(list.files("fcs_aligned\\", full.names = TRUE)))
        unlink("fcs_aligned\\DistributionPlots", recursive = TRUE)
      }

      hjk <- file.path("fcs_untransformed", list.files("fcs_untransformed\\")[1])
      # hjk<- "fcs_untransformed\\Batch_1_anchorstim.fcs"
      hjk <- flowCore::read.FCS(hjk, truncate_max_range = FALSE)
      for (i in attributes(colnames(flowCore::exprs(hjk))))
      {
        jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
      }

      choice2 <- c()
      for (i in jk)
      {
        jkl <- flowCore::description(hjk)[[i]]
        choice2 <- append(choice2, jkl)
      }
      jjk1 <- choice2[as.numeric(input$variable7_cyto)]
      channels <- jjk1
      file_list22 <- paste0("fcs_untransformed\\", file_list2)
      where_the_pattern_is <- stringr::str_which(file_list22, "anchorstim")
      file_list222 <- file_list22[where_the_pattern_is]

      fsom <- CytoNorm::prepareFlowSOM(file_list222,
        channels,
        nCells = input$integer1_cyto_cell_number,
        FlowSOM.params = list(
          xdim = 10,
          ydim = 10,
          nClus = input$integer1_cyto,
          scale = FALSE
        ),
        transformList = NULL,
        seed = 1
      )
      cvss <- CytoNorm::testCV(fsom, cluster_values = c(3, 5, 7, 10, input$integer1_cyto, 15))
      plot3 <- PlotOverviewCV2(fsom, cvss)
      output$cdk1_cyto <- renderPlot({
        plot3
      })
    })
  })

  observeEvent(input$action20_preprocess, {
    try({
      antigen_numbers <- Panel_making_function("fcs_untransformed", "CytofRUV_output")
      Metadata_making_function("fcs_untransformed", "CytofRUV_output")
      list_in_fcs_untransformed <- list.files("fcs_untransformed", full.names = TRUE)
      file.copy(list_in_fcs_untransformed, "CytofRUV_output")



      setwd("CytofRUV_output")
      rep_samples <- c()
      for (i in c(1:6)) {
        anchor_vali <- read_excel("Metadata.xlsx")$sample_id[str_detect(read_excel("Metadata.xlsx")$sample_id, "anchor|vali")]
        results <- anchor_vali[str_detect(anchor_vali, as.character(i))]
        if (!identical(results, character(0))) {
          rep_samples <- append(rep_samples, results[1])
        }
      }
      updatePickerInput(session, inputId = "choose_anchor", choices = rep_samples)
      updateSliderInput(session, "svd_d_number", max = antigen_numbers - 5)


      setwd("..")
    })
  })

  observeEvent(input$action21_Normalize, {
    try({
      source("mmd.R")
      setwd("CytofRUV_output")
      clusters_nb <- as.numeric(input$cluster_cytofRUV)
      seed <- as.numeric(input$seed_numberID)
      rep_samples <- c()
      for (i in c(1:6)) {
        anchor_vali <- read_excel("Metadata.xlsx")$sample_id[str_detect(read_excel("Metadata.xlsx")$sample_id, "anchor|vali")]
        results <- anchor_vali[str_detect(anchor_vali, as.character(i))]
        if (!identical(results, character(0))) {
          rep_samples <- append(rep_samples, results[1])
        }
      }

      setwd("..")
      result <- append(input$choose_anchor, paste0("vali", str_match(input$choose_anchor, "\\d")[, 1]))
      result <- rep_samples
      loading_data(
        clusters_nb = clusters_nb, seed = seed, rep_samples = result,
        k_value = input$svd_d_number, copy_results_here = "res_of_cytofRUV_before_confirmed"
      )
    })
    mmd_result <- mmd_heatmap("res_of_cytofRUV_before_confirmed",
      transform = "no"
    )

    mmd_result <- mmd_result$size
    rownames(mmd_result) <- colnames(mmd_result)
    output$table_outcome_cytofRuv <- renderTable({
      mmd_result
    })
  })

  observeEvent(input$action_confirm_cytoRUV, {
    try({
      setwd("CytofRUV_output")
      clusters_nb <- as.numeric(input$cluster_cytofRUV)
      seed <- as.numeric(input$seed_numberID)
      rep_samples <- c()
      for (i in c(1:6)) {
        anchor_vali <- read_excel("Metadata.xlsx")$sample_id[str_detect(read_excel("Metadata.xlsx")$sample_id, "anchor|vali")]
        results <- anchor_vali[str_detect(anchor_vali, as.character(i))]
        if (!identical(results, character(0))) {
          rep_samples <- append(rep_samples, results[1])
        }
      }

      setwd("..")
      result <- append(input$choose_anchor, paste0("vali", str_match(input$choose_anchor, "\\d")[, 1]))
      result <- rep_samples
      loading_data(clusters_nb = clusters_nb, seed = seed, rep_samples = result, k_value = input$svd_d_number)
    })
  })

  observeEvent(input$action6_cyto, {
    source("mmd.R")
    if (length(list.files("fcs_untransformed\\")) != 0) {
      file_list <- list.files("fcs_untransformed\\")
      file_list2 <- stringr::str_subset(file_list, "Batch_")
      do.call(file.remove, list(list.files("fcs_aligned\\", full.names = TRUE)))
      unlink("fcs_aligned\\DistributionPlots", recursive = TRUE)
    }


    hjk <- "fcs_untransformed\\Batch_1_anchorstim.fcs"
    hjk <- flowCore::read.FCS(hjk, truncate_max_range = FALSE)
    for (i in attributes(colnames(flowCore::exprs(hjk))))
    {
      jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
    }

    choice2 <- c()
    for (i in jk)
    {
      jkl <- flowCore::description(hjk)[[i]]
      choice2 <- append(choice2, jkl)
    }
    jjk1 <- choice2[as.numeric(input$variable7_cyto)]
    channels <- jjk1
    file_list22 <- paste0("fcs_untransformed\\", file_list2)
    where_the_pattern_is <- stringr::str_which(file_list22, "anchorstim")
    file_list222 <- file_list22[where_the_pattern_is]

    where_the_pattern_iss <- stringr::str_which(file_list22, "vali")
    file_list223 <- file_list22[where_the_pattern_iss]
    #### copy to the root
    for (iii in file_list22)


    {
      file.copy(iii, basename(iii), overwrite = TRUE)
    }
    model <- CytoNorm::CytoNorm.train(
      files = file_list222,
      outputDir = ".",
      labels = as.numeric(gsub("\\D", "", file_list222)),
      channels,
      transformList = NULL,
      FlowSOM.params = list(nCells = input$integer1_cyto_cell_number, xdim = 5, ydim = 5, nClus = input$integer1_cyto, scale = FALSE),
      normMethod.train = CytoNorm::QuantileNorm.train,
      normParams = list(
        nQ = as.numeric(input$quantiles),
        goal = input$integer2_cyto
      ),
      seed = 1,
      verbose = TRUE
    )
    CytoNorm::CytoNorm.normalize(
      model = model,
      files = c(basename(file_list223), basename(file_list222)),
      labels = as.numeric(gsub("\\D", "", c(basename(file_list223), basename(file_list222)))),
      transformList = NULL,
      transformList.reverse = NULL,
      normMethod.normalize = CytoNorm::QuantileNorm.normalize,
      outputDir = ".",
      prefix = "Norm_",
      clean = TRUE,
      verbose = TRUE
    )

    list_file <- list.files(path = ".", full.names = TRUE)
    file_list_after <- stringr::str_starts(list_file, "./Norm_")
    file_list_after <- list_file[file_list_after]
    file_list_after_1 <- stringr::str_ends(file_list_after, "anchorstim.fcs")
    file_list_after_2 <- stringr::str_ends(file_list_after, "vali.fcs")
    file_list_after2 <- file_list_after[file_list_after_2]

    file_list_after <- file_list_after[file_list_after_1]
    file_list_before <- stringr::str_starts(list_file, "./Batch_")
    file_list_before <- list_file[file_list_before]
    file_list_before_1 <- stringr::str_ends(file_list_before, "anchorstim.fcs")
    file_list_before_2 <- stringr::str_ends(file_list_before, "vali.fcs")
    file_list_before2 <- file_list_before[file_list_before_2]
    file_list_before <- file_list_before[file_list_before_1]


    for (iii in file_list_after)


    {
      file.copy(iii, paste0("fcs_untransformed\\", basename(iii)), overwrite = TRUE)
      file.remove(iii)
    }
    for (iii in file_list_before)
    {
      file.remove(iii)
    }


    for (iii in file_list_after2)


    {
      file.copy(iii, paste0("fcs_untransformed\\", basename(iii)), overwrite = TRUE)
      file.remove(iii)
    }
    for (iii in file_list_before2)
    {
      file.remove(iii)
    }

    ### removing the before for the non-vali non-anchor
    file_list_before <- stringr::str_starts(list_file, "./Batch_")
    file_list_before <- list_file[file_list_before]
    file_list_before_2 <- stringr::str_ends(file_list_before, "nonanchor1.fcs")
    file_list_before2 <- file_list_before[file_list_before_2]
    for (iii in file_list_before2)
    {
      file.remove(iii)
    }



    display <- mmd_heatmap2("fcs_untransformed/")
    output$cdk15_cyto <- renderTable({
      display
    })


    list_file <- list.files(path = "fcs_untransformed/", , full.names = TRUE)

    file_list_after <- stringr::str_starts(list_file, "fcs_untransformed/Norm_")
    file_list_after <- list_file[file_list_after]


    for (iii in file_list_after)


    {
      file.remove(iii)
    }
  })

  observeEvent(
    input$action15_cyto,
    try(
      {{ if (length(list.files("fcs_untransformed\\")) != 0) {
        file_list <- list.files("fcs_untransformed\\")
        file_list2 <- stringr::str_subset(file_list, "Batch_")
        do.call(file.remove, list(list.files("fcs_aligned\\", full.names = TRUE)))
        unlink("fcs_aligned\\DistributionPlots", recursive = TRUE)
      }


      hjk <- "fcs_untransformed\\Batch_1_anchorstim.fcs"
      hjk <- flowCore::read.FCS(hjk, truncate_max_range = FALSE)
      for (i in attributes(colnames(flowCore::exprs(hjk))))
      {
        jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
      }

      choice2 <- c()
      for (i in jk)
      {
        jkl <- flowCore::description(hjk)[[i]]
        choice2 <- append(choice2, jkl)
      }
      jjk1 <- choice2[as.numeric(input$variable7_cyto)]
      channels <- jjk1
      file_list22 <- paste0("fcs_untransformed\\", file_list2)
      where_the_pattern_is <- stringr::str_which(file_list22, "anchorstim")
      file_list222 <- file_list22[where_the_pattern_is]

      where_the_pattern_iss <- stringr::str_which(file_list22, "nonanchor1")
      file_list223 <- file_list22[where_the_pattern_iss]
      file_list22 <- paste0("fcs_untransformed\\", file_list2)
      vali_files_to_be_corrected <- stringr::str_which(file_list22, "nonanchor_vali")
      vali_files_to_be_corrected <- file_list22[vali_files_to_be_corrected]


      #### copy to the root
      for (iii in file_list22)


      {
        file.copy(iii, basename(iii), overwrite = TRUE)
      }
      model <- CytoNorm::CytoNorm.train(
        files = file_list222,
        outputDir = ".",
        labels = as.numeric(gsub("\\D", "", file_list222)),
        channels,
        transformList = NULL,
        FlowSOM.params = list(nCells = input$integer1_cyto_cell_number, xdim = 5, ydim = 5, nClus = input$integer1_cyto, scale = FALSE),
        normMethod.train = CytoNorm::QuantileNorm.train,
        normParams = list(
          nQ = as.numeric(input$quantiles),
          goal = input$integer2_cyto
        ),
        seed = 1,
        verbose = TRUE
      )

      ### generate labels for the nonanchors non-vali
      want <- grep(pattern = "(^\\D+\\_\\d)", x = basename(file_list223), value = TRUE)
      want <- sub(pattern = "(^\\D+\\_\\d).*", replacement = "\\1", x = want)
      want <- as.numeric(gsub("\\D", "", want))

      CytoNorm::CytoNorm.normalize(
        model = model,
        files = c(basename(file_list223), basename(file_list222), basename(vali_files_to_be_corrected)),
        labels = c(want, as.numeric(gsub("\\D", "", c(basename(file_list222), basename(vali_files_to_be_corrected))))),
        transformList = NULL,
        transformList.reverse = NULL,
        normMethod.normalize = CytoNorm::QuantileNorm.normalize,
        outputDir = ".",
        prefix = "Norm_",
        clean = TRUE,
        verbose = TRUE
      )

      list_file <- list.files(path = ".", full.names = TRUE)

      file_list_after <- stringr::str_starts(list_file, "./Norm_")
      file_list_after <- list_file[file_list_after]
      file_list_after_1 <- stringr::str_ends(file_list_after, "anchorstim.fcs")
      file_list_after_2 <- stringr::str_ends(file_list_after, "nonanchor1.fcs")
      vali_files_corrected <- stringr::str_ends(file_list_after, "vali.fcs")

      file_list_after2 <- file_list_after[file_list_after_2]
      vali_files_corrected <- file_list_after[vali_files_corrected]

      file_list_after <- file_list_after[file_list_after_1]

      file_list_before <- stringr::str_starts(list_file, "./Batch_")
      file_list_before <- list_file[file_list_before]
      file_list_before_1 <- stringr::str_ends(file_list_before, "anchorstim.fcs")
      file_list_before_2 <- stringr::str_ends(file_list_before, "nonanchor1.fcs")
      file_list_before2 <- file_list_before[file_list_before_2]
      file_list_before <- file_list_before[file_list_before_1]

      for (iii in file_list_after)


      {
        file.copy(iii, paste0("fcs_aligned\\", basename(iii)), overwrite = TRUE)
        file.remove(iii)
      }
      for (iii in file_list_before)
      {
        file.remove(iii)
      }


      for (iii in file_list_after2)


      {
        file.copy(iii, paste0("fcs_aligned\\", basename(iii)), overwrite = TRUE)
        file.remove(iii)
      }

      for (iii in vali_files_corrected)


      {
        file.copy(iii, paste0("fcs_aligned\\", basename(iii)), overwrite = TRUE)
        file.remove(iii)
      }

      for (iii in file_list_before2)
      {
        file.remove(iii)
      }

      file_list_before <- stringr::str_starts(list_file, "./Batch_")
      file_list_before <- list_file[file_list_before]
      file_list_before_2 <- stringr::str_ends(file_list_before, "vali.fcs")
      file_list_before2 <- file_list_before[file_list_before_2]
      for (iii in file_list_before2)
      {
        file.remove(iii)
      }


      path <- "fcs_aligned//"

      list_in_fcs_untransformed <- list.files(path, full.names = TRUE)
      list2 <- list_in_fcs_untransformed[str_detect(list_in_fcs_untransformed, ".fcs")]
      file.rename(list2, paste0(path, str_match(basename(list2), "(Batch_.+)")[, 1]))
      list_in_fcs_untransformed <- list.files(path, full.names = TRUE)
      list2 <- list_in_fcs_untransformed[str_detect(list_in_fcs_untransformed, ".fcs")]
      file.copy(list2, "fcs_untransformed", overwrite = TRUE) }}
    )
  )

  observeEvent(input$action6, {
    try({
      if (length(list.files("fcs_untransformed\\")) != 0) {
        file_list <- list.files("fcs_untransformed\\")
        file_list2 <- stringr::str_subset(file_list, "Batch_")
        do.call(file.remove, list(list.files("fcs_aligned\\", full.names = TRUE)))
        unlink("fcs_aligned\\DistributionPlots", recursive = TRUE)
      }
      hjk <- "fcs_untransformed\\Batch_1_anchorstim.fcs"
      hjk <- flowCore::read.FCS(hjk, truncate_max_range = FALSE)
      for (i in attributes(colnames(flowCore::exprs(hjk))))
      {
        jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
      }

      choice2 <- c()
      for (i in jk)
      {
        jkl <- flowCore::description(hjk)[[i]]
        choice2 <- append(choice2, jkl)
      }
      jjk1 <- choice2[as.numeric(input$variable7)]


      sink("channelstoadjust.txt")
      for (channel in jjk1)
      {
        cat(channel)
        cat("\n")
      }
      sink()

      source("BatchAdjust.R")
      BatchAdjust(
        basedir = "fcs_untransformed\\",
        outdir = "fcs_aligned\\",
        channelsFile = "channelstoadjust.txt",
        batchKeyword = "Batch_",
        anchorKeyword = "anchorstim",
        method = input$integer2,
        transformation = FALSE,
        addExt = NULL,
        plotDiagnostics = FALSE
      )
      my_txt <- readLines("channelstoadjust.txt")
      df1 <- rep(my_txt, each = 2)
      df2 <- rep(c("source", "validation"), times = length(my_txt))
      rom22 <- list.files(path = "fcs_untransformed/", pattern = "*.fcs", full.names = TRUE)

      length <- c()
      for (xyz in rom22)
      {
        yz <- flowCore::read.FCS(xyz, truncate_max_range = FALSE)
        kkk <- dim(flowCore::exprs(yz))[1]
        length <- append(length, kkk)
      }
      minimum1 <- min(length)
      minimum2 <- min(minimum1, 2000)
      where_the_pattern_is <- stringr::str_which(rom22, "anchorstim")
      rom_source <- rom22[where_the_pattern_is]
      BB1 <- flowCore::read.FCS(rom_source[1], truncate_max_range = FALSE)
      BB2 <- flowCore::read.FCS(rom_source[2], truncate_max_range = FALSE)
      newBBT1 <- flowCore::exprs(BB1)
      newBBT2 <- flowCore::exprs(BB2)
      newTT1 <- newBBT1[sample(nrow(newBBT1), minimum2), ]
      newTT2 <- newBBT2[sample(nrow(newBBT2), minimum2), ]
      newTT1 <- as.data.frame(newTT1)
      newTT2 <- as.data.frame(newTT2)
      where_the_pattern_iss <- stringr::str_which(rom22, "nonanchor_vali")
      rom_validation <- rom22[where_the_pattern_iss]
      BB3 <- flowCore::read.FCS(rom_validation[1], truncate_max_range = FALSE)
      BB4 <- flowCore::read.FCS(rom_validation[2], truncate_max_range = FALSE)
      newBBT3 <- flowCore::exprs(BB3)
      newBBT4 <- flowCore::exprs(BB4)
      newTT3 <- newBBT3[sample(nrow(newBBT3), minimum2), ]
      newTT4 <- newBBT4[sample(nrow(newBBT4), minimum2), ]
      newTT3 <- as.data.frame(newTT3)
      newTT4 <- as.data.frame(newTT4)
      df3 <- c()
      for (parameters in my_txt)
      {
        AAA1 <- newTT1[[parameters]]
        AAA2 <- newTT2[[parameters]]
        AAA3 <- newTT3[[parameters]]
        AAA4 <- newTT4[[parameters]]
        temporal <- c()
        suppressWarnings(e_source <- dgof::ks.test(AAA1, AAA2))
        attributes(e_source$statistic) <- NULL
        suppressWarnings(e_validation <- dgof::ks.test(AAA3, AAA4))
        attributes(e_validation$statistic) <- NULL
        temporal <- append(temporal, e_source$statistic)
        temporal <- append(temporal, e_validation$statistic)

        df3 <- append(df3, temporal)
      }
      rom22 <- list.files(path = "fcs_aligned/", pattern = "*.fcs", full.names = TRUE)
      where_the_pattern_is <- stringr::str_which(rom22, "anchorstim")
      rom_source <- rom22[where_the_pattern_is]
      BB1 <- flowCore::read.FCS(rom_source[1], truncate_max_range = FALSE)
      BB2 <- flowCore::read.FCS(rom_source[2], truncate_max_range = FALSE)
      newBBT1 <- flowCore::exprs(BB1)
      newBBT2 <- flowCore::exprs(BB2)
      newTT1 <- newBBT1[sample(nrow(newBBT1), minimum2), ]
      newTT2 <- newBBT2[sample(nrow(newBBT2), minimum2), ]
      newTT1 <- as.data.frame(newTT1)
      newTT2 <- as.data.frame(newTT2)
      where_the_pattern_iss <- stringr::str_which(rom22, "nonanchor_vali")
      rom_validation <- rom22[where_the_pattern_iss]
      BB3 <- flowCore::read.FCS(rom_validation[1], truncate_max_range = FALSE)
      BB4 <- flowCore::read.FCS(rom_validation[2], truncate_max_range = FALSE)
      newBBT3 <- flowCore::exprs(BB3)
      newBBT4 <- flowCore::exprs(BB4)
      newTT3 <- newBBT3[sample(nrow(newBBT3), minimum2), ]
      newTT4 <- newBBT4[sample(nrow(newBBT4), minimum2), ]
      newTT3 <- as.data.frame(newTT3)
      newTT4 <- as.data.frame(newTT4)
      df4 <- c()
      for (parameters in my_txt)
      {
        AAA1 <- newTT1[[parameters]]
        AAA2 <- newTT2[[parameters]]
        AAA3 <- newTT3[[parameters]]
        AAA4 <- newTT4[[parameters]]
        temporal <- c()
        suppressWarnings(e_source <- dgof::ks.test(AAA1, AAA2))
        attributes(e_source$statistic) <- NULL
        suppressWarnings(e_validation <- dgof::ks.test(AAA3, AAA4))
        attributes(e_validation$statistic) <- NULL
        temporal <- append(temporal, e_source$statistic)
        temporal <- append(temporal, e_validation$statistic)

        df4 <- append(df4, temporal)
      }
      show_table <- c()

      show_table <- rbind(show_table, df1)
      show_table <- rbind(show_table, df2)
      show_table <- rbind(show_table, df3)
      show_table <- rbind(show_table, df4)
      display_table <- show_table


      output$cdk1 <- renderTable({
        display_table
      })
    })
  })

  observeEvent(input$action3, {
    rom6 <- list.files(path = "fcs_untransformed/", pattern = "*anchorstim.fcs", full.names = TRUE)
    empty <- c()


    length <- c()
    #
    Batches_num <- c()
    for (xyz in rom6)
    {
      yz <- flowCore::read.FCS(xyz, truncate_max_range = FALSE)
      kkk <- dim(flowCore::exprs(yz))[1]
      length <- append(length, kkk)
      #
      Base_name <- basename(xyz)
      #
      batch_name <- strsplit(Base_name, "")[[1]][7]
      #
      Batches_num <- append(Batches_num, batch_name)
    }
    minimum1 <- min(length)


    minimum2 <- min(minimum1, 1000)
    #### to get the batch_number
    for (jj in rom6)
    {
      BB <- flowCore::read.FCS(jj, truncate_max_range = FALSE)
      BBT <- data.frame(flowCore::exprs(BB))
      set.seed(45)
      subset <- BBT[sample(nrow(BBT), minimum2), ]
      empty <- rbind(empty, subset)
    }
    set.seed(4)
    umapframe <- umap::umap(empty)

    colnames(umapframe$layout) <- c("UMAP_1", "UMAP_2")

    umapframe$layout <- data.frame(umapframe$layout)
    #
    # Batches_num= paste(c(1:length(rom6)))

    Batches <- rep(Batches_num, each = minimum2)

    plot1 <- ggplot2::ggplot(data = umapframe$layout, ggplot2::aes(x = UMAP_1, y = UMAP_2, group = Batches, col = Batches)) +
      # ggplot2::scale_color_brewer(palette = "PuOr")+
      ggplot2::geom_point(shape = 18, size = 3) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black")) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold", colour = "black")) +
      ggplot2::theme(axis.text = ggplot2::element_text(face = "bold", colour = "black", size = 12)) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(face = "bold", colour = "black")) +
      ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1))


    output$oid2 <- renderPlot({
      plot1
    })
    output$foo <- downloadHandler(
      filename = "test_name.jpeg",
      content = function(file) {
        ggplot2::ggsave(file, plot = plot1, device = "jpeg", antialias = "none", dpi = 300)
      }
    )
  })

  observeEvent(input$action7, {
    source("mmd.R")
    mmd_result <- mmd_heatmap("fcs_untransformed",
      transform = "no"
    )

    mmd_result <- mmd_result$size
    mmd_result<- as.data.frame(mmd_result)
    mmd_result_for_plot <- mmd_result %>%
    rownames_to_column() %>%
    gather(colname, value, -rowname)
    # plot1<- grid.table(mmd_result)
    mmd_plot<- ggplot(mmd_result_for_plot, aes(x = rowname, y = colname, fill = value)) +
    geom_tile()
    output$oid2 <- renderPlot({
      #grid.table(mmd_result)
      mmd_plot
    })
  })

  observeEvent(input$action15, {
    # if (file.exists())
    try({
      if (length(list.files("fcs_untransformed\\")) != 0) {
        file_list <- list.files("fcs_untransformed\\")
        file_list2 <- stringr::str_subset(file_list, "Batch_")
        # for (iii in file_list2)


        #  {file.copy(iii, paste("fcs_untransformed\\", basename(iii)), overwrite = TRUE)
        # }
        do.call(file.remove, list(list.files("fcs_aligned\\", full.names = TRUE)))
        # do.call(file.remove, unlink("fcs_aligned\\DistributionPlots", recursive = TRUE))})
        unlink("fcs_aligned\\DistributionPlots", recursive = TRUE)

        # file.remove("fcs_aligned")
      }

      # files <- list.files(path="to_be_corrected\\", pattern="*.fcs", full.names=TRUE, recursive=FALSE)

      hjk <- "fcs_untransformed\\Batch_1_anchorstim.fcs"
      # hjk<- flowCore::read.FCS(jj)
      # rom <- list.files(path = "C:/Users/oaona/app_batch_correction/www/csv_untransformed/", pattern= "*.fcs", full.names = TRUE)
      hjk <- flowCore::read.FCS(hjk, truncate_max_range = FALSE)
      for (i in attributes(colnames(flowCore::exprs(hjk))))
      {
        jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
      }

      choice2 <- c()
      for (i in jk)
      {
        jkl <- flowCore::description(hjk)[[i]]
        choice2 <- append(choice2, jkl)
      }
      jjk1 <- choice2[as.numeric(input$variable7)]


      sink("channelstoadjust.txt")
      for (channel in jjk1)
      {
        cat(channel)
        cat("\n")
      }
      sink()

      source("BatchAdjust.R")
      BatchAdjust(
        basedir = "fcs_untransformed\\",
        outdir = "fcs_aligned\\",
        channelsFile = "channelstoadjust.txt",
        batchKeyword = "Batch_",
        anchorKeyword = "anchorstim",
        method = input$integer2,
        transformation = FALSE,
        addExt = NULL,
        plotDiagnostics = FALSE
      )
      my_txt <- readLines("channelstoadjust.txt")
      df1 <- rep(my_txt, each = 2)
      df2 <- rep(c("source", "validation"), times = length(my_txt))
      rom22 <- list.files(path = "fcs_untransformed/", pattern = "*.fcs", full.names = TRUE)

      length <- c()
      for (xyz in rom22)
      {
        yz <- flowCore::read.FCS(xyz, truncate_max_range = FALSE)
        kkk <- dim(flowCore::exprs(yz))[1]
        length <- append(length, kkk)
      }
      minimum1 <- min(length)

      minimum2 <- min(minimum1, 2000)
      where_the_pattern_is <- stringr::str_which(rom22, "anchorstim")
      rom_source <- rom22[where_the_pattern_is]
      BB1 <- flowCore::read.FCS(rom_source[1], truncate_max_range = FALSE)
      BB2 <- flowCore::read.FCS(rom_source[2], truncate_max_range = FALSE)
      newBBT1 <- flowCore::exprs(BB1)
      newBBT2 <- flowCore::exprs(BB2)
      newTT1 <- newBBT1[sample(nrow(newBBT1), minimum2), ]
      newTT2 <- newBBT2[sample(nrow(newBBT2), minimum2), ]
      newTT1 <- as.data.frame(newTT1)
      newTT2 <- as.data.frame(newTT2)
      where_the_pattern_iss <- stringr::str_which(rom22, "nonanchor_vali")
      rom_validation <- rom22[where_the_pattern_iss]
      BB3 <- flowCore::read.FCS(rom_validation[1], truncate_max_range = FALSE)
      BB4 <- flowCore::read.FCS(rom_validation[2], truncate_max_range = FALSE)
      newBBT3 <- flowCore::exprs(BB3)
      newBBT4 <- flowCore::exprs(BB4)
      newTT3 <- newBBT3[sample(nrow(newBBT3), minimum2), ]
      newTT4 <- newBBT4[sample(nrow(newBBT4), minimum2), ]
      newTT3 <- as.data.frame(newTT3)
      newTT4 <- as.data.frame(newTT4)
      df3 <- c()
      for (parameters in my_txt)
      {
        AAA1 <- newTT1[[parameters]]
        AAA2 <- newTT2[[parameters]]
        AAA3 <- newTT3[[parameters]]
        AAA4 <- newTT4[[parameters]]
        temporal <- c()
        suppressWarnings(e_source <- dgof::ks.test(AAA1, AAA2))
        attributes(e_source$statistic) <- NULL
        suppressWarnings(e_validation <- dgof::ks.test(AAA3, AAA4))
        attributes(e_validation$statistic) <- NULL
        temporal <- append(temporal, e_source$statistic)
        temporal <- append(temporal, e_validation$statistic)

        df3 <- append(df3, temporal)
      }
      rom22 <- list.files(path = "fcs_aligned/", pattern = "*.fcs", full.names = TRUE)
      where_the_pattern_is <- stringr::str_which(rom22, "anchorstim")
      rom_source <- rom22[where_the_pattern_is]
      BB1 <- flowCore::read.FCS(rom_source[1], truncate_max_range = FALSE)
      BB2 <- flowCore::read.FCS(rom_source[2], truncate_max_range = FALSE)
      newBBT1 <- flowCore::exprs(BB1)
      newBBT2 <- flowCore::exprs(BB2)
      newTT1 <- newBBT1[sample(nrow(newBBT1), minimum2), ]
      newTT2 <- newBBT2[sample(nrow(newBBT2), minimum2), ]
      newTT1 <- as.data.frame(newTT1)
      newTT2 <- as.data.frame(newTT2)
      where_the_pattern_iss <- stringr::str_which(rom22, "nonanchor_vali")
      rom_validation <- rom22[where_the_pattern_iss]
      BB3 <- flowCore::read.FCS(rom_validation[1], truncate_max_range = FALSE)
      BB4 <- flowCore::read.FCS(rom_validation[2], truncate_max_range = FALSE)
      newBBT3 <- flowCore::exprs(BB3)
      newBBT4 <- flowCore::exprs(BB4)
      newTT3 <- newBBT3[sample(nrow(newBBT3), minimum2), ]
      newTT4 <- newBBT4[sample(nrow(newBBT4), minimum2), ]
      newTT3 <- as.data.frame(newTT3)
      newTT4 <- as.data.frame(newTT4)
      df4 <- c()
      for (parameters in my_txt)
      {
        AAA1 <- newTT1[[parameters]]
        AAA2 <- newTT2[[parameters]]
        AAA3 <- newTT3[[parameters]]
        AAA4 <- newTT4[[parameters]]
        temporal <- c()
        suppressWarnings(e_source <- dgof::ks.test(AAA1, AAA2))
        attributes(e_source$statistic) <- NULL
        suppressWarnings(e_validation <- dgof::ks.test(AAA3, AAA4))
        attributes(e_validation$statistic) <- NULL
        temporal <- append(temporal, e_source$statistic)
        temporal <- append(temporal, e_validation$statistic)

        df4 <- append(df4, temporal)
      }
      show_table <- c()

      show_table <- rbind(show_table, df1)
      show_table <- rbind(show_table, df2)
      show_table <- rbind(show_table, df3)
      show_table <- rbind(show_table, df4)
      display_table <- show_table
      print("hello")
      print(display_table)
      output$cdk1 <- renderTable({
        display_table
      })

      file_list <- list.files("fcs_untransformed\\")
      rom <- list.files(path = "fcs_untransformed/", pattern = "*.fcs", full.names = TRUE)
      length <- c()
      for (xyz in rom)
      {
        yz <- flowCore::read.FCS(xyz, truncate_max_range = FALSE)
        kkk <- dim(flowCore::exprs(yz))[1]
        length <- append(length, kkk)
      }
      kkj <- dim(flowCore::exprs(yz))[2]
      pp <- cumsum(length)
      qq <- pp
      qq <- R.utils::insert(qq, 1, 1)
      qq[1] <- 0


      concatBBT <- c()
      for (jj in rom)
      {
        BB <- flowCore::read.FCS(jj, truncate_max_range = FALSE)
        newBBT <- flowCore::exprs(BB)[, c(1:kkj)]

        concatBBT <- rbind(concatBBT, newBBT)
      }
      number_of_markers <- dim(concatBBT)[2]
      num_of_fcs <- length(qq) - 1
      kkk <- c()
      for (i in c(1:(num_of_fcs - 1)))
      {
        jjj <- c()
        j <- c((i + 1):num_of_fcs)
        for (q in j)
        {
          assign(paste("S", i, q, sep = ""), c())
          jjj <- append(jjj, paste("S", i, q, sep = ""))
        }
        kkk <- append(kkk, jjj)
      }

      distance <- c()
      for (i in c(1:(num_of_fcs - 1)))
      {
        j <- c((i + 1):num_of_fcs)
        for (q in j)
        {
          for (w in c(1:number_of_markers))
          {
            bbc <- qq[i] + 1
            cnn <- qq[i + 1]
            nta <- qq[q] + 1
            fox <- qq[q + 1]

            rommX <- concatBBT[c(bbc:cnn), ]
            rommY <- concatBBT[c(nta:fox), ]
            suppressWarnings(e <- dgof::ks.test(
              rommX[sample(nrow(rommX), minimum2), ][, w],
              rommY[sample(nrow(rommY), minimum2), ][, w]
            ))
            attributes(e$statistic) <- NULL
            assign(paste("S", i, q, sep = ""), append(get(paste("S", i, q, sep = "")), e$statistic))
          }
          distance <- cbind(distance, get(paste("S", i, q, sep = "")))
        }
      }

      file_list <- list.files("fcs_aligned\\")
      file_list2 <- stringr::str_subset(file_list, "Batch_")
      for (iii in file_list2)
      {
        file.copy(paste("fcs_aligned\\", iii, sep = ""), paste("fcs_untransformed\\", basename(iii), sep = ""), overwrite = TRUE)
      }
      do.call(file.remove, list(list.files("fcs_aligned\\", full.names = TRUE)))
      # do.call(file.remove, unlink("fcs_aligned\\DistributionPlots", recursive = TRUE))})
      unlink("fcs_aligned\\DistributionPlots", recursive = TRUE)

      # file.remove("fcs_aligned")
    })
  })

  observeEvent(input$action, {
    files <- list.files(path = "to_be_corrected\\", pattern = "*.fcs", full.names = TRUE, recursive = FALSE)
    #### abc
    # hjk<- file.path("to_be_corrected", "Batch_1_anchorstim.fcs")
    # hjk<- file.path("to_be_corrected", "Batch_1_anchorstim.fcs")
    #
    hjk <- file.path("to_be_corrected", list.files("to_be_corrected")[1])


    # hjk<- flowCore::read.FCS(jj)
    # rom <- list.files(path = "C:/Users/oaona/app_batch_correction/www/csv_untransformed/", pattern= "*.fcs", full.names = TRUE)
    hjk <- flowCore::read.FCS(hjk, truncate_max_range = FALSE)
    for (i in attributes(colnames(flowCore::exprs(hjk))))
    {
      jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
    }

    choice2 <- c()
    for (i in jk)
    {
      jkl <- flowCore::description(hjk)[[i]]
      choice2 <- append(choice2, jkl)
    }
    jjk1 <- choice2[as.numeric(input$variable)]

    ###### obtain the metal description of the batch1_anchor
    translist <- flowCore::transformList(colnames(flowCore::exprs(hjk)), CytoNorm::cytofTransform)
    hjk <- flowCore::transform(hjk, translist)
    jjk <- flowCore::exprs(hjk)[, as.numeric(input$variable)]
    our_head <- colnames(jjk)
    ####
    for (file in files)
    {
      jjk <- flowCore::exprs(flowCore::read.FCS(file, truncate_max_range = FALSE))[, our_head]
      jjk <- jjk[rowSums(jjk[]) > 0, ]
      jjk <- na.omit(jjk)
      jjk <- jjk[!is.infinite(rowSums(jjk)), ]

      colnames(jjk) <- jjk1



      colnames(jjk) <- jjk1
      if (grepl("1_anchor", basename(file), fixed = TRUE)) {
        write.table(jjk, file.path("vae", "target_train_data.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
      } else if (grepl("2_anchor", basename(file), fixed = TRUE)) {
        write.table(jjk, file.path("vae", "source_train_data.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
      } else if (grepl("1_nonanchor_vali", basename(file), fixed = TRUE)) {
        write.table(jjk, file.path("vae", "target_test_data.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
      } else if (grepl("2_nonanchor_vali", basename(file), fixed = TRUE)) {
        write.table(jjk, file.path("vae", "source_test_data.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
      } else {}

      fcs_file_name <- paste(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file)), ".fcs", sep = "")
      make_fcs_from_df(jjk, "fcs_untransformed", fcs_file_name)
    }
  })

  observeEvent(input$action4, {
    files <- list.files(path = file.path("to_be_corrected"), pattern = "*.fcs", full.names = TRUE, recursive = FALSE)

    hjk <- files[1]
    hjk <- flowCore::read.FCS(hjk, truncate_max_range = FALSE)
    for (i in attributes(colnames(flowCore::exprs(hjk))))
    {
      jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
    }
    choice2 <- c()
    for (i in jk)
    {
      jkl <- flowCore::description(hjk)[[i]]
      choice2 <- append(choice2, jkl)
    }
    jjk1 <- choice2[as.numeric(input$variable)]
    ###### obtain the metal description of the batch1_anchor
    translist <- flowCore::transformList(colnames(flowCore::exprs(hjk)), CytoNorm::cytofTransform)
    hjk <- flowCore::transform(hjk, translist)
    jjk <- flowCore::exprs(hjk)[, as.numeric(input$variable)]
    our_head <- colnames(jjk)
    for (file in files)
    {
      BB <- flowCore::read.FCS(file, truncate_max_range = FALSE)
      translist <- flowCore::transformList(colnames(flowCore::exprs(BB)), CytoNorm::cytofTransform)
      BB <- flowCore::transform(BB, translist)
      jjk <- flowCore::exprs(BB)[, our_head]
      jjk <- jjk[rowSums(jjk[]) > 0, ]
      jjk <- na.omit(jjk)
      jjk <- jjk[!is.infinite(rowSums(jjk)), ]
      colnames(jjk) <- jjk1
      if (grepl("1_anchor", basename(file), fixed = TRUE)) {
        write.table(jjk, file.path("vae", "target_train_data.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
      } else if (grepl("2_anchor", basename(file), fixed = TRUE)) {
        write.table(jjk, file.path("vae", "source_train_data.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
      } else if (grepl("1_nonanchor_vali", basename(file), fixed = TRUE)) {
        write.table(jjk, file.path("vae", "target_test_data.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
      } else if (grepl("2_nonanchor_vali", basename(file), fixed = TRUE)) {
        write.table(jjk, file.path("vae", "source_test_data.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
      }
      fcs_file_name <- paste(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file)), ".fcs", sep = "")
      make_fcs_from_df(jjk, "fcs_untransformed", fcs_file_name)
    }
  })

  observeEvent(input$action20, {
    if (length(list.files("fcs_untransformed\\")) != 0) {
      jm <- file.path("fcs_untransformed", list.files("fcs_untransformed\\")[1])
      hjk <- flowCore::read.FCS(jm, truncate_max_range = FALSE)
      for (i in attributes(colnames(flowCore::exprs(hjk))))
      {
        jk <- paste(substr(i, 1, nchar(i) - 1), "S", sep = "")
      }
      choice <- c()
      for (i in jk)
      {
        jkl <- flowCore::description(hjk)[[i]]
        choice <- append(choice, jkl)
      }
      choicelist2 <- c()
      j <- 0
      for (i in choice)
      {
        j <- j + 1
        choicelist2[[i]] <- j
      }
      if (length(list.files("fcs_untransformed\\")) != 0) {
        file_list <- list.files("fcs_untransformed\\")
        file_list2 <- stringr::str_subset(file_list, "Batch_")
        do.call(file.remove, list(list.files("fcs_aligned\\", full.names = TRUE)))
        unlink("fcs_aligned\\DistributionPlots", recursive = TRUE)
        file_list22 <- paste0("fcs_untransformed\\", file_list2)
        where_the_pattern_is <- stringr::str_which(file_list22, "anchorstim")
        file_list222 <- file_list22[where_the_pattern_is]
        goal_list <- as.numeric(gsub("\\D", "", file_list222))
        goal_list <- append(goal_list, "mean")
        updatePickerInput(session, inputId = "integer2_cyto", choices = goal_list)
      }
      updatePickerInput(session, inputId = "variable7", choices = choicelist2)
      updatePickerInput(session, inputId = "variable7_cyto", choices = choicelist2)
      ###################### update my goal here
    }
  })

  observeEvent(input$csv_button, {
    do.call(file.remove, list(list.files("vae\\", full.names = TRUE)))
    # it will work when fcss are already either transformed or untransformed
    if (length(list.files("fcs_untransformed\\")) != 0) { # get all the files in the fcs_untransformed folder
      jm <- file.path("fcs_untransformed", list.files("fcs_untransformed\\"))
      ############################ get the column names to be used
      hjk <- jm[1]
      hjk <- flowCore::read.FCS(hjk, truncate_max_range = FALSE)
      jjk <- flowCore::exprs(hjk)
      our_head <- colnames(jjk)
      for (j in jm)
      {
        jjk <- flowCore::exprs(flowCore::read.FCS(j, truncate_max_range = FALSE))[, our_head]
        colnames(jjk) <- our_head
        jkl <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(j))
        write.table(jjk, file.path("vae", paste0(jkl, ".csv")), row.names = FALSE, col.names = FALSE, sep = ",")
      }
    }
  })

  observeEvent(input$action500, {
    system(paste("python C:/Users/oaona/Downloads/BatchEffectVaev2/calibrate.py --data_path vae --experiment_name c15_beta_delta.05_cytof_mlp",
      # "python C:/Users/oaona/app_batch_correction/python_reticlulate/BatchEffectRemoval2018-master/calibrate.py --data_path C:/Users/oaona/OneDrive/Desktop/JAWAD/uri --experiment_name c15_beta_delta.05_cytof_mlp",
      paste("--n_epochs", input$n_epochs, sep = " "),
      paste("--code_dim", input$code_dim, sep = " "),
      paste("--delta", input$delta, sep = " "),
      paste("--beta", input$beta, sep = " "),
      paste("--gamma", input$gamma, sep = " "),
      paste("--batch_size", input$batch_size, sep = " "),
      sep = " "
    ))

    # system("python C:/Users/oaona/app_batch_correction/python_reticlulate/BatchEffectRemoval2018-master/calibrate.py --n_epochs 20 --delta 0.1 --data_path C:/Users/oaona/OneDrive/Desktop/JAWAD/uri --experiment_name c15_beta.2_gamma10.0_delta.05_cytof_mlp")
  })

  observeEvent(input$action2, {
    rom <- list.files(path = "fcs_untransformed/", pattern = ("*stim.fcs|*nonanchor_vali.fcs"), full.names = TRUE)
    Batches_num <- c()
    for (xyz in rom)
    {
      Base_name <- basename(xyz)
      #
      batch_name <- strsplit(Base_name, "")[[1]][7]
      batch_name <- as.numeric(batch_name)
      #
      Batches_num <- append(Batches_num, batch_name)
    }
    length <- c()
    for (xyz in rom)
    {
      yz <- flowCore::read.FCS(xyz)
      kkk <- dim(flowCore::exprs(yz))[1]
      length <- append(length, kkk)
    }
    kkj <- dim(flowCore::exprs(yz))[2]
    pp <- cumsum(length)
    qq <- pp
    qq <- R.utils::insert(qq, 1, 1)
    qq[1] <- 0

    minimum1 <- min(length)
    concatBBT <- c()
    for (jj in rom)
    {
      BB <- flowCore::read.FCS(jj, truncate_max_range = FALSE)
      newBBT <- flowCore::exprs(BB)[, c(1:kkj)]

      concatBBT <- rbind(concatBBT, newBBT)
    }




    for (i in Batches_num[-length(Batches_num)])
    {
      jjj <- c()
      j <- Batches_num[-1]
      for (q in j)
      {
        assign(paste("Batch_", i, "_vs_", q, sep = ""), c())
        jjj <- append(jjj, paste("Batch_", i, "_vs_", q, sep = ""))
      }




      number_of_markers <- dim(concatBBT)[2]
      num_of_fcs <- length(qq) - 1


      players <- as.character(Batches_num)
      rankings <- c(1:length(rom))
      league <- setNames(as.list(rankings), players)




      kkk <- c()
      for (i in Batches_num[-length(Batches_num)])
      {
        jjj <- c()
        j <- Batches_num[-1]
        for (q in j)
        {
          assign(paste("Batch_", i, "_vs_", q, sep = ""), c())
          jjj <- append(jjj, paste("Batch_", i, "_vs_", q, sep = ""))
        }
        kkk <- append(kkk, jjj)
      }
    }
    minimum2 <- min(c(minimum1, 5000))
    distance <- c()
    for (i in Batches_num[-length(Batches_num)])
    {
      j <- Batches_num[c((league[[as.character(i)]] + 1):num_of_fcs)]
      for (q in j)
      {
        for (w in c(1:number_of_markers))
        {
          bbc <- qq[league[[as.character(i)]]] + 1
          cnn <- qq[league[[as.character(i)]] + 1]
          nta <- qq[league[[as.character(q)]]] + 1
          fox <- qq[league[[as.character(q)]] + 1]
          rommX <- concatBBT[c(bbc:cnn), ]
          rommY <- concatBBT[c(nta:fox), ]
          suppressWarnings(e <- dgof::ks.test(
            rommX[sample(nrow(rommX), minimum2), ][, w],
            rommY[sample(nrow(rommY), minimum2), ][, w]
          ))
          attributes(e$statistic) <- NULL
          assign(paste("Batch_", i, "_vs_", q, sep = ""), append(get(paste("Batch_", i, "_vs_", q, sep = "")), e$statistic))
        }
        distance <- cbind(distance, get(paste("Batch_", i, "_vs_", q, sep = "")))
      }
    }




    jk <- (colnames(concatBBT))
    attributes(jk) <- NULL
    colnames(distance) <- kkk
    rownames(distance) <- jk
    xdata <- data.frame(distance)
    data.table::setDT(xdata, keep.rownames = "ID")[]
    names(xdata)[1] <- "ID"
    melted <- reshape2::melt(xdata, id.vars = "ID")

    plot1 <- ggplot2::ggplot(data = melted, ggplot2::aes(x = value, y = ID, group = variable, col = variable)) +
      ggplot2::geom_point(shape = 21, size = 3, stroke = 2) +
      ggplot2::labs(x = "KS_Distance") +
      # ggplot2::labs(x="KS_Distance", y= "Markers", title="KS_Distance for each marker", fill= "Batches")+

      # ggplot2::guides(fill=ggplot2::guide_legend(title="Batches"))+
      ggplot2::theme(legend.position = "none") +
      # ggplot2::theme(legend.title = ggplot2::element_text(color = "blue", size = 10, face="bold"))+
      # ggplot2::theme(legend.text = ggplot2::element_text(color = "black", size = 10, face="bold"))+
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black")) +
      # ggplot2::theme(plot.title= ggplot2::element_text(hjust=0.5, face="bold"))+
      ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold", colour = "black")) +
      ggplot2::ylab(NULL) +
      ggplot2::theme(axis.text = ggplot2::element_text(face = "bold", colour = "black", size = 12)) +

      # ggplot2::theme(axis.title.y = ggplot2::element_text(face= "bold", colour="black"))+

      ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 2))


    # ggplot2::geom_point(shape= 18, size= 3)+ ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = 'black'))


    output$oid2 <- renderPlot({
      plot1
    })
    output$foo <- downloadHandler(
      filename = "test.jpeg",
      content = function(file) {
        ggplot2::ggsave(file, plot = plot1, device = "jpeg", antialias = "none", dpi = 300)
      }
    )
  })

  observeEvent(input$action5, {
    hjk <- "to_be_corrected\\Batch_1_anchorstim.fcs"
    rom <- list.files(path = "fcs_untransformed/", pattern = ("*stim.fcs|*nonanchor_vali.fcs"), full.names = TRUE)
    length <- c()
    for (xyz in rom)
    {
      yz <- flowCore::read.FCS(xyz, truncate_max_range = FALSE)
      kkk <- dim(flowCore::exprs(yz))[1]
      length <- append(length, kkk)
    }
    kkj <- dim(flowCore::exprs(yz))[2]
    pp <- cumsum(length)
    qq <- pp
    qq <- R.utils::insert(qq, 1, 1)
    qq[1] <- 0



    concatBBT <- c()
    for (jj in rom)
    {
      BB <- flowCore::read.FCS(jj, truncate_max_range = FALSE)
      newBBT <- flowCore::exprs(BB)[, c(1:kkj)]

      concatBBT <- rbind(concatBBT, newBBT)
    }
    make_fcs_from_df(concatBBT, "concat", "concatBBT1.fcs")
    

    set.seed(45)
    fSOM <- FlowSOM::FlowSOM("concat/concatBBT1.fcs",
      # Input options:
      compensate = FALSE, transform = FALSE,
      scale = FALSE,
      # SOM options:
      colsToUse = c(1:length(input$variable)), xdim = 10, ydim = 10,
      # Metaclustering options:
      nClus = 10
    )


    xjx <- 10
    yyy <- fSOM[[2]]
    attributes(yyy) <- NULL
    for (i in c(1:xjx))
    {
      assign(paste("yyy", i, sep = ""), which(i == yyy))
    }
    vv <- c()
    for (x in c(1:xjx))
    {
      yy <- paste("yyy", x, sep = "")
      vv <- append(vv, yy)
    }

    pqr <- c()
    for (i in c(2:(length(rom) + 1)))
    {
      new <- c()
      kjk <- FlowSOM::GetMetaclusters(fSOM)[(qq[i - 1] + 1):qq[i]]
      lengths <- length(kjk)
      kjk <- table(kjk)
      kjk <- as.numeric(kjk)
      new2 <- (kjk / lengths) * 100
      pqr <- rbind(pqr, new2)
    }

    # generate the labels for proportion
    batch_name <- c()
    for (batchsname in rom)
    {
      blm <- substring(basename(batchsname), 1, 7)
      batch_name <- append(batch_name, blm)
    }
    rownames(pqr) <- NULL
    rownames(pqr) <- batch_name


    xdata <- pqr

    xdata <- t(xdata)
    xdata <- data.frame(xdata)
    xdata <- cbind(id <- c(1:10), xdata)
    names(xdata)[1] <- "ID"

    melted <- reshape2::melt(xdata, id.vars = "ID")
    melted$ID <- as.factor(melted$ID)
    melted_pct <- melted %>%
      group_by(ID) %>%
      mutate(pct = prop.table(value))
    plot2 <- ggplot2::ggplot(data = melted_pct, ggplot2::aes(x = ID, y = pct, group = variable, fill = variable)) +
      ggplot2::geom_bar(stat = "identity", color = "black") +
      ggplot2::labs(x = "cluster_number", y = "proportions", title = "Proportion of each cluster") +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Batches")) +
      ggplot2::theme(legend.title = ggplot2::element_text(color = "blue", size = 10)) +
      # ggplot2::theme(panel.background = ggplot2::element_rect(fill= "grey"))+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      # ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%"))
      ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold", colour = "black")) +
      ggplot2::theme(axis.text = ggplot2::element_text(face = "bold", colour = "black", size = 12)) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(face = "bold", colour = "black")) +
      ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1))





    output$oid2 <- renderPlot({
      plot2
    })
    output$foo <- downloadHandler(
      filename = "test.jpeg",
      content = function(file) {
        ggplot2::ggsave(file, plot = plot2, device = "jpeg", antialias = "none", dpi = 300)
      }
    )
  })
}

shinyApp(ui = ui, server = server)
