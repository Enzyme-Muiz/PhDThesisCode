BatchVizUI <- function() {
  tabPanel(
    title = "Batch_effect_viz",
    fluidRow
    (
      column(12, align = "center", tags$div(
        class = "fade",
        tags$em(tags$h1("Batch Effect Visualization"))
      )),
    ),
    sidebarLayout
    (
      sidebarPanel(
        width = 7,
        fluidRow(
          column(6,
            center = "align",
            fileInput(
              inputId = "file1",
              label = NULL,
              multiple = FALSE,
              accept = ".FCS",
              width = "250px",
              buttonLabel = "Browse",
              placeholder = "anchor_sam_1"
            )
          ),
          column(6,
            align = "center",
            fileInput(
              inputId = "file2",
              label = NULL,
              multiple = FALSE,
              accept = ".FCS",
              width = "250px",
              buttonLabel = "Browse",
              placeholder = "anchor_sam_2"
            )
          )
        ),
        fluidRow(
          column(6,
            center = "align",
            fileInput(
              inputId = "file3",
              label = NULL,
              multiple = FALSE,
              accept = ".FCS",
              width = "250px",
              buttonLabel = "Browse",
              placeholder = "anchor_sam_3"
            )
          ),
          column(6,
            align = "center",
            fileInput(
              inputId = "file4",
              label = NULL,
              multiple = FALSE,
              accept = ".FCS",
              width = "250px",
              buttonLabel = "Browse",
              placeholder = "anchor_sam_4"
            )
          )
        ),
        fluidRow(
          column(6,
            center = "align",
            fileInput(
              inputId = "file5",
              label = NULL,
              multiple = FALSE,
              accept = ".FCS",
              width = "250px",
              buttonLabel = "Browse",
              placeholder = "anchor_sam_5"
            )
          ),
          column(6,
            align = "center",
            fileInput(
              inputId = "file6",
              label = NULL,
              multiple = FALSE,
              accept = ".FCS",
              width = "250px",
              buttonLabel = "Browse",
              placeholder = "anchor_sam_6"
            )
          )
        )
      ),
      mainPanel(
        width = 4,
        align = "center", plotOutput("oid2", height = "300px")
      )
    ),
    fluidRow(
      column(
        2,
        pickerInput("variable", choices = c(), options = list(`actions-box` = TRUE), multiple = TRUE)
      ),
      column(
        2,
        dropdown(
          actionButton("action", "untransformed"),
          actionButton("action4", "arcsinh_trans"),
          # actionButton("action_untansform_", "arcsinh_trans"),
          label = "Preprocess",
          inputId = "dataset_types"
        )
        # actionButton(inputId = "action",label = "column_selection")
      ),
      column(
        2,
        dropdown(
          actionButton("action20", "update"),
          actionButton("csv_button", "csv_button"),
          actionButton("action2", "KS_batch"),
          actionButton("action3", "UMAP_scatter"),
          actionButton("action5", "clusters_proportion"),
          actionButton("action7", "MMD"),
          # actionButton("action_marker_hist", "histogram"),
          label = "Batch_viz",
          inputId = "Batch_viz2"
        )

        # actionButton(inputId = "action2",label = "KS_batch")
      ),
      column(
        2,
        downloadButton("foo")
      ),
      column(
        2,
        dropdown(
          label = "FILES_INPUT",
          inputId = "files_input",
          fileInput(
            inputId = "file90AA",
            label = NULL,
            multiple = FALSE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "Validation_1"
          ),
          fileInput(
            inputId = "file90",
            label = NULL,
            multiple = TRUE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "batch_1"
          ),
          fileInput(
            inputId = "file90A",
            label = NULL,
            multiple = TRUE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "batch_1_2"
          ),
          fileInput(
            inputId = "file91AA",
            label = NULL,
            multiple = FALSE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "Validation_2"
          ),
          fileInput(
            inputId = "file91",
            label = NULL,
            multiple = TRUE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "batch_2"
          ),
          fileInput(
            inputId = "file91A",
            label = NULL,
            multiple = TRUE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "batch_2_2"
          ),
          fileInput(
            inputId = "file92",
            label = NULL,
            multiple = FALSE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "validation_3"
          ),
          fileInput(
            inputId = "file92A",
            label = NULL,
            multiple = TRUE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "batch_3_2"
          ),
          fileInput(
            inputId = "file93",
            label = NULL,
            multiple = FALSE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "validation_4"
          ),
          fileInput(
            inputId = "file93A",
            label = NULL,
            multiple = TRUE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "batch_4_2"
          ),
          fileInput(
            inputId = "file94",
            label = NULL,
            multiple = FALSE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "validation_5"
          ),
          fileInput(
            inputId = "file94A",
            label = NULL,
            multiple = TRUE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "batch_5_2"
          ),
          fileInput(
            inputId = "file95",
            label = NULL,
            multiple = FALSE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "validation_6"
          ),
          fileInput(
            inputId = "file95A",
            label = NULL,
            multiple = TRUE,
            accept = ".FCS",
            width = "200px",
            buttonLabel = "Browse",
            placeholder = "batch_6_2"
          )
        )
      )
    ),
    tableOutput("contents")
  )
}
