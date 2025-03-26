cytofRUVUI <- function() {
  tabPanel(
    title = "CytofRUV",
    fluidRow
    (
      column(12, align = "center", tags$div(
        class = "cytofRUV_channel",
        tags$em(tags$h1("CytofRUV Normalization"))
      )),
    ),
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          sliderInput("cluster_cytofRUV", "cluster_number:",
            min = 3, max = 20,
            value = 5
          ),
          numericInput(
            "seed_numberID",
            "seed_number",
            5,
            min = 1,
            max = 2000,
            step = 1,
            width = "200px"
          ),
          sliderInput("svd_d_number", "svd_d_number:",
            min = 3, max = 20,
            value = 5
          ),
          actionButton("action20_preprocess", "Preprocess"),
          pickerInput("choose_anchor",
            label = "anchors_used", choices = c(),
            # options = list(`actions-box` = TRUE),
            # multiple = TRUE
          ),
          actionButton("action21_Normalize", "Normalize", title = "this button will cytofUV-correct batch effect"),
          actionButton("action_confirm_cytoRUV", "confirm_CytofRUV_Adjust")
        )
      ),
      mainPanel(
        tableOutput("table_outcome_cytofRuv"),
      )
    )
  )
}
