cytoNormUI<- function(){
	tabPanel(title= "CytoNorm",
fluidRow (
               column(12,align="center", tags$div(class= "cytonorm_channel",
                                                 tags$em(tags$h1("CytoNorm Normalization")))),), 
sidebarLayout(
sidebarPanel(   
fluidRow            (
               column(4,align="center", tags$div(class= "parameter_cyto",
                                                 tags$em(tags$h4("Parameter_selection_cyto"))))),
               sliderInput("integer1_cyto", "cluster_number:",
                  min = 1, max = 20,
                  value = 5),
sliderInput("integer1_cyto_cell_number", "number_of_cells:",
                  min = 1000, max = 20000,
                  value = 2000, step = 500),
sliderInput("quantiles", "quantile_to_use:",
                  min = 1, max = 101,
                  value = 101, step = 1),
pickerInput("variable7_cyto", label= "channelstoAdjust", choices=c(),  options = list(`actions-box` = TRUE), multiple = TRUE),
pickerInput("integer2_cyto", label= "goal" ,choices=c("mean"),  options = list(`actions-box` = TRUE), multiple = FALSE),
actionButton("action66_cyto", "cluster_number_check"),
actionButton("action6_cyto", "CytoNorm_Adjust"),
actionButton("action15_cyto", "confirm_CytoNorm_Adjust")
  ),
mainPanel(
  plotOutput("cdk1_cyto"),
  tableOutput("cdk15_cyto"),
        )
  )
  )
}