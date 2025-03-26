BatchAdjustUI<- function(){
	tabPanel(title= "Batch_Adjust",

fluidRow
             (
               column(12,align="center", tags$div(class= "per_channel",
                                                 tags$em(tags$h1("Per_channel Normalization")))),
               
             ), 

sidebarLayout(
sidebarPanel(   

fluidRow
             (
               column(4,align="center", tags$div(class= "parameter",
                                                 tags$em(tags$h4("Parameter_selection"))))),
               


sliderInput("integer1", "Percentile:",
                  min = 1, max = 100,
                  value = 5),
pickerInput("variable7", label= "channelstoAdjust", choices=c(),  options = list(`actions-box` = TRUE), multiple = TRUE),

 pickerInput("integer2", label= "parameters" ,choices=c("SD"),  options = list(`actions-box` = TRUE), multiple = FALSE),

 actionButton("action6", "Batch_Adjust"),

 actionButton("action15", "confirm_Adjust")
 
     

  ),


mainPanel(
  tableOutput("cdk1"),
        )
  )
  )
}