vaeUI<- function(){
tabPanel(title= "VAE",
  fluidRow
             (
              column(12,align="center", tags$div(class= "vae",
                                                 tags$em(tags$h1("variational_autoencoder correction")))),
               
             ), 

sidebarLayout(
sidebarPanel(   

fluidRow
             (
               column(4,align="center", tags$div(class= "parameters",
                                                 tags$em(tags$h4("Parameter_selection"))))),
               


sliderInput("n_epochs", "n_epochs:",
                  min = 1, max = 100,
                  value = 5),
sliderInput("code_dim", "code_dim:",
                  min = 5, max = 20,
                  value = 5),
sliderInput("delta", "delta:",
                  min = 0, max = 1,
                  step = 0.05, 
                  value = 0.1),
sliderInput("beta", "beta:",
                  min = 0, max = 2,
                  step = 0.2, 
                  value = 1),
sliderInput("gamma", "gamma:",
                  min = 0, max = 20,
                  step = 2, 
                  value = 5),
sliderInput("batch_size", "batch_size:",
                  min = 50, max = 500,
                  step = 50, 
                  value = 50),


 actionButton("action500", "Batch_Adjust"),

 actionButton("action55", "confirm_Adjust")
 
     

  ),


mainPanel(
  tableOutput("vae_result"),
        )
  )



  )
}