
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      width=3,
      h4("mediation-net"),
      tags$hr(),
      h5("Explore Mediation Results from $FuturePaperID"),
      tags$hr(),
      numericInput(inputId = "slider_num", 
                   label = "Required Proportion Mediated (PM):",
                   value = 0.25, 
                   min = 0, 
                   step = 0.01,
                   max = 1),
      sliderInput(inputId = "slider",
                  label = NULL,
                  min = 0,
                  max = 1,
                  value = 0.25,
                  step = 0.05,
                  sep = "."),
      textOutput("error"),
      tags$hr(),
      h5("Show Mediations For These:"),
      selectInput(
        inputId = "selected_genes",
        label = "Selected Genes:", 
        choices = "All", 
        multiple = T, 
        selected = "All"
          ),
      selectInput(
        inputId = "selected_metabolites",
        label = "Selected Metabolites:", 
        choices = c("All"),
        multiple = T, 
        selected = "All"
          ),
      actionButton(inputId = "calculate_network", 
                   label = "Show Network"),
      tags$hr()
    ),
    mainPanel(shiny::plotOutput("network"))
  )
)

  
  
