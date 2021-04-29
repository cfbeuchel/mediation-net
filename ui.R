
ui <- fluidPage(
  titlePanel("Mediation Explorer"),
  sidebarLayout(
    sidebarPanel(width=3,
      helpText("Explore Mediation Results from "),
      tags$hr(),
      sliderInput("slider",
                  label = "Required Proportion mediated of genes:",
                  min = 0,
                  max = 1,
                  value = 0.2,
                  step = 0.05,
                  sep = "."),
      textOutput("error"),
      selectInput(
        inputId = "selected_genes",
        label = "Show Mediations For These Genes", 
        choices = "All", 
        multiple = T, 
        selected = "All"
          ),
      selectInput(
        inputId = "selected_metabolites",
        label = "Show Mediations For These Metabolites", 
        choices = c("All"),
        multiple = T, 
        selected = "All"
          ),
      actionButton(inputId = "calculate_network", 
                   label = "Show Network"),
      tags$hr(),
      downloadButton("download_data", "Download Data", class = "butt")
    ),
    # mainPanel(O visNetworkOutput("network"))
    mainPanel(shiny::plotOutput("network")
              # DT::dataTableOutput("preview")
              )
  )
)

  
  
