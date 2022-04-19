# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # load dependencies
  source(file = "functions/app_dependencies.R")
  
  # for storage
  values <- reactiveValues()
  
  # at startup
  observeEvent("", {
    
    # load mediation data and store it
    values$dat <- fread("data/plotData.tsv")
    
    # set executable default PM (will build network, but is fairly large)
    values$pm <- 0.25
    
    # filter possible genes/metabolites based on PM filter
    values$possible_genes <- values$dat[quot.dir.vs.raw.metab >= values$pm | quot.dir.vs.raw.gx >= values$pm, unique(gene)]
    values$possible_metabolites <- values$dat[gene %in% values$possible_genes & (quot.dir.vs.raw.metab >= values$pm | quot.dir.vs.raw.gx >= values$pm), unique(metabolite)]
    
    # set default selection to everything
    values$selected_genes <- "All"
    values$selected_metabolites <- "All"
    
    # supply slider with possible genes/metabolites, but keep selection to "All"
    updateSelectInput(session, "selected_genes", choices = c("All", values$possible_genes), selected = "All")
    updateSelectInput(session, "selected_metabolites", choices = c("All", values$possible_metabolites), selected = "All")
    
    values$error_messages <- c(
      "Error! Unknown life-form detected on board!",
      "Error! Radiation leak detected!",
      "Error! Critical core-temperature has been exceeded!",
      "Error! Self-destruction imminent!",
      "Error! Proximity altert! Collision imminent!",
      "Error! Oxygen levels reaching critical levels!",
      "Error! Hull-breach detected on Deck 4!",
      "Error! Intrusion altert!",
      "Error! Manual override required!"
    )
    
  })
  
  # Input Sync ----
  observeEvent({input$slider_num}, {
    if(input$slider_num != input$slider){
      values$pm <- input$slider_num
      updateSliderInput(session, "slider",value = values$pm)
    }
  })
  observeEvent({input$slider}, {
    if(input$slider_num != input$slider){
      values$pm <- input$slider
      updateNumericInput(session, "slider_num",value = values$pm)
    }
    
  })
  
  
  # PM Input Slider ----
  observeEvent({input$slider | input$slider_num}, {
      
    validate(
      need(
        expr = input$slider == input$slider_num,
        message = "Syncing sliders..."
        )
    )
    
    # debug
    # pm <-.6 
    # possible_genes <- dat[quot.dir.vs.raw.metab >= pm | quot.dir.vs.raw.gx >= pm, unique(gene)]
    # possible_metabolites <- dat[gene %in% possible_genes, unique(metabolite)]
    
    # values$pm <- input$slider
    values$possible_genes <- values$dat[quot.dir.vs.raw.metab >= values$pm | quot.dir.vs.raw.gx >= values$pm, unique(gene)]
    values$selected_genes <- values$selected_genes[values$selected_genes %in% c("All", values$possible_genes)]
    
    # if "All" is in selected genes, get all genes and disregard filter
    if("All" %in% values$selected_genes){
      
      values$possible_metabolites <- values$dat[(quot.dir.vs.raw.metab >= values$pm | quot.dir.vs.raw.gx >= values$pm), unique(metabolite)]
      
    } else {
      
      # filter possible metabolites by mediating with selected genes and PM filter
      values$possible_metabolites <- values$dat[gene %in% values$selected_genes & (quot.dir.vs.raw.metab >= values$pm | quot.dir.vs.raw.gx >= values$pm), unique(metabolite)]
      
    }
    
    # Update Gene/Metabolite Input ----
    updateSelectInput(
      session, 
      "selected_genes",
      choices = c("All", values$possible_genes),
      selected = values$selected_genes[values$selected_genes %in% c("All", values$possible_genes)]
    )
    updateSelectInput(
      session, 
      "selected_metabolites",
      choices = c("All", values$possible_metabolites), 
      selected = values$selected_metabolites[values$selected_metabolites %in% c("All",values$possible_metabolites)]
    ) 
    
  })
  
  # Metabolite Input ----
  observeEvent({
    input$selected_metabolites
  }, {
    
    output$error <- renderText(
      validate(
        need(
          expr = uniqueN(input$selected_metabolites) >=1,
          message = "Nothing selected! Reduce PM filter and/or select a metabolite."))
    )
    
    values$selected_metabolites <- input$selected_metabolites
    
    # in addition, filter selection of genes based on selected metabolites
    
    # ----
    if("All" %in% values$selected_metabolites){
      
      values$possible_genes <- values$dat[(quot.dir.vs.raw.metab >= values$pm | quot.dir.vs.raw.gx >= values$pm), unique(gene)]
      
    } else {
      
      # filter possible metabolites by mediating with selected genes and PM filter
      values$possible_genes <- values$dat[metabolite %in% values$selected_metabolites & (quot.dir.vs.raw.metab >= values$pm | quot.dir.vs.raw.gx >= values$pm), unique(gene)]
      
    }
    
    # update metabolite slider with new possible metabolites and filter selection if applies
    updateSelectInput(session, 
                      "selected_genes",
                      choices = c("All", values$possible_genes), 
                      selected = values$selected_genes[values$selected_genes %in% c("All", values$possible_genes)]
    ) 
    
    # ----
    
  })
  
  # Gene Input ----
  observeEvent({
    input$selected_genes
  }, {
    
    output$error <- renderText(
      validate(
        need(
          expr = uniqueN(input$selected_genes) >= 1,
          message = "Nothing selected! Reduce PM filter and/or select a gene."))
    )
    
    # set new selected genes
    values$selected_genes <- input$selected_genes
    
    
    # if "All" is in selected genes, get all genes and disregard filter
    if("All" %in% values$selected_genes){
      
      values$possible_metabolites <- values$dat[(quot.dir.vs.raw.metab >= values$pm | quot.dir.vs.raw.gx >= values$pm), unique(metabolite)]
      
    } else {
      
      # filter possible metabolites by mediating with selected genes and PM filter
      values$possible_metabolites <- values$dat[gene %in% values$selected_genes & (quot.dir.vs.raw.metab >= values$pm | quot.dir.vs.raw.gx >= values$pm), unique(metabolite)]
      
    }
    
    # update metabolite slider with new possible metabolites and filter selection if applies
    updateSelectInput(session, 
                      "selected_metabolites",
                      choices = c("All", values$possible_metabolites), 
                      selected = values$selected_metabolites[values$selected_metabolites %in% c("All", values$possible_metabolites)]
    ) 
    
  })
  
  observeEvent(input$calculate_network, {
    
    error_messages <- isolate(values$error_messages)
    selected_genes <- isolate(input$selected_genes)
    selected_metabolites <- isolate(input$selected_metabolites)
    pm <- isolate(values$pm)
    dat <- isolate(values$dat)
    
    # ignore previous selection and choose all genes/metabolites when "All" is selected (in addition to other features)
    if("All" %in% selected_genes){
      # selected_genes <- "SERPINA13P"
      selected_genes <- dat[quot.dir.vs.raw.metab >= pm | quot.dir.vs.raw.gx >= pm, unique(gene)]
    }
    if("All" %in% selected_metabolites){
      # selected_metabolites <- "C4OH"
      selected_metabolites <- dat[gene %in% selected_genes & (quot.dir.vs.raw.metab >= pm | quot.dir.vs.raw.gx >= pm), unique(metabolite)]
    }
    
    output$error <- renderText(
      validate(
        need(expr = nrow(dat) >=1, message = "Not enough data to build network! Please choose another filter!"),
        need(expr = length(selected_genes) >= 1 | length(selected_metabolites) >= 1, message = "No mediations for selection!")
      )
    )
    
    # validate(
    #   need(expr = nrow(dat) >=1, message = "Not enough data to build network! Please choose another filter!"),
    #   need(expr = !is.null(dat), message = "Error! Please Retry"),
    #   need(expr = length(selected_genes) >= 1 | length(selected_metabolites) >= 1, message = "No mediations for selection!")
    # )
    
    # create network data based on selected genes & metabolites as well as on selected PM
    network_dat <- create_mediation_network_data(
      d2 = dat[quot.dir.vs.raw.metab >= pm | quot.dir.vs.raw.gx >= pm],
      gene.filter = selected_genes,
      show.weak = TRUE,
      save.table = FALSE)
    
    network_dat <- network_dat[metabolite %in% c(selected_metabolites, NA), ]
    
    # debug messages
    # message(paste(dim(network_dat), collapse = ", "))
    # message(paste(dim(dat), collapse = ", "))
    # message(paste(selected_genes, collapse = ", "))
    # message(paste(selected_metabolites, collapse = ", "))
    
    output$network <- renderPlot({
      
      validate(
        need(expr = length(selected_genes) >= 1 | length(selected_metabolites) >= 1, message = "No genes/metabolites selected!"),
        need(expr = nrow(network_dat) >= 1, message = "No mediations for selection! Change filter and/or metabolite/gene selection!"),
        need(expr = nrow(network_dat) <= 600, message = "Network too large for this puny application! Please apply PM filter or select specific genes/metabolites!"),
        need(expr = (uniqueN(network_dat$gene) >=1) & (uniqueN(network_dat$metabolite) >=1), message = sample(error_messages, 1)),
        need(expr = !is.null(network_dat) & !is.null(dat), message = sample(error_messages, 1))
      )
      mediation_network(
        d2 = dat,
        d2f = network_dat
      )
    },
    width = 1024,
    height = 1024
    )
    
  })
  
} # end of server