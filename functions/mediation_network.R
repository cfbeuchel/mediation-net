mediation_network <- function(d2f,d2) {
  
  
  
  # Create Nodes ----
  
  # get nodes
  d2fm <- melt(d2f, 
               measure.vars = c("metabolite","pheno","gene"))
  d2fm[variable !="metabolite", m.class1:="None"]
  d2fm[, r2.cat.in.bmi := ifelse(variable=="metabolite", r2.cat.metab, r2.cat.gx)]
  d2fm[, beta.reduction := ifelse(variable=="metabolite",quot.dir.vs.raw.metab, quot.dir.vs.raw.gx)]
  d2.nodes <- unique(d2fm[,.(label = value, 
                             group = variable, 
                             r2.cat.in.bmi
  )])
  rm(d2fm)
  d2.nodes[label=="C16"]
  d2.nodes[is.na(r2.cat.in.bmi), r2.cat.in.bmi := "None"]
  d2.nodes[group=="pheno",r2.cat.in.bmi:="BMI"]
  
  # create col for beta-reduction/r2 in node
  d2.nodes <- d2.nodes[,.(label,group,r2.cat.in.bmi)] %>% unique
  
  
  # remove double vertices
  if(any(duplicated(d2.nodes$label))){
    
    # reorder to set the highest always on top
    d2.nodes$r2.cat.in.bmi <- factor(d2.nodes$r2.cat.in.bmi, 
                                     ordered = T, 
                                     levels = c("None", "Low","Middle","High", "Highest", "BMI")[
                                       c("None", "Low","Middle","High", "Highest", "BMI") %in% unique(d2.nodes$r2.cat.in.bmi)
                                     ])
    d2.nodes <- d2.nodes[order(-r2.cat.in.bmi)]
    
    # select only the top (max assoc)
    d2.nodes <- rbindlist(list(
      d2.nodes[!allDuplicatedEntries(label), ],
      d2.nodes[allDuplicatedEntries(label), head(.SD,1), by=label]
    ))
    
    d2.nodes$r2.cat.in.bmi <- as.character(d2.nodes$r2.cat.in.bmi)
    
  }
  
  # some probes with various r2 in bmi are duplicated -> take the mean
  setorder(d2.nodes, -group)
  d2.nodes[,id:=1:.N]
  
  # Create Edges ----
  
  # get edges
  d2.edges <- unique(rbindlist(
    list(
      # pheno~metab
      melt(d2f[mediation.group %in% c(
        "Strongly mediated gene expression effects",
        "Strongly mediated metabolite effects", 
        "Strong bi-directional mediation", 
        "Weak mediation"), 
        .(metabolite,
          pheno,
          gene,
          mediation.group, 
          beta.sign = ifelse(is.na(beta.raw.metab), "None", 
                             ifelse(sign(beta.raw.metab)== 1, "Positive","Negative")),
          beta.strength = cat.beta.raw.metab,
          r2 = r2.raw.metab
        )],
        id.vars = c("metabolite","mediation.group","beta.sign","beta.strength", "r2"),
        measure.vars = c("pheno"))[
          ,.(from = metabolite, 
             to = value, 
             mediation.group = "End", 
             mediation.direction = "None",
             beta.sign, 
             beta.reduction = beta.strength
          )],
      
      # gx~metab
      melt(d2f[,.(mediation.group,
                  mediation.direction,
                  beta.sign = sign.beta.mediation,
                  beta.strength = cat.beta.mediation,
                  beta.reduction
      ),by=.(metabolite,
             pheno,
             gene)],
      id.vars = c("gene","mediation.group","mediation.direction", "beta.sign","beta.strength", "beta.reduction"),
      measure.vars = c("metabolite"))[
        ,.(from = gene, 
           to = value, 
           mediation.group,
           mediation.direction,
           beta.sign,
           # beta.strength,
           beta.reduction
           # r2,
           # r2.quot.dir.vs.raw
           # r2.diff.raw.vs.dir
        )],
      
      # pheno~gene
      melt(d2f[mediation.group %in% c(
        "Strongly mediated gene expression effects",
        "Strongly mediated metabolite effects", 
        "Strong bi-directional mediation", 
        "Weak mediation"),
        .(metabolite,
          pheno,
          gene,
          mediation.group,
          mediation.direction,
          beta.sign = ifelse(is.na(beta.raw.gx), "None", 
                             ifelse(sign(beta.raw.gx)== 1, "Positive","Negative")),
          beta.strength = cat.beta.raw.gx
          # r2 = r2.raw.gx
        )],
        id.vars = c("gene","mediation.group","mediation.direction","beta.sign","beta.strength"),
        measure.vars = c("pheno"))[
          ,.(from = gene, 
             to = value, 
             mediation.group = "End", 
             mediation.direction = "None",
             beta.sign, 
             beta.reduction = beta.strength
          )]
    )
  ))
  
  # create d2 network ----
  # create net object
  d2.edges <- d2.edges[!is.na(from)]
  d2.net <- graph_from_data_frame(d=d2.edges, 
                                  vertices=d2.nodes, 
                                  directed=T) 
  
  # add colors
  colrs <- c("pheno" = "khaki1",
             "gene" = "#2A5676",
             "metabolite" = "#EF4868")
  V(d2.net)$color <- colrs[V(d2.net)$group]
  
  # shape of node
  node.shape <- c("pheno" = "none",
                  "gene" = "circle",
                  "metabolite" = "square")
  V(d2.net)$shape <- node.shape[V(d2.net)$group]
  
  # font of node
  node.font <- c("pheno" = 2,
                 "gene" = 1,
                 "metabolite" = 1)
  V(d2.net)$font <- node.font[V(d2.net)$group]
  
  
  # TODO: edge and node alpha by quot - high quotient, high transparency ----
  # edge.alpha <- c("None" = ,
  #                 "Low" = ,
  #                 "Middle" = ,
  #                 "High" = 
  #                 )
  # d2.edges$beta.strength %>% table
  # d2.edges$beta.reduction %>% table
  
  # edge width based on beta reduction
  edge.width <- c("Low" =0.8,
                  "Middle"=1.2,
                  "High"=1.6,
                  "Highest" = 2,
                  "None"=0,
                  "End" = 0.5) + 1.1 # "neutrale" verbindung zur basis
  E(d2.net)$width <- edge.width[as.character(d2.edges$beta.reduction)] #  + 0.6 # TAG TOPGENES
  
  # arrow direction
  arrow.dir <- c("Metab" = 1,
                 "Gx" = 2,
                 "Both" = 3,
                 "None" = 0)
  E(d2.net)$arrow.mode <- arrow.dir[E(d2.net)$mediation.direction]
  
  # edge color
  edge.col <- c("Strongly mediated metabolite effects" = "#EF4868",
                "Bi-directional mediation" = "palegreen3",
                "Strongly mediated gene expression effects" = "#2A5676",
                "Weak mediation" = "grey35",
                "End" = alpha("grey70",.5))
  E(d2.net)$edge.color <- edge.col[E(d2.net)$mediation.group]
  
  # edge line type - effect sign
  linetypes <- c("None" = "solid", 
                 "Positive" = "longdash",
                 "Both" = "dotdash",
                 "Negative" = "dotted")
  E(d2.net)$edge.lty <- linetypes[E(d2.net)$beta.sign]
  
  
  # size by total r2 in pheno ----
  node.size <- c(
    "None" = 0.5,
    "Low"=1,
    "Middle"=2,
    "High"=4,
    "Highest"=6,
    "BMI"=15
  )
  V(d2.net)$size <- node.size[d2.nodes$r2.cat.in.bmi]
  
  # node font size
  node.cex<- c(
    "gene"=.8,
    "metabolite"=1,
    "pheno"=1.1
  )
  V(d2.net)$vertex.label.cex <- node.cex[V(d2.net)$group]
  
  # node label dist
  node.dist<- c(
    "gene"=1,
    "metabolite"=-1,
    "pheno"=0
  )
  V(d2.net)$node.label.dist <- node.dist[V(d2.net)$group]
  
  # Labels
  if(any(d2.nodes$label=="log.bmi")) {
    d2.nodes[label == "log.bmi",label := "Log-BMI"]
  } else if(any(d2.nodes$label=="diabetes.status.tri")){
    d2.nodes[label == "diabetes.status.tri",label := "T2D"]
  }
  
  # vertex label color based on group
  label_col <- c(
    "gene" = "#2A5676",
    "metabolite" = "#EF4868",
    "pheno" = "grey10"
  )
  V(d2.net)$vertex.label.color <- label_col[V(d2.net)$group]
  
  # network layout
  # https://kateto.net/networks-r-igraph
  # my.layout <- layout_with_kk(d2.net) # good seperation
  # my.layout <- layout_with_lgl(d2.net,root = c("log.bmi")) # chooses root node
  my.layout <- layout_nicely(d2.net)
  
  # legend parameters
  
  lgd.cex <- 1.2
  
  # PLOT ----
  
  plot.igraph(d2.net,
              
              edge.arrow.size=.9,
              edge.curved = F,
              edge.color= E(d2.net)$edge.color,
              edge.lty=E(d2.net)$edge.lty,
              
              
              vertex.shape=V(d2.net)$shape,
              frame.color = V(d2.net)$vertex.label.color,
              vertex.color = adjustcolor(V(d2.net)$color, alpha.f = 0.85),
              
              vertex.label.dist = V(d2.net)$node.label.dist,
              vertex.label=d2.nodes$label,
              vertex.label.color=V(d2.net)$vertex.label.color,
              vertex.label.family="Helvetica",
              vertex.label.font=V(d2.net)$font,
              vertex.label.cex = V(d2.net)$vertex.label.cex +0.7, # + 0.7 # TAG TOPGENES
              
              layout = my.layout
  )
  
  # Legend: Node
  legend(
    x="bottom", #-1,5
    c("Gene","Metabolite"), 
    pch=c(21, 22),
    col="grey30", 
    pt.bg=colrs[-1], 
    pt.cex=2, 
    cex=lgd.cex, 
    bty="n", 
    ncol=1)
  
  # Legend: Edge Color
  legend(x="bottomleft",
         
         c(
           "Strongly mediated metabolite effects",
           "Strong bi-directional mediation",
           "Strongly mediated gene expression effects", 
           "Weak mediation", 
           "No mediation tested"),
         lty="solid",
         col=edge.col,
         lwd=1.5,
         cex=lgd.cex, 
         bty="n", 
         ncol=1)
  
  # Legend: Edge line type
  legend(x="bottomright", 
         c(
           "No sign. effect estimate",
           "Positive standardized effect estimate",
           "Both, positive and negative\n effect estimate",
           "Negative standardized effect estimate"
         ), 
         lty=linetypes,
         col="grey35",
         lwd=1.5,
         cex=lgd.cex,
         bty="n", 
         ncol=1)
  
  # Legend: Edge line thickness
  # beta.red <- c(d2f[sig.mediation.metab==TRUE, 1- quot.dir.vs.raw.metab],d2f[sig.mediation.gx==TRUE, 1- quot.dir.vs.raw.gx])
  # quantile(beta.red, c(0,0.25,0.5,0.75,1)) %>% round(3)
  legend(x="topright",
         title = c("Edges between metabolites and genes:\nReduction in effect estimate by mediation"),
         c("0%–25%",
           "25%–50%",
           "50–75%",
           "75%–100%"),
         lwd=c("Low" =0.8,
               "Middle"=1.2,
               "High"=1.6,
               "Highest"=2) +0.9, # + 0.6 # TAG TOPGENES
         bty="n",
         ncol=1,
         cex=lgd.cex
  )
  
  node <- with(d2, quantile(c(r2.raw.gx, r2.raw.metab), c(0,0.25,0.5,0.75,1),na.rm=T)) %>% round(4)
  
  # Legend: Node thickness
  legend(x="topleft",
         title = c("Raw explained variance of node in BMI"),
         c(
           "None",
           paste(node[1], node[2], sep = "–"),
           paste(node[2], node[3], sep = "–"),
           paste(node[3], node[4], sep = "–"),
           paste(node[4], node[5], sep = "–")
         ),
         pch = c(1),
         pt.cex = c("None" = 0.5,
                    "Low" =1,
                    "Middle"=2,
                    "High"=4,
                    "Highest"=6),
         bty="n",
         ncol=1,
         cex=lgd.cex
  )
  
}