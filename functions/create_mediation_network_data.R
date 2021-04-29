create_mediation_network_data <- function(
  d2, 
  gene.filter, 
  show.weak = FALSE,
  save.table = FALSE) {
  
  # summary of obesity gene assocs
  d2[,.N]
  d2[,.(uniqueN(gene), uniqueN(metabolite))]
  d2[,.(unique(m.class1),unique(m.class2)),by=metabolite][,table(V1)]
  d2[,.(unique(m.class1),unique(m.class2)),by=metabolite][,table(V2)]
  
  # Create Additional Columns ----
  
  # get sign of mediation effect
  # d2$sign.beta.mediation <-  NULL
  d2$sign.beta.mediation <-  character()
  d2[mediation.group=="Strong bi-directional mediation", 
     sign.beta.mediation :=   ifelse(sign(beta.mediation.gx) == sign(beta.mediation.metab),
                                     sign(beta.mediation.gx),
                                     "Both"),]
  d2[mediation.group=="Strongly mediated metabolite effects", sign.beta.mediation := sign(beta.mediation.metab)]
  d2[mediation.group=="Strongly mediated gene expression effects", sign.beta.mediation := sign(beta.mediation.gx)]
  d2[sign.beta.mediation == "-1", sign.beta.mediation := "Negative"]
  d2[sign.beta.mediation == "1", sign.beta.mediation := "Positive"]
  
  # get beta direction of "Weak mediation" group
  d2[mediation.group=="Weak mediation" & ((sig.mediation.gx) == F | is.na(beta.mediation.gx)), sign.beta.mediation := sign(beta.mediation.metab)]
  d2[mediation.group=="Weak mediation" & ((sig.mediation.metab) == F | is.na(beta.mediation.metab)), sign.beta.mediation := sign(beta.mediation.gx)]
  d2[mediation.group=="Weak mediation" & !is.na(beta.mediation.metab) & !is.na(beta.mediation.gx) & 
       sign(beta.mediation.gx) == sign(beta.mediation.metab), # are signs equal?
     sign.beta.mediation := sign(beta.mediation.gx), # they are equal, just choose any
     ]
  d2[mediation.group == "Weak mediation" & sign.beta.mediation == "-1", sign.beta.mediation := "Negative"]
  d2[mediation.group == "Weak mediation" & sign.beta.mediation == "1", sign.beta.mediation := "Positive"]
  d2[mediation.group=="Weak mediation" & !is.na(beta.mediation.metab) & !is.na(beta.mediation.gx) & 
       sign(beta.mediation.gx) != sign(beta.mediation.metab), # are signs unequal?
     sign.beta.mediation := "Both" # they are different, show "both"
     ]
  
  # check
  # d2[sign.beta.mediation == "-1", ]
  
  # check
  # d2[mediation.group=="Weak mediation",.(sign.beta.mediation, 
  #                               beta.mediation.metab, beta.mediation.gx, 
  #                               sig.mediation.gx, sig.mediation.metab)] %>% 
  #   as.data.frame() %>% head(20)
  
  # get size of (strongest) mediation effect beta
  d2[, strength.beta.mediation := ifelse(
    !is.na(beta.mediation.metab) & !is.na(beta.mediation.gx),
    ifelse(abs(beta.mediation.gx) > abs(beta.mediation.metab),
           abs(beta.mediation.gx),
           abs(beta.mediation.metab)),
    ifelse(is.na(beta.mediation.metab),
           abs(beta.mediation.gx), 
           abs(beta.mediation.metab) 
    ))]
  
  # d2$strength.beta.mediation %>% is.na %>% table
  
  # check
  # # d2[,.(mediation.group,
  #       sign.beta.mediation, 
  #       beta.mediation.metab, 
  #       beta.mediation.gx,
  #       strength.beta.mediation)] %>% as.data.frame() %>% head(15)
  
  # mediation direction indication
  d2[, mediation.direction := ifelse(
    sig.mediation.gx == T & sig.mediation.metab == T,
    "Both",
    ifelse(sig.mediation.gx == F & sig.mediation.metab == F,
           "None", ifelse(sig.mediation.gx == T, "Gx", 
                          ifelse(sig.mediation.metab == T,
                                 "Metab", "ERROR")
           )))]
  
  # d2[, .(sig.mediation.gx, sig.mediation.metab, mediation.direction)] %>% as.data.frame %>% head(15)
  # d2$mediation.direction %>% table
  
  # check
  # d2[,.(mediation.group,
  #       sign.beta.mediation, 
  #       beta.mediation.metab, 
  #       beta.mediation.gx,
  #       strength.beta.mediation,
  #       mediation.direction)] %>% as.data.frame() %>% head(15)
  
  
  # get terciles of effect strength and categorize beta strength
  d2$cat.beta.mediation <- with(d2, cut(strength.beta.mediation,
                                        breaks=quantile(strength.beta.mediation, c(0,1/3,2/3,1)),
                                        include.lowest = T,
                                        labels=c("Low","Middle","High")
  ))
  d2[,range(strength.beta.mediation),by=cat.beta.mediation]
  
  # categories for gx
  # try to reduce avaiable categories by using unique
  d2[sig.raw.gx==TRUE, cat.beta.raw.gx:=cut(abs(beta.raw.gx),
                                            breaks=unique(quantile(abs(beta.raw.gx), c(0,1/3,2/3,1))),
                                            include.lowest = T,
                                            labels=c("Low","Middle","High")
  )]
  d2[,range(beta.raw.gx),by=cat.beta.raw.gx]
  d2[sig.raw.gx==FALSE,cat.beta.raw.gx := "None"]
  
  # categories for metab
  if(uniqueN(d2$metabolite==1)){
    
    d2[sig.raw.metab==TRUE,cat.beta.raw.metab:="Singular"]
    
  } else{
    
  
  d2[sig.raw.metab==TRUE, cat.beta.raw.metab:=cut(abs(beta.raw.metab),
                                                  breaks=quantile(abs(beta.raw.metab), c(0,1/3,2/3,1)),
                                                  include.lowest = T,
                                                  labels=c("Low","Middle","High")
  )]
    
  }
  
  d2[,range(beta.raw.metab),by=cat.beta.raw.metab]
  d2[sig.raw.metab==FALSE,cat.beta.raw.metab := "None"]
  
  # categories for beta reduction of gx
  d2[sig.mediation.gx==TRUE, cat.beta.reduction.gx:=cut(
    quot.dir.vs.raw.gx,
    breaks=c(0, 0.25,0.5,0.75, 1),
    # breaks=c(min(quot.dir.vs.raw.gx), 0.25,0.5,0.75, max(quot.dir.vs.raw.gx)),            
    # breaks=quantile(1-quot.dir.vs.raw.gx, c(0,0.25,0.5,0.75,1)),
    include.lowest = T,
    labels=c("Low","Middle","High", "Highest")
  )]
  d2$cat.beta.reduction.gx <- as.character(d2$cat.beta.reduction.gx)
  d2[is.na(cat.beta.reduction.gx), cat.beta.reduction.gx := "None"]
  d2[,.(cat.beta.reduction.gx,sig.mediation.gx)] %>% as.data.frame %>% head(15)
  
  # categories for beta reduction of metab
  d2[sig.mediation.metab==TRUE, cat.beta.reduction.metab:=cut(
    quot.dir.vs.raw.metab,
    breaks=c(0, 0.25,0.5,0.75, 1),
    # breaks=c(min(quot.dir.vs.raw.metab), 0.25,0.5,0.75, max(quot.dir.vs.raw.metab)),
    # breaks=quantile(1-quot.dir.vs.raw.metab, c(0,0.25,0.5,0.75,1)),
    include.lowest = T,
    labels=c("Low","Middle","High", "Highest")
  )]
  d2$cat.beta.reduction.metab <- as.character(d2$cat.beta.reduction.metab)
  d2[is.na(cat.beta.reduction.metab), cat.beta.reduction.metab := "None"]
  d2[,.(quot.dir.vs.raw.metab,cat.beta.reduction.metab,sig.mediation.metab)] %>% as.data.frame %>% head(15)
  
  # fill all the NAs to some "None" Category
  d2[is.na(cat.beta.reduction.metab), cat.beta.reduction.metab := "None"]
  d2[is.na(cat.beta.reduction.gx), cat.beta.reduction.gx := "None"]
  d2[is.na(cat.beta.raw.gx), cat.beta.raw.gx := "None"]
  d2[is.na(cat.beta.raw.metab), cat.beta.raw.metab := "None"]
  
  # beta reduction group
  d2[, beta.reduction := 
       ifelse(
         mediation.group == "Strongly mediated metabolite effects",
         cat.beta.reduction.metab, ifelse(
           mediation.group == "Strongly mediated gene expression effects",
           cat.beta.reduction.gx, ifelse(
             mediation.group == "Strong bi-directional mediation" & quot.dir.vs.raw.gx < quot.dir.vs.raw.metab,
             cat.beta.reduction.gx, 
             ifelse(mediation.group == "Weak mediation" & mediation.direction == "Gx", cat.beta.reduction.gx,
                    ifelse(mediation.group == "Weak mediation" & mediation.direction == "Metab",cat.beta.reduction.metab,
                           ifelse(mediation.group == "Weak mediation" & mediation.direction == "Both" & (cat.beta.reduction.gx == cat.beta.reduction.metab),
                                  cat.beta.reduction.gx, "CHECK!"))))))]
  
  d2[,.(cat.beta.reduction.gx,
        cat.beta.reduction.metab, 
        beta.reduction,
        mediation.group,
        mediation.direction)] %>% as.data.frame() %>% head(15)
  
  #check
  stopifnot(d2[beta.reduction == "CHECK!", .N] == 0 )
  
  # r2 in bmi category metab/gx
  
  if(any(duplicated(quantile(c(d2$r2.raw.gx, d2$r2.raw.metab), c(0,0.25,0.5,0.75,1),na.rm=T)))){
    
  d2$r2.cat.gx <- "Singular"
    d2$r2.cat.metab <- "Singular"
    
    
  } else{
    
    d2$r2.cat.gx <- with(d2, 
                         cut(r2.raw.gx,
                             breaks=quantile(c(r2.raw.gx, r2.raw.metab), c(0,0.25,0.5,0.75,1),na.rm=T),
                             include.lowest = T,
                             labels=c("Low","Middle","High","Highest")))
    d2[is.na(r2.cat.gx), r2.cat.gx := "None"]
    d2$r2.cat.metab <- with(d2, 
                            cut(r2.raw.metab,
                                breaks=quantile(c(r2.raw.gx, r2.raw.metab), c(0,0.25,0.5,0.75,1),na.rm=T),
                                include.lowest = T,
                                labels=c("Low","Middle","High","Highest")))
    d2[is.na(r2.cat.metab), r2.cat.metab := "None"]
    
  }
  
  
  d2$r2.cat.gx <- as.character(d2$r2.cat.gx)
  d2$r2.cat.metab <- as.character(d2$r2.cat.metab)
  
  # Filter Data to BMI Genes ----
  
  # netzwerk based on known obesity genes
  d2f <- d2[(sig.mediation.gx == T | sig.mediation.metab == T) & gene %in% c(gene.filter), ]
  
  if(save.table==TRUE){
    
    d2f[, .(pheno, metabolite, gene, 
            gx.probe, mediation.group, 
            beta.mediation.gx, 
            se.mediation.gx = "PLACEHOLDER",
            p.mediation.gx = "PLACEHOLDER",
            r2.raw.gx,r2.raw.metab, m.path, 
            m.class1,m.class2, g.path)] %>% 
      fwrite(dpx("mediationBMIGenesTable.tsv","res/"),sep="\t", dec = ",")
    
  }
  
  return(d2f)
  
}