
#' @title MMRFRG_SurvivalKM
#' @description Creates a survival plot from MMRF-RG patient clinical data
#' using survival library. It uses the fields D_PT_deathdy, D_PT_PRIMARYREASON and D_PT_lstalive
#' columns for groups.
#' @param patient is the data.frame of the patient clinical data downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA14_PER_PATIENT file) and imported into environment.
#' @param trt is the data.frame of the patient clinical data (i.e. treatment-response) downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP file) and imported into environment.
#' @param FilterBy Column with groups to plot. This is a mandatory field.
#' Example:
#' \tabular{ll}{
#'race \tab Race \cr
#'stage \tab ISS Stage \cr
#'treatment \tab  Treatment class \cr
#'bestresp \tab Best overall response 	\cr
#'gender \tab gender 	\cr
#'}

#' @param risk.table show or not the risk table
#' @param legend Legend title of the figure
#' @param xlim x axis limits e.g. xlim = c(0, 1000). Present narrower X axis, but not affect survival estimates.
#' @param main main title of the plot
#' @param labels labels of the plot
#' @param ylab y axis text of the plot
#' @param xlab x axis text of the plot
#' @param filename The name of the pdf file.
#' @param color Define the colors/Pallete for lines.
#' @param width Image width
#' @param height Image height
#' @param pvalue show p-value of log-rank test
#' @param conf.range  show confidence intervals for point estimates of survival curves.
#' @param dpi Figure quality
#' @import survminer
#' @import survival
#' @import gridExtra
#' @import ggplot2
#' @import stringr
#' @export
#' @return Survival plot
#' @examples
#' 
#' patient <- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                                    "MMRF_0002","MMRF_0003",
#'                                    "MMRF_0004","MMRF_0005",
#'                                    "MMRF_0006","MMRF_0007",
#'                                    "MMRF_0008","MMRF_0009"),
#'                       D_PT_PRIMARYREASON = c("Death",NA,NA,"Death","Death", 
#'                                              "Death","NA","NA","Death","Death"),
#'                       D_PT_deathdy =  c(NA,NA,2226,172,NA,NA,1328,681,786,NA),
#'                       D_PT_lstalive = c(250,356,NA,NA,1814,223,NA,NA,NA,1450),
#'                       DEMOG_GENDER = c(rep(1,5),rep(2,5)), #2 stands for female, 1 standas for male#'                       
#'                       DEMOG_ETHNICITY=c(rep("Hispanic or Latino",5),rep("Not Hispanic or Latino",5)),
#'                       D_PT_issstage_char=c(rep("Stage III",3),rep("Stage II",2),rep("Stage I",5))
#'  )
#' 
#' 
#' trt<- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                               "MMRF_0002","MMRF_0003",
#'                               "MMRF_0004","MMRF_0005",
#'                               "MMRF_0006","MMRF_0007",
#'                               "MMRF_0008","MMRF_0009"),
#'                  trtclass=c(rep("Bortezomib-based",2),rep("IMIDs-based",5),rep("combined bortezomib/IMIDs-based",3)),                                                    
#'                  bestresp=c(rep("Partial Responsed",5),rep("Very Good Partial Response",5))               
#'                                    
#'  )
#' 
#' 
#' 
#'    MMRFRG_SurvivalKM(patient,
#'                      trt,
#'                      FilterBy="treatment", 
#'                      filename=NULL,
#'                      xlim = c(100,3000),
#'                      conf.range = FALSE,
#'                      color = c("Dark2"))


MMRFRG_SurvivalKM <- function(
  patient,
  trt,
  risk.table = FALSE,
  FilterBy = NULL,
  legend = "Legend",
  labels = NULL,
  xlim = NULL,
  main = "Kaplan-Meier Survival Curve",
  ylab = "Probability of survival",
  xlab = "Time since diagnosis (days)",
  filename = "survival.pdf",
  color = NULL,
  height = 8,
  width = 12,
  dpi = 300,
  pvalue = TRUE,
  conf.range = TRUE) {
  
  
  
  if (!all(c("D_PT_PRIMARYREASON", "D_PT_lstalive","D_PT_deathdy") %in% colnames(patient)))
    stop(
      "Missing Columns D_PT_PRIMARYREASON, D_PT_lstalive and  D_PT_deathdy in survival dataframe"
    )
  
  
  
  if(is.null(trt) || (is.null(patient))){
    stop("Please provide the file of treatment-response and/or patient.")
  }
  
  
  condition <- c("race","stage","treatment","bestresp","gender")
  parameter <- c("DEMOG_ETHNICITY", "D_PT_issstage_char","trtclass","bestresp","DEMOG_GENDER")
  tab.condition <- data.frame(condition, parameter)
  
  
  names(patient)[1]<-"public_id"
  
  patient<-select(patient,public_id,D_PT_PRIMARYREASON,
                  D_PT_lstalive,D_PT_deathdy,
                  DEMOG_ETHNICITY,D_PT_issstage_char,DEMOG_GENDER
  )

  
  df.merge<-merge(x = patient, y = trt, by = "public_id", type=left)
  
  
  
  if (is.null(color)) {
    color <- rainbow(length(unique(patient[, FilterBy])))
  }
  
  group <- NULL
  
  
  
  
  if (is.null(FilterBy)) {
    stop("Please provide a value for FilterBy parameter")
  } else {
    par<-tab.condition[tab.condition$condition==FilterBy,]$parameter  #check tab.condition 
    FilterBy<-par
    
    if (length(unique(df.merge[, FilterBy])) == 1) {
      stop(
        paste0(
          "Only this group found:\n",
          unique(df.merge[, FilterBy])
        )
      )
    }
  }#end
  
  
  
  notDead <- is.na(df.merge$D_PT_deathdy)
  
  if (any(notDead == TRUE)) {
    df.merge[notDead, "D_PT_deathdy"] <- df.merge[notDead, "D_PT_lstalive"]
  }
  
  #TRUE(DEAD)/FALSE (ALIVE)
  df.merge$s <- grepl("Death", df.merge$D_PT_PRIMARYREASON, ignore.case = TRUE)
  
  
  
  df.merge$type <- as.factor(df.merge[, FilterBy])
  df.merge <-  df.merge[, c("D_PT_deathdy", "s", "type")]
  #   formula 
  f.m <-formula(survival::Surv(as.numeric(df.merge$D_PT_deathdy), event = df.merge$s) ~ df.merge$type)
  fit <- do.call(survival::survfit, list(formula = f.m, data = df.merge))
  
  label.add.n <- function(x) {
    na.idx <- is.na(df.merge[, "D_PT_deathdy"])
    negative.idx <- df.merge[, "D_PT_deathdy"] < 0
    idx <- !(na.idx | negative.idx)
    return(paste0(x, " (n = ",
                  sum(df.merge[idx, "type"] == x), ")"))
  }
  
  if (is.null(labels)) {
    d <- survminer::surv_summary(fit, data = df.merge)
    order <-
      unname(sapply(levels(d$strata), function(x)
        unlist(str_split(x, "="))[2]))
    labels <- sapply(order, label.add.n)
  }
  
  
  
  
  if (length(xlim) == 1) {
    xlim <- c(0, xlim)
  }
  
  
  suppressWarnings({
    
    
    surv <- survminer::ggsurvplot( 
      fit,
      risk.table = risk.table,
      df.merge,
      pval = pvalue,
      conf.range = conf.range,
      xlim = xlim,
      main = main,
      xlab = xlab,
      legend.title = legend,
      palette =  color,
      legend = "right",
      legend.labs = levels(FilterBy)
      
    )
  })
  
  
  
  
  if (!is.null(filename)) {
    
    filenm<-paste0(filename,".pdf")
    path<-file.path(getwd())
    path<-paste0(path,"/","ResultsPlot","/",filenm)
    
    ggsave(
      surv$plot,
      filename = path,
      device = pdf,
      width = width,
      height = height,
      units = "in"
    )
    message(paste0("File saved as: ", path))
    if (risk.table) {
      g1 <- ggplotGrob(surv$plot)
      g2 <- ggplotGrob(surv$table)
      min_ncol <- min(ncol(g2), ncol(g1))
      g <-
        gridExtra::gtable_rbind(g1[, 1:min_ncol], g2[, 1:min_ncol], size = "last")
      ggsave(
        g,
        filename = filename,
        width = width,
        height = height,
        dpi = dpi
      )
    }
  } else {
    return(surv)
  }
}






#--------------------------------------------------------------------------------

#' @title survival analysis (SA) univariate with Kaplan-Meier (KM) method.
#' @description MMRFGDC_SurvivalKM performs an univariate Kaplan-Meier (KM) survival analysis (SA)
#' between High and Low groups. 
#' @param clinical is the clinical data.frame obtained from MMRFGDC_QueryClinic
#' @param MatrixGE is a matrix of Gene expression (genes in rows, samples in cols) from MMRFGDC_prepare
#' @param ListGenes is a list of gene symbols.
#' @param Results is a parameter (default = FALSE) if is TRUE the KM plot is shown.
#' @param High is a quantile threshold to identify samples with high expression of a gene
#' @param Low is a quantile threshold to identify samples with low expression of a gene
#' @param p.cut p.values threshold. Default: 0.05
#' @param group1 a string containing the barcode list of the samples in in control group
#' @param group2 a string containing the barcode list of the samples in in disease group
#' @export
#' @return table with survival genes pvalues from KM.
#' @examples

#'
#' clin.mm <- MMRFGDC_QueryClinic(type = "clinical")
#'  
#' query.exp.count <- GDCquery(project = "MMRF-COMMPASS", 
#'                             data.category = "Transcriptome Profiling",
#'                             data.type = "Gene Expression Quantification",
#'                             experimental.strategy = "RNA-Seq",
#'                             workflow.type="HTSeq - Counts")
#'                             
#'   GDCdownload(query.exp.count)
#'  
#'  MMRFdata.prep.sub <- MMRFGDC_prepare(query.exp.count,
#'                                       save = TRUE ,
#'                                       save.filename = "data/RNASeqSE.rda" ,
#'                                       directory = "GDCdata",
#'                                       summarizedExperiment = TRUE)
#'  
#' ListBarcode1<-MMRFGDC_QuerySamples(query=query.exp.count,typesample="TBM") 
#' ListBarcode2<-MMRFGDC_QuerySamples(query=query.exp.count,typesample="TRBM")
#' 
#' 
#' 
#' MMRFdata.prep.sub.df<-assay(MMRFdata.prep.sub)
#' temp<-MMRFdata.prep.sub.df[, unique(colnames(MMRFdata.prep.sub.df))]
#' 
#'  Selecting only 20 genes for example
#'  dataMM.log <- log2(temp[1:20,] + 1)
#' 
#' SurvKM <- MMRFGDC_SurvivalKM(clin.mm,
#'                              dataMM.log,
#'                              ListGenes = rownames(dataMM.log),
#'                              Results = FALSE,
#'                              High = 0.67,
#'                              Low = 0.33,
#'                              p.cut = 0.05,
#'                              ListBarcode1,
#'                              ListBarcode2)
#'  
  
MMRFGDC_SurvivalKM <- function(
  clinical,
  MatrixGE,
  ListGenes,
  Results = FALSE,
  High = 0.67,
  Low = 0.33,
  p.cut = 0.05,
  group1,
  group2
) {
  
  # check_package("survival")
  # Check which genes we really have in the matrix
  ListGenes <- intersect(rownames(MatrixGE), ListGenes)
  
  # Split gene expression matrix btw the groups
  dataGroup1 <- MatrixGE[ListGenes, group1, drop = FALSE]
  dataGroup2 <- MatrixGE[ListGenes, group2, drop = FALSE]
  
  colnames(dataGroup1)  <- substr(colnames(dataGroup1), 1, 12)
  
  cfu <-  clinical[clinical[, "bcr_patient_barcode"] %in% substr(colnames(dataGroup1), 1, 12), ]
  if ("days_to_last_followup" %in% colnames(cfu))
    colnames(cfu)[grep("days_to_last_followup", colnames(cfu))] <-  "days_to_last_follow_up"
  cfu <-as.data.frame(subset(
    cfu,
    select = c(
      "bcr_patient_barcode",
      "days_to_death",
      "days_to_last_follow_up",
      "vital_status"
    )
  ))
  
  # Set alive death to inf
  if (length(grep("alive", cfu$vital_status, ignore.case = TRUE)) > 0)
    cfu[grep("alive", cfu$vital_status, ignore.case = TRUE), "days_to_death"] <-
    "-Inf"
  
  # Set dead follow up to inf
  if (length(grep("dead", cfu$vital_status, ignore.case = TRUE)) > 0)
    cfu[grep("dead", cfu$vital_status, ignore.case = TRUE), "days_to_last_follow_up"] <-
    "-Inf"
  
  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]
  
  followUpLevel <- FALSE
  
  #FC_FDR_table_mRNA
  tabSurv_Matrix <-
    matrix(0, nrow(as.matrix(rownames(dataGroup2))), 8)
  colnames(tabSurv_Matrix) <- c(
    "mRNA",
    "pvalue",
    "Cancer Deaths",
    "Cancer Deaths with High",
    "Cancer Deaths with Low",
    "Mean Tumor High",
    "Mean Tumor Low",
    "Mean Normal"
  )
  
  tabSurv_Matrix <- as.data.frame(tabSurv_Matrix)
  
  cfu$days_to_death <- as.numeric(as.character(cfu$days_to_death))
  cfu$days_to_last_follow_up <-  as.numeric(as.character(cfu$days_to_last_follow_up))
  rownames(cfu) <- cfu[, "bcr_patient_barcode"] #mod1
  
  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]
  
  cfu_complete <- cfu
  ngenes <- nrow(as.matrix(rownames(dataGroup2)))
  
  # Evaluate each gene
  for (i in 1:nrow(as.matrix(rownames(dataGroup2))))  {
    cat(paste0((ngenes - i), "."))
    mRNAselected <- as.matrix(rownames(dataGroup2))[i]
    mRNAselected_values <-  dataGroup1[rownames(dataGroup1) == mRNAselected, ]
    mRNAselected_values_normal <- dataGroup2[rownames(dataGroup2) == mRNAselected, ]
    if (all(mRNAselected_values == 0))
      next # All genes are 0
    tabSurv_Matrix[i, "mRNA"] <- mRNAselected
    
    
    # Get Thresh values for cancer expression
    mRNAselected_values_ordered <-  sort(mRNAselected_values, decreasing = TRUE)
    mRNAselected_values_ordered_High <-  as.numeric(quantile(as.numeric(mRNAselected_values_ordered), High)[1])
    mRNAselected_values_ordered_Low <- as.numeric(quantile(as.numeric(mRNAselected_values_ordered), Low)[1])
    
    mRNAselected_values_newvector <- mRNAselected_values
    
    
    if (!is.na(mRNAselected_values_ordered_High)) {
      # How many samples do we have
      numberOfSamples <- length(mRNAselected_values_ordered)
      
      
      skip_to_next <- FALSE 
      
      tryCatch({   
        
        # High group (above High)
        lastelementHigh <-
          max(which(
            mRNAselected_values_ordered > mRNAselected_values_ordered_High
          ))
        
      }, error = function(e) { skip_to_next <- TRUE} )     
      
      
      if(skip_to_next) { next } 
      
      
      
      
      # Low group (below Low)
      firstelementLow <-
        min(
          which(
            mRNAselected_values_ordered <= mRNAselected_values_ordered_Low
          )
        )
      
      skip_to_next <- FALSE
      tryCatch(samples_High_mRNA_selected <- names(mRNAselected_values_ordered[1:lastelementHigh]), 
               error = function(e) { skip_to_next <- TRUE})
      
      if(skip_to_next) { next } 
      
      
      
      samples_Low_mRNA_selected <- names(mRNAselected_values_ordered[firstelementLow:numberOfSamples])
      
      # Which samples are in the intermediate group (above ThreshLow and below High)
      samples_UNCHANGED_mRNA_selected <-
        names(mRNAselected_values_newvector[which((mRNAselected_values_newvector) > mRNAselected_values_ordered_Low & mRNAselected_values_newvector < mRNAselected_values_ordered_High
        )])
      
      cfu_onlyHigh <-cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_High_mRNA_selected, ]
      cfu_onlyLow <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_Low_mRNA_selected, ]
      cfu_onlyUNCHANGED <-cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_UNCHANGED_mRNA_selected, ]
      
      cfu_ordered <- NULL
      cfu_ordered <- rbind(cfu_onlyHigh, cfu_onlyLow)
      cfu <- cfu_ordered
      
      ttime <- as.numeric(cfu[, "days_to_death"])
      
      sum(status <- ttime > 0) # morti
      deads_complete <- sum(status <- ttime > 0)
      
      ttime_only_High <- cfu_onlyHigh[, "days_to_death"]
      deads_High <- sum(ttime_only_High > 0)
      
      if (dim(cfu_onlyLow)[1] >= 1) {
        ttime_only_Low <- cfu_onlyLow[, "days_to_death"]
        deads_Low <- sum(ttime_only_Low > 0)
      } else {
        deads_Low <- 0
      }
      
      tabSurv_Matrix[i, "Cancer Deaths"] <- deads_complete
      tabSurv_Matrix[i, "Cancer Deaths with High"] <- deads_High
      tabSurv_Matrix[i, "Cancer Deaths with Low"] <-  deads_Low
      tabSurv_Matrix[i, "Mean Normal"] <-  mean(as.numeric(mRNAselected_values_normal))
      dataGroup1_onlyHigh_sample <-dataGroup1[, samples_High_mRNA_selected, drop = FALSE]
      dataGroup1_onlyHigh_sample_mRNASelected <- dataGroup1_onlyHigh_sample[rownames(dataGroup1_onlyHigh_sample) == mRNAselected, ]
      dataGroup1_onlyLow_sample <- dataGroup1[, samples_Low_mRNA_selected, drop = FALSE]
      dataGroup1_onlyLow_sample_mRNASelected <-dataGroup1_onlyLow_sample[rownames(dataGroup1_onlyLow_sample) == mRNAselected, ]
      tabSurv_Matrix[i, "Mean Tumor High"] <-  mean(as.numeric(dataGroup1_onlyHigh_sample_mRNASelected))
      tabSurv_Matrix[i, "Mean Tumor Low"] <- mean(as.numeric(dataGroup1_onlyLow_sample_mRNASelected))
      
      ttime[!status] <- as.numeric(cfu[!status, "days_to_last_follow_up"])
      ttime[which(ttime == -Inf)] <- 0
      
      ttime <- survival::Surv(ttime, status)
      rownames(ttime) <- rownames(cfu)
      legendHigh <- paste(mRNAselected, "High")
      legendLow  <- paste(mRNAselected, "Low")
      
      tabSurv_pvalue <- tryCatch({
        tabSurv <-
          survival::survdiff(ttime  ~ c(rep(
            "High", nrow(cfu_onlyHigh)
          ), rep(
            "Low", nrow(cfu_onlyLow)
          )))
        tabSurv_chis <- unlist(tabSurv)$chisq
        tabSurv_pvalue <- as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
      }, error = function(e) {
        return(Inf)
      })
      tabSurv_Matrix[i, "pvalue"] <- tabSurv_pvalue
      
      if (Results == TRUE) {
        titlePlot <-
          paste("Kaplan-Meier Survival analysis, pvalue=",
                tabSurv_pvalue)
        plot(
          survival::survfit(ttime ~ c(
            rep("low", nrow(cfu_onlyHigh)), rep("high", nrow(cfu_onlyLow))
          )),
          col = c("green", "red"),
          main = titlePlot,
          xlab = "Days",
          ylab = "Survival"
        )
        legend(
          100,
          1,
          legend = c(legendLow, legendHigh),
          col = c("green", "red"),
          text.col = c("green", "red"),
          pch = 15
        )
        print(tabSurv)
      }
    } #end if
    
  } #end for
  
  tabSurv_Matrix[tabSurv_Matrix == "-Inf"] <- 0
  
  tabSurvKM <- tabSurv_Matrix
  
  # Filtering by selected pvalue < 0.01
  tabSurvKM <- tabSurvKM[tabSurvKM$mRNA != 0, ]
  tabSurvKM <- tabSurvKM[tabSurvKM$pvalue < p.cut, ]
  tabSurvKM <- tabSurvKM[!duplicated(tabSurvKM$mRNA), ]
  rownames(tabSurvKM) <- tabSurvKM$mRNA
  tabSurvKM <- tabSurvKM[, -1]
  tabSurvKM <-
    tabSurvKM[order(tabSurvKM$pvalue, decreasing = FALSE), ]
  
  colnames(tabSurvKM) <-  gsub("Cancer", "Group2", colnames(tabSurvKM))
  colnames(tabSurvKM) <- gsub("Tumor", "Group2", colnames(tabSurvKM))
  colnames(tabSurvKM) <-  gsub("Normal", "Group1", colnames(tabSurvKM))
  
  
  return(tabSurvKM)
}










