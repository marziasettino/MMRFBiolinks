
#' @title Array Array Intensity correlation (AAIC) and correlation boxplot to define outlier
#' @description MMRFanalyzeGDC_Preprocessing perform Array Array Intensity correlation (AAIC).
#' It defines a square symmetric matrix of spearman correlation among samples.
#' According this matrix and boxplot of correlation samples by samples it is possible
#' to find samples with low correlation that can be identified as possible outliers.
#' @param object of gene expression of class RangedSummarizedExperiment from TCGAprepare
#' @param cor.cut is a threshold to filter samples according their spearman correlation in
#' samples by samples. default cor.cut is 0
#' @param filename Filename of the image file
#' @param width Image width
#' @param height Image height
#' @param datatype is a string from RangedSummarizedExperiment assay
#' @importFrom grDevices dev.list
#' @importFrom SummarizedExperiment assays
#' @importFrom TCGAbiolinks TCGAanalyze_Preprocessing
#' @export
#' @return Plot with array array intensity correlation and boxplot of correlation samples by samples

MMRFanalyzeGDC_Preprocessing <- function(object,
                                      cor.cut = 0,
                                      filename = NULL,
                                      width = 1000,
                                      height = 1000,
                                      datatype = names(assays(object))[1]){
  
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("Package \"TCGAbiolinks\" is needed. Please install it.",
         call. = FALSE)
  }
  
 TCGAanalyze_Preprocessing(object,cor.cut, filename, width,height,datatype)
  
}








#' @title survival analysis (SA) univariate with Kaplan-Meier (KM) method.
#' @description MMRFanalyzeGDC_SurvivalKM perform an univariate Kaplan-Meier (KM) survival analysis (SA).
#' It performed Kaplan-Meier survival univariate using complete follow up with all days
#' taking one gene a time from Genelist of gene symbols.
#' For each gene according its level of mean expression in cancer samples,
#' defining two thresholds for quantile
#' expression of that gene in all samples (default ThreshTop=0.67,ThreshDown=0.33) it is possible
#' to define a threshold of intensity of gene expression to divide the samples in 3 groups
#' (High, intermediate, low).
#' MMRFanalyzeGDC_SurvivalKM performs SA between High and low groups using following functions
#' from survival package
#' \enumerate{
#' \item survival::Surv
#' \item survival::survdiff
#' \item survival::survfit
#' }
#' @param clinical_patient is a data.frame using function 'clinic' with information
#' related to barcode / samples such as bcr_patient_barcode, days_to_death ,
#' days_to_last_follow_up , vital_status, etc
#' @param dataGE is a matrix of Gene expression (genes in rows, samples in cols) from MMRFGDC_prepare
#' @param Genelist is a list of gene symbols where perform survival KM.
#' @param Survresult is a parameter (default = FALSE) if is TRUE will show KM plot and results.
#' @param ThreshTop is a quantile threshold to identify samples with high expression of a gene
#' @param ThreshDown is a quantile threshold to identify samples with low expression of a gene
#' @param p.cut p.values threshold. Default: 0.05
#' @param group1 a string containing the barcode list of the samples in in control group
#' @param group2 a string containing the barcode list of the samples in in disease group
#' @importFrom survival Surv survdiff survfit
#' @importFrom TCGAbiolinks TCGAanalyze_SurvivalKM
#' @export
#' @return table with survival genes pvalues from KM.
#' @examples
#'  # Selecting only 20 genes for example
#'  dataMMcomplete <- log2(dataMM[1:20,] + 1)
#'  clinical_patient_Cancer <- MMRFqueryGDC_clinic("clinical")
#'  
#'  tabSurvKM <- MMRFanalyzeGDC_SurvivalKM(clinical_patient_Cancer,
#'                                      dataMMcomplete,
#'                                      Genelist = rownames(dataMMcomplete),
#'                                      Survresult = TRUE,
#'                                      p.cut = 0.2,
#'                                      ThreshTop = 0.67,
#'                                      ThreshDown = 0.33)
#' 
MMRFanalyzeGDC_SurvivalKM <- function(clinical_patient,
                                   dataGE,
                                   Genelist,
                                   Survresult = FALSE,
                                   ThreshTop = 0.67,
                                   ThreshDown = 0.33,
                                   p.cut = 0.05,
                                   group1,
                                   group2){
  
 
  
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("Package \"TCGAbiolinks\" is needed. Please install it.",
         call. = FALSE)
  }
  
  colnames(dataGE) <- substr(colnames(dataGE),1,9)
  tabSurvKM<-TCGAanalyze_SurvivalKM(clinical_patient,dataGE,Genelist, Survresult, ThreshTop, ThreshDown, p.cut, group1, group2)
  
  return(tabSurvKM)
  
  
}



#' @title Draw plot of the Best Overall Response to the Treatment
#' @description
#' Draw plot of the Best Overall Response to the Treatment (top 20 to make faster)
#' @param therapyname Therapy name
#' @param treat.resp is a data.frame of clinical information downloaded from MMRF-Commpass Researcher Gateway 
#' and imported into environment
#' @param topN top number of case count
#' @param dpi Image dpi
#' @param filename The name of the png file
#' @param width Image width
#' @param height Image height
#' @param dpi Image dpi
#' @import ggplot2
#' @import dplyr 
#' @examples
#' MMRFgetGateway_BestOverallResponsePlot(clinMMGateway,"Bortezomib",height=5, width=8)
#' MMRFgetGateway_BestOverallResponsePlot(clinMMGateway,topN=40, height=15, width=15)
#' @export
#' @return table with the case count of the Best overall response to treatments




MMRFgetGateway_BestOverallResponsePlot<- function(treat.resp,therapyname=NULL,topN=20,dpi=100,filename="BestOverall_responsePlot", height=20, width=20){
  
  if(!is.null(therapyname)){
    
    for (theapy.i in 1:length(therapyname)) {
     
    
    
    temp.ther<-paste0("\\",therapyname[theapy.i],"\\b", sep="")
    index.therapy<-grep(temp.ther,treat.resp$trtname,ignore.case = TRUE)
    filt<- treat.resp[index.therapy,]
    filt.group <- tally(group_by(filt,bestresp,trtshnm))
    filt.group<-filt.group[order(filt.group$n, decreasing = TRUE),]
    filt.group<-subset(filt.group[1:topN,], bestresp!="")
  
    pplot<- ggplot(filt.group, aes(bestresp,trtshnm)) +
            geom_tile(aes(fill = n), color = "steelblue") +
            scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar") +
            ylab(paste0("Treatments","-",therapyname[theapy.i])) +
            xlab("Best Overall Response") +
            theme_bw() +
            theme(text = element_text(size=11),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12),
            plot.title = element_text(size=11),
            axis.title=element_text(size=12,face="bold"),
            axis.text.x = element_text(angle = 90, hjust = 1)) +
            labs(fill = "Case Count") +
            ggtitle(paste0(therapyname," - ","Plot of Best Overall Response ", "(top case count=",topN,")"))
      
   
   filenm = paste0(filename,"_", theapy.i, ".png")
   path<-file.path(getwd())
   path<-paste0(path,"/",filenm)
   ggsave(filename =path, width = width, height = height, dpi = dpi) #save the last drawn plot
   message(paste("Plot saved in: ", file.path(getwd(),filenm)))
   
    }
   return(filt.group)
                                                                                              
  }
  
  else{
    group<- tally(group_by(treat.resp, bestresp,trtshnm))
    group <- subset(group, (bestresp!="") )
    group.count<-group[order(group$n, decreasing = TRUE),]
    group.count<-group.count[1:topN,]
    
    pplot<- ggplot(group.count, aes(bestresp,trtshnm)) +
                geom_tile(aes(fill = n), color = "steelblue") +
                scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar") +
               # scale_fill_gradient(low = "white", high = "steelblue") +
                ylab("Treatments") +
                xlab("Best Overall Response") +
                theme_bw() +
                theme(text = element_text(size=11),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 12),
                plot.title = element_text(size=11),
                axis.title=element_text(size=12,face="bold"),
                axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(fill = "Case Count") + 
                ggtitle(paste0("Plot of Best Overall Response ", "(top case count=",topN,")"))
    
    path<-file.path(getwd())
    path<-paste0(path,"/",filename,".png")
    
    ggsave(filename = path, width = width, height = height, dpi = dpi)
   
    ggsave(filename = path)
    message(paste("Plot saved in: ", file.path(getwd(),filename)))
    
    return(group.count)
    
  }
  
 
 
}



