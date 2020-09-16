
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
#' @import ggplot2
#' @import dplyr 
#' @examples
#' MMRFGetGateway_BOresponsePlot(clinMMGateway,"Bortezomib",height=5, width=8)
#' MMRFgetGateway_BOresponsePlot(clinMMGateway,topN=40, height=15, width=15)
#' @export
#' @return table with the case count of the Best overall response to treatments




MMRFGetGateway_BOresponsePlot<- function(treat.resp,therapyname=NULL,topN=20,dpi=100,filename="BestOverall_responsePlot", height=20, width=20){
  
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
      
   
   filenm <- paste0(filename,"_", theapy.i, ".png")
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





<<<<<<< HEAD
#' @title draws plot correlating the time to Best Overall Response  (BO) leveraging the BO classification
=======
#' @title Draw plot of Time to Best Overall Response
>>>>>>> d011280e9c88104b3de929e05b2ffe323910ee19
#' @description
#' Draw plot of Time to the Best Overall Response
#' @param therapyname Therapy name
#' @param ttime cycles/days
#' @param treat.resp is a data.frame of clinical information downloaded from MMRF-Commpass Researcher Gateway 
#' and imported into environment
#' @param dpi Image dpi
#' @param filename The name of the png file
#' @param width Image width
#' @param height Image height
<<<<<<< HEAD
#' @import ggplot2
#' @import dplyr 
#' @examples
#' MMRFGetGateway_TimeBOresponsePlot(clinMMGateway,c("Bortezomib","Dexamethasone"))
#' MMRFGetGateway_TimeBOresponsePlot(clinMMGateway,c("Bortezomib","Dexamethasone"),"days")
#' MMRFGetGateway_TimeBOresponsePlot(clinMMGateway,"Bortezomib","days")
=======
#' @param dpi Image dpi
#' @import ggplot2
#' @import dplyr 
#' @examples
#' MMRFgetGateway_TimeBestOverallResponsePlot(clinMMGateway,c("Bortezomib","Dexamethasone"))
#' MMRFgetGateway_TimeBestOverallResponsePlot(clinMMGateway,c("Bortezomib","Dexamethasone"),"days")
#' MMRFgetGateway_TimeBestOverallResponsePlot(clinMMGateway,"Bortezomib","days")
>>>>>>> d011280e9c88104b3de929e05b2ffe323910ee19
#' @export
#' @return table with the case count of the Best overall response to treatments




<<<<<<< HEAD
MMRFGetGateway_TimeBOresponsePlot<- function(treat.resp,therapyname=NULL,ttime="cycles", dpi=100, filename="TimeBestOverall_responsePlot", height=10, width=10){
=======
MMRFgetGateway_TimeBestOverallResponsePlot<- function(treat.resp,therapyname=NULL,ttime="cycles", dpi=100, filename="TimeBestOverall_responsePlot", height=10, width=10){
>>>>>>> d011280e9c88104b3de929e05b2ffe323910ee19
  
  if((ttime!="cycles") & (ttime!="days")){
    
    
    stop("Please set a valid argument for time parameter: cycles/days")
  }
  
  treat.resp.sel<-subset(treat.resp, select=c("public_id","bestresp","trtname","trtshnm","ttbrespcyc", "ttbrespdy"))   
  
  
  
  
  
  
  if(!is.null(therapyname)){
    
    for (therapy.i in 1:length(therapyname)) {
      
      
      
      temp.ther<-paste0("\\",therapyname[therapy.i],"\\b", sep="")
      index.therapy<-grep(temp.ther,treat.resp$trtname,ignore.case = TRUE)
      filt<- treat.resp[index.therapy,]
      
      if(ttime=="cycles"){     
        ggplot(data = filt, aes(x=as.character(bestresp), y=ttbrespcyc)) +
          geom_boxplot(fill="steelblue") +
          labs(title=paste0(therapyname[therapy.i]," - ","Time to best overall response (cycles) by Best Overall Response"), x="Best Overall Response", y="Time to best overall response (cycles)")  + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
        
      }else {
        ggplot(data = filt, aes(x=as.character(bestresp), y=ttbrespdy)) +
          geom_boxplot(fill="steelblue") +
          labs(title=paste0(therapyname[therapy.i]," - ","Time to best overall response (days) by Best Overall Response"), x="Best Overall Response", y="Time to best overall response (cycles)")  + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
      }
      
      
      
      
      
      filenm <- paste0(filename,"_", therapy.i, ".png")
      path<-file.path(getwd())
      path<-paste0(path,"/",filenm)
      ggsave(filename =path, width = width, height = height, dpi = dpi) #save the last drawn plot
      message(paste("Plot saved in: ", file.path(getwd(),filenm)))
      
    }  
    
  } 
  else  stop("Please set a valid argument for therapyname parameter")
  
}



<<<<<<< HEAD
#' @title Draw plot of Treatment duration for the Best Overall (BO) Response filtered by Therapy classification
=======
#' @title Draw plot of Time to Best Overall Response
>>>>>>> d011280e9c88104b3de929e05b2ffe323910ee19
#' @description
#' Draw plot of the Treatment duration cycle or days
#' @param therapyname Therapy name
#' @param ttime cycles/days
#' @param line Line of therapy
#' @param treat.resp is a data.frame of clinical information downloaded from MMRF-Commpass Researcher Gateway 
#' and imported into environment
<<<<<<< HEAD
#' @param bor is the type of BO Response 
=======
#' @param bor is the type of BOR#' 
>>>>>>> d011280e9c88104b3de929e05b2ffe323910ee19
#' Example:
#' \tabular{ll}{
#'CR \tab   Complete Response \cr
#'PR \tab   Partial Response \cr
#'VGPR \tab   Very Good Partial Response \cr
#'SD \tab Stable Disease \cr
#'PD \tab  Progressive Disease \cr
#'sCR \tab   Stringent Complete Response \cr
#'}
#' @param dpi Image dpi
#' @param filename The name of the png file
#' @param width Image width
#' @param height Image height
#' @param dpi Image dpi
#' @import ggplot2
#' @import dplyr 
#' @examples
<<<<<<< HEAD
#' MMRFgetGateway_TrtBOduration(clinMMGateway,"Bortezomib",ttime="cycles",bor="PR",height=10, width=10)
#' MMRFgetGateway_TrtBOduration(clinMMGateway,"Bortezomib",ttime="days",bor="VGPR",height=10, width=10)
#' MMRFgetGateway_TrtBOduration(clinMMGateway,c("Bortezomib","Lenalidomide"),ttime="days",bor="VGPR",height=10, width=10)
=======
#' MMRFgetGateway_TrtdurationBO(clinMMGateway,"Bortezomib",ttime="cycles",bor="PR",height=10, width=10)
#' MMRFgetGateway_TrtdurationBO(clinMMGateway,"Bortezomib",ttime="days",bor="VGPR",height=10, width=10)
#' MMRFgetGateway_TrtdurationBO(clinMMGateway,c("Bortezomib","Lenalidomide"),ttime="days",bor="VGPR",height=10, width=10)
>>>>>>> d011280e9c88104b3de929e05b2ffe323910ee19
#' @export





<<<<<<< HEAD
MMRFgetGateway_TrtBOduration<- function(treat.resp,therapyname=NULL,ttime="cycles", line=1, bor="CR", dpi=100, filename="Trt_DurationPlot", height=8, width=8){
=======
MMRFgetGateway_TrtdurationBO<- function(treat.resp,therapyname=NULL,ttime="cycles", line=1, bor="CR", dpi=100, filename="Trt_DurationPlot", height=8, width=8){
>>>>>>> d011280e9c88104b3de929e05b2ffe323910ee19
  
  if(ttime!="cycles" & ttime!="days"){
    
    
    stop("Please set a valid argument for time parameter: cycles/days")
  }
  
  
  code <- c("CR","PR","VGPR","SD","PD","sCR")
  
  bestresp<-c("Complete Response","Partial Response", "Very Good Partial Response", 
              "Stable Disease", "Progressive Disease",  "Stringent Complete Response")
  
  table.bor <- data.frame(code, bestresp)
  
  resp<-table.bor[table.bor$code==bor,]
  
  if(!bor %in% table.bor$code) 
    stop("Please set a valid argument for bor parameter: CR, PR, VGPR, SD, PD, sCR")  
  
  
  
  
  
  if(!is.null(therapyname)){
    
    for (therapy.i in 1:length(therapyname)) {
      
      
      
      
      temp.ther<-paste0("\\",therapyname[therapy.i],"\\b", sep="")
      index.therapy<-grep(temp.ther,treat.resp$trtname,ignore.case = TRUE)
      filt<- treat.resp[index.therapy,]
      filt<- filt[filt$line==1 & filt$bestrespsh==bor,]
      
      
      if(ttime=="cycles"){     
        pplot<-ggplot(data = filt, aes(x=as.character(therclass), y=trtdurcyc)) +
          geom_boxplot(fill="steelblue") +
          labs(title=paste0(therapyname[therapy.i]," - "," Treatment duration (cycle) for Best Overall Response Vs Therapy classification"," (",resp$bestresp,")"), x="Therapy classification", y="Treatment duration (cycles)")  + 
          theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
        
        
        
      }else if (ttime=="days"){
        pplot<- ggplot(data = filt, aes(x=as.character(therclass), y=trtdurdy)) +
          geom_boxplot(fill="steelblue") +
          labs(title=paste0(therapyname[therapy.i]," - ","Treatment duration (days) for Best Overall Response Vs Therapy classification"," (",resp$bestresp,")"), x="Therapy classification", y="Treatment duration (days)")  + 
          theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
      }
      
      
      
      filenm <- paste0(filename,"_", therapy.i, ".png")
      path<-file.path(getwd())
      path<-paste0(path,"/",filenm)
      ggsave(filename=filenm, width = width, height = height, dpi = dpi)
      message(paste("Plot saved in: ", file.path(getwd(),filenm)))
      
      
      
      
    }  
    
  } 
  
  else  
    stop("Please set a valid argument for therapyname parameter.")
  
}













