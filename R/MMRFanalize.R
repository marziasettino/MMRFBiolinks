
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
#' MMRFGetGateway_BOresponsePlot(clinMMGateway,topN=40, height=13, width=13)
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
            axis.text.x = element_text(angle = 45, hjust = 1)) +
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
                theme(text = element_text(size=10),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 12),
                plot.title = element_text(size=11),
                axis.title=element_text(size=12,face="bold"),
                axis.text.x = element_text(angle = 45, hjust = 1)) +
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





#' @title draws plot correlating the time for the Best Overall Response  (BO) to BO classification
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
#' @import ggplot2
#' @import dplyr 
#' @examples
#' MMRFGetGateway_TimeBOresponsePlot(clinMMGateway,c("Bortezomib","Dexamethasone"))
#' MMRFGetGateway_TimeBOresponsePlot(clinMMGateway,c("Bortezomib","Dexamethasone"),"days")
#' MMRFGetGateway_TimeBOresponsePlot(clinMMGateway,"Bortezomib","days")
#' @export
#' @return table with the case count of the Best overall response to treatments




MMRFGetGateway_TimeBOresponsePlot<- function(treat.resp,therapyname=NULL,ttime="cycles", dpi=100, filename="TimeBestOverall_responsePlot", height=10, width=10){
  
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
          labs(title=paste0(therapyname[therapy.i]," - ","Time to best overall response by Best Overall Response"), x="Best Overall Response", y="Time to best overall response (cycles)")  + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
        
      }else {
        ggplot(data = filt, aes(x=as.character(bestresp), y=ttbrespdy)) +
          geom_boxplot(fill="steelblue") +
          labs(title=paste0(therapyname[therapy.i]," - ","Time to best overall response by Best Overall Response"), x="Best Overall Response", y="Time to best overall response (days)")  + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
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



#' @title Draw plot of Treatment duration filtered by Therapy classification
#' @description
#' Draw plot of the Treatment duration cycle or days
#' @param therapyname Therapy name
#' @param ttime cycles/days
#' @param line Line of therapy
#' @param treat.resp is a data.frame of clinical information downloaded from MMRF-Commpass Researcher Gateway 
#' and imported into environment
#' @param bor is the type of BO Response 
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
#' MMRFGetGateway_TrtBOdurationPlot(clinMMGateway,"Bortezomib",ttime="cycles",bor="PR",height=10, width=10)
#' MMRFGetGateway_TrtBOdurationPlot(clinMMGateway,"Bortezomib",ttime="days",bor="VGPR",height=10, width=10)
#' MMRFGetGateway_TrtBOdurationPlot(clinMMGateway,c("Bortezomib","Lenalidomide"),ttime="days",bor="VGPR",height=10, width=10)
#' @export





MMRFGetGateway_TrtBOdurationPlot<- function(treat.resp,therapyname=NULL,ttime="cycles", line=1, bor="CR", dpi=100, filename="Trt_DurationPlot", height=8, width=8){
  
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













