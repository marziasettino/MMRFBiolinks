#variant.ann<-MMRF_CoMMpass_IA14a_All_Canonical_Variants
#patient<-MMRF_CoMMpass_IA14_PER_PATIENT
#trt<-MMRF_CoMMpass_IA14_STAND_ALONE_TRTRESP
#topN<-20
#filenm<-"VariantCountPlot"

#width <-10
#height <- 10
#topN<-50
#filenm<-NULL



#-------------
# variant <- c("rs377332977", "rs372745823","rs186634824","rs371179860")




#-----------------------------------------------------
#' @title draws plot of annoteted variants by Best Overall Response and Therapy class
#' @description
#' Draw heatmap of annotated count of variants (i.e. count of dbSNP variants for each sample)
#' @param topN is the top number of variant count
#' @param variant.ann is the data.frame of annotated variants downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA14a_All_Canonical_Variants file) and imported into environment
#' @param filenm is the name of the png file. If filenm is Null, the plot is draw but it is not saved.
#' @param width Image width
#' @param height Image height
#' @import ggplot2
#' @import dplyr 
#' @examples
#' variant.ann<- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                                      "MMRF_0002","MMRF_0003",
#'                                      "MMRF_0004","MMRF_0005",
#'                                      "MMRF_0006","MMRF_0007",
#'                                      "MMRF_0008","MMRF_0009"),
#'                  ID=c(rep("rs755588843",2),rep("rs569344016",5),rep("rs2066497",2),rep(".",1)),                                                    
#'                  mutType=c(rep("intragenic_variant",3),
#'                            rep("missense_variant",2),
#'                            rep("intron_variant",1),
#'                            rep("5_prime_UTR_variant",4))             
#'                                  
#'  )
#' trt<- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                               "MMRF_0002","MMRF_0003",
#'                               "MMRF_0004","MMRF_0005",
#'                               "MMRF_0006","MMRF_0007",
#'                               "MMRF_0008","MMRF_0009"),
#'                  trtclass=c(rep("Bortezomib-based",2),rep("IMIDs-based",5),rep("combined bortezomib/IMIDs-based",3)),                                                    
#'                  bestresp=c(rep("Partial Responsed",2),rep("Very Good Partial Response",2),
#'                             rep("Progressive Disease",3),rep("Stable Disease",2),rep("Stringent Complete Response",1))               
#'                                    
#'  )
#' 
#' 
#' 
#' 
#' 
#' summary.var<-MMRF_RG_VariantCountPlot(variant.ann,trt,topN=50,filenm=NULL)
#' @export
#' @return plot with the top count of the dbSNP variant categorized by Best Overall Response and Therapy class




MMRF_RG_VariantCountPlot<- function(variant.ann, trt, topN=20,filenm="VariantCountPlot", height=10, width=10){
  
  if(is.null(variant.ann) || is.null(trt)){
    stop("Please provide the file of the annotated variants and treatment-response files.")
  }
  
 
  names(variant.ann)[1]<-"public_id"
  variant.ann$public_id<-substr(variant.ann$public_id,1,9)
  
 if("ANN....EFFECT" %in% names(variant.ann)){
   names(variant.ann)[50]<-"mutType"
 }
  
  
  variant.ann<-select(variant.ann,public_id,ID,mutType)
  variant.ann<-subset(variant.ann, variant.ann$ID != "." & !is.nan(variant.ann$ID))
  variant.ann<-unique(variant.ann)
  
  df.merge<-merge(x = variant.ann, y = trt, by = "public_id", type=left)
 
  
  variant.summary<-df.merge %>% group_by(ID,trtclass,bestresp) %>% summarize(n())
  names(variant.summary)[4]<-"count"
 
 
  
  variant.summary<-variant.summary[order(variant.summary$count, decreasing = TRUE),]
 
  
  plt<-ggplot(head(variant.summary,topN), aes(x = count, y = ID,shape=trtclass)) + 
    geom_point() + facet_grid(~bestresp) +
    theme(text = element_text(size=9),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.title = element_text(size=11),
          axis.title=element_text(size=12,face="bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle(paste0("Variants count by Best Overall Response and  \n Therapy class ", "(top Variant count=",topN,")"))+
    xlab(label = "N variant")+
    ylab("dbSNP ID")
  
  
  #filenm<-filenm
  
                 
  if (!is.null(filenm)) {
    
    filenm<-paste0(filenm,".pdf")
    path<-file.path(getwd())
    path<-paste0(path,"/","ResultsPlot","/",filenm)
    
    ggsave(path, device = pdf, width = width, height = height, units = "in")
    message(paste0("File saved as: ", path))
   
  } else {
    print(plt)
  }
 
return(variant.summary)
}










#' @title Filter patient information  by dbSNP variant 
#' @description
#' Filter patient information  by dbSNP variant 
#' @param variant is the vector of dbSNP ID
#' @param variant.ann is the data.frame of annotated variants downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA14a_All_Canonical_Variants file) and imported into environment
#' @param patient is the data.frame of the patient clinical data downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA14_PER_PATIENT file) and imported into environment.
#' @import dplyr 
#' @examples
#'  variant <- c("rs755588843", "rs569344016","rs2066497")
#'  patient.var<-MMRF_RG_GetIDSamplebyVariant(variant.ann,patient,variant)
#' @export
#' @return dataframe of patient information filtered by dbSNP variant 

MMRF_RG_GetIDSamplebyVariant<- function(variant.ann, patient, variant){
  
  if(is.null(variant) || is.null(variant.ann) || is.null(patient)){
    stop("Please provide the patient  or variant file.")
  }else {
      if(is.null(variant)){
        stop("Please provide a valid list of dbSNP ID.")
    }
    
}

  id_samples<-NULL
  
  names(variant.ann)[1]<-"public_id"
  variant.ann$public_id<-substr(variant.ann$public_id,1,9)
  
  for (rs.i in 1:length(variant)) {
    
    var<-variant[rs.i]
    id<-variant.ann[variant.ann$ID==var,]$public_id
    id_samples<-union(id_samples,id)
    
  }    
  
  
  names(patient)[1]<-"public_id"
  
  df<-variant.ann[variant.ann$public_id %in% id_samples,]
  
  
  df.merge<-merge(x = patient, y = df, by = "public_id", type=left)
  
  
  return(df.merge)
}























