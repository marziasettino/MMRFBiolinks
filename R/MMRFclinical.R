#' @title Get GDC clinical data
#' @description
#' MMRFqueryGDC_clinic will download and prepare all clinical information from the API
#' as the one with using the button from each project
#' @param type A valid type. Options "clinical", "Biospecimen"  (see list with getGDCprojects()$project_id)]
#' @param save.csv Write clinical information into a csv document
#' @export
#' @importFrom data.table rbindlist as.data.table
#' @importFrom jsonlite fromJSON
#' @importFrom TCGAbiolinks GDCquery_clinic
#' @examples
#' clin.mm<-MMRFqueryGDC_clinic(type = "clinical")
#' clin.mm<-MMRFqueryGDC_clinic(type = "biospecimen")
#' @return A data frame with the clinical information




MMRFqueryGDC_clinic<- function(type = "clinical", save.csv = FALSE){
  
  clin<-GDCquery_clinic(project="MMRF-COMMPASS", type, save.csv)
  
  #names(clin)[names(clin) == 'submitter_id'] <- 'bcr_patient_barcode'
  
  #if clinical
  
  if (type!="clinical" & type!="Biospecimen" ){
    return("Error message: type does not exist")
  }
  
  if (type=="clinical"){
    names(clin)[1] <- "bcr_patient_barcode"
    
  }
  else if  (type=="Biospecimen"){
    bcr_patient_barcode.sub<-clin[colnames(clin)=='submitter_id']
    colnames(bcr_patient_barcode.sub)[1]<-"bcr_patient_barcode"
    bcr_patient_barcode.sub$bcr_patient_barcode<-substr(bcr_patient_barcode.sub$bcr_patient_barcode,1,9)
    clin<-cbind(bcr_patient_barcode.sub,clin)
    
  }
  
  
  return(clin)
  
}






#' @title Get list of treatments 
#' @description
#' get patient clinical information filtered by therapy name
#' @param clin.mm is a data.frame using function 'clinic' with information
#' related to barcode / samples such as bcr_patient_barcode, days_to_death ,
#' days_to_last_follow_up , vital_status, etc
#' @export
#' @return a data.frame 




MMRFgetGDC_Treatments<- function(clin.mm){ 
  
  df<-NULL
  
  treat.list<-clin.mm$treatments
  for(i in 1:length(treat.list)){
    treat.aux<-as.data.frame(treat.list[[i]])
    barcode<-substring(treat.aux$submitter_id,1,9)
    treat<-treat.aux$therapeutic_agents
    line<-treat.aux$regimen_or_line_of_therapy
    df<-rbind(df, data.frame(barcode,treat,line))
    
  }  
  
  
  return(unique(df))
}





#' @title Search patient clinical information filtered by therapy name
#' @description
#' Search patient clinical information filtered by therapy name
#' @param therapyname Therapy name
#' @param clin.mm is a data.frame using function 'clinic' with information
#' related to barcode / samples such as bcr_patient_barcode, days_to_death ,
#' days_to_last_follow_up , vital_status, etc
#' @examples
#' bar.dexa<-MMRFgetGDC_BarcodeTherapy("Dexamethasone",clin.mm)
#' @export
#' @return a character vector with samples identifiers




MMRFgetGDC_BarcodeTherapy<- function(therapyname,clin.mm){  
  
  treat.tab<-MMRFgetGDC_Treatments(clin.mm)
  barcode<-NULL
  
  for (i in 1:length(therapyname)) {
    df<-filter(treat.tab,treat==therapyname[i])
    barcode<-union(barcode,df$barcode)
  }
  
  return(unique(barcode))
}




  




  
#' @title Get Best Overall response 
#' @description
#' filter trt.resp by samples identifiers
#' @param identifier is a vector of samples identifiers
#' @param trt.resp is a data.frame of clinical information downloaded from MMRF-Commpass Researcher Gateway 
#' and imported into R environment
#' @examples
#' listSamples <- c("MMRF_001","MMRF_002",
#'                  "MMRF_003","MMRF_003",
#'                  "MMRF_004","MMRF_005",
#'                  "MMRF_006","MMRF_007",
#'                  "MMRF_008","MMRF_009")
#'                  
#' bestOveall<-MMRFgetGateway_BestOverallResponse(listSamples, clinMMGateway)              
#' @export
#' @return a dataframe



MMRFgetGateway_BestOverallResponse<- function(identifier,treat.resp){ 
  inter<-intersect(identifier,treat.resp$public_id)  
  
  filt<-treat.resp[treat.resp$public_id %in% inter,]
  filt<-filt[,c("public_id","trtname","trtshnm","bestresp","bestrespsh")] 
  
  
  
  return(filt)
}






#' @title Get Best Overall response type 
#' @description
#' filter trt.resp by Best Overall Response (BOR) type 
#' @param bor is the type of BOR#' 
#' Example:
#' \tabular{ll}{
#'CR \tab   Complete Response \cr
#'PR \tab   Partial Response \cr
#'VGPR \tab   Very Good Partial Response \cr
#'SD \tab Stable Disease \cr
#'PD \tab  Progressive Disease \cr
#'sCR \tab   Stringent Complete Response \cr
#'}
#' @param trt.resp is a data.frame of clinical information downloaded from MMRF-Commpass Researcher Gateway 
#' and imported into R environment
#' @examples
#' bestOveallType<-MMRFgetGateway_BestOverallResponseType(clinMMGateway,"PR" )              
#' @export
#' @return a dataframe





MMRFgetGateway_BestOverallResponseType<- function(treat.resp, bor){ 
  
  
  code <- c("CR","PR","VGPR","SD","PD","sCR")
  
  bestresp<-c("Complete Response","Partial Response", "Very Good Partial Response", 
              "Stable Disease", "Progressive Disease",  "Stringent Complete Response")
  
  table.bor <- data.frame(code, bestresp)
  
 if(!bor %in% table.bor$code) 
   stop("Please set a valid argument for bor parameter: CR, PR, VGPR, SD, PD, sCR")
  
  
  
 
  
 filt<-treat.resp[treat.resp$bestrespsh==bor,]
 return(filt)
 
 
  
 
}








