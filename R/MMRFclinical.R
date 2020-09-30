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
#' clin.mm<-MMRFqueryGDC_clinic(type = "Biospecimen")
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
#' get treatment list from clinical data
#' @param clin.mm is a data.frame containing clinical information from GDC Data Portal 
#' (e.g.data days_to_death ,' days_to_last_follow_up , vital_status, etc)
#' @examples
#' treat.list<-MMRFGetGDC_Treatments(clin.mm)
#' @export
#' @return a data.frame 




MMRFGetGDC_Treatments<- function(clin.mm){ 
  
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



#' @title Get information about Therapy filtered by sample identifier
#' @description
#' get Therapy list from clinical data
#' @param clin.mm is a data.frame containing clinical information from GDC Data Portal 
#' (e.g.data days_to_death ,' days_to_last_follow_up , vital_status, etc)
#' @examples
#' listSamples <- c("MMRF_1951","MMRF_1474",
#'                  "MMRF_2709","MMRF_2199",
#'                  "MMRF_1216","MMRF_2119",
#'                  "MMRF_2546","MMRF_2613",
#'                  "MMRF_1647","MMRF_2170")
#' therapy.info<-MMRFGetGDC_TherapybyIdentifier(listSamples, clin.mm)
#' @export
#' @return a data.frame 




MMRFGetGDC_TherapyByIdentifier<- function(listSamples,clin.mm){ 
  
  
 
 
  clin.mm<-clin.mm[,c("bcr_patient_barcode","treatments")]
  
  
  clin.mm.filt<-NULL
  clin.mm.trt<-NULL
  
  
  for(i in 1:length(listSamples)){
    clin.mm.filt.aux<-clin.mm[clin.mm$bcr_patient_barcode==listSamples[i],]
     
     clin.mm.filt<-rbind(clin.mm.filt, clin.mm.filt.aux)
  } 
  
  
  for(i in 1:length(clin.mm.filt)){
    
   # clin.mm.bar<-clin.mm.aux$bcr_patient_barcode
    clin.mm.trt.aux<-clin.mm.filt$treatments[[i]]
   
    
    clin.mm.trt<-rbind(clin.mm.trt,clin.mm.trt.aux)
    
  }  
  
  clin.mm.trt<-clin.mm.trt[,c("submitter_id","therapeutic_agents","regimen_or_line_of_therapy","days_to_treatment_start","days_to_treatment_end")]
  
  
  return(clin.mm.trt)
}





#' @title Search patient clinical information filtered by therapy name
#' @description
#' Search patient clinical information filtered by therapy name
#' @param therapyname Therapy name
#' @param clin.mm is a data.frame using function 'clinic' with information
#' related to identifier / samples such as bcr_patient_barcode, days_to_death ,
#' days_to_last_follow_up , vital_status, etc
#' @examples
#' identifier.dexa<-MMRFGetGDC_IdentifierByTherapy("Dexamethasone",clin.mm)
#' @export
#' @return a character vector with samples identifiers




MMRFGetGDC_IdentifierByTherapy<- function(therapyname,clin.mm){  
  
  treat.tab<-MMRFGetGDC_Treatments(clin.mm)
  identifier<-NULL
  
  for (i in 1:length(therapyname)) {
   
    df<-filter(treat.tab,treat==therapyname[i])
    identifier<-union(identifier,df$barcode)
 }
  
  return(unique(identifier))
}




  




  
#' @title Get Best Overall response filtered by sample identifier
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
#' bestOverall<-MMRFGetGateway_BOresponse(listSamples, clinMMGateway)              
#' @export
#' @return a dataframe



MMRFGetGateway_BOresponse<- function(identifier,treat.resp){ 
  inter<-intersect(identifier,treat.resp$public_id)  
  
  filt<-treat.resp[treat.resp$public_id %in% inter,]
  filt<-filt[,c("public_id","trtname","trtshnm","bestresp","bestrespsh")] 
  
  
  
  return(filt)
}






#' @title filter clinical data by Best Overall Response (BOR) type 
#' @description
#' filter trt.resp by Best Overall Response (BOR) type 
#' @param bor is the type of Best Overall Response (BOR)
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
#' bestOverallType<-MMRFGetGateway_BOresponseType(clinMMGateway,"PR")              
#' @export
#' @return a dataframe





MMRFGetGateway_BOresponseType<- function(treat.resp, bor){ 
  
  
  code <- c("CR","PR","VGPR","SD","PD","sCR")
  
  bestresp<-c("Complete Response","Partial Response", "Very Good Partial Response", 
              "Stable Disease", "Progressive Disease",  "Stringent Complete Response")
  
  table.bor <- data.frame(code, bestresp)
  
 if(!bor %in% table.bor$code) 
   stop("Please set a valid argument for bor parameter: CR, PR, VGPR, SD, PD, sCR")
  
  
  
 
  
 filt<-treat.resp[treat.resp$bestrespsh==bor,]
 return(filt)
 
 
  
 
}


