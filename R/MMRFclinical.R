#' @title MMRFGDC_QueryClinic
#' @description
#' MMRFGDC_QueryClinic downloads and prepares all clinical information from API
#' 
#' @param type A valid type. Options "clinical", "Biospecimen"  (see list with getGDCprojects()$project_id)]
#' @param save.csv Write clinical information into a csv document
#' @export
#' @importFrom data.table rbindlist as.data.table
#' @importFrom jsonlite fromJSON
#' @importFrom TCGAbiolinks GDCquery_clinic
#' @examples
#' clin.mm<-MMRFGDC_QueryClinic(type = "clinical")
#' clin.mm<-MMRFGDC_QueryClinic(type = "Biospecimen")
#' @return A data frame with the clinical information




MMRFGDC_QueryClinic<- function(type = "clinical", save.csv = FALSE){
  
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



#' @title MMRFGDC_GetTherapyByID
#' @description
#'  Get information about Therapy filtered by sample identifier
#' @param clin.mm is a data.frame containing clinical information from GDC Data Portal 
#' (e.g.data days_to_death ,' days_to_last_follow_up , vital_status, etc)
#' @examples
#' listSamples <- c("MMRF_1951","MMRF_1474",
#'                  "MMRF_2709","MMRF_2199",
#'                  "MMRF_1216","MMRF_2119",
#'                  "MMRF_2546","MMRF_2613",
#'                  "MMRF_1647","MMRF_2170")
#' therapy.info<-MMRFGDC_GetTherapyByID(listSamples, clin.mm)
#' @export
#' @return a data.frame 




MMRFGDC_GetTherapyByID<- function(listSamples,clin.mm){ 
  
  
 
 
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







  




 


#----------------------------------------


#' @title MMRFGDC_QuerySampleTypes  
#' @description
#'   Retrieve samples identifiers from GDCquery output for filtering them by the selected type sample 
#' @param query the resuting dataframe of GDCquery
#' @param typesample a character vector indicating tissue type to query.
#' Example:
#' \tabular{ll}{
#'TBM \tab  Primary Blood Derived Cancer-Bone Marrow \cr
#'TRBM \tab Recurrent Blood Derived Cancer-Bone Marrow \cr
#'TB \tab   Primary Blood Derived Cancer-Peripheral Blood \cr
#'TRB \tab Recurrent Blood Derived Cancer - Peripheral Blood 	\cr
#'}

#' @examples
#' \dontrun{
#' query <- GDCquery(project = "MMRF-COMMPASS", 
#'                              data.category = "Transcriptome Profiling",
#'                              data.type = "Gene Expression Quantification",
#'                              experimental.strategy = "RNA-Seq",
#'                              workflow.type="HTSeq - FPKM")
#'
#'
#'tsample<-"TRBM"
#'MMRFGDC_QuerySampleTypes(query,tsample)
#'MMRFGDC_QuerySampleTypes(query,c("TB","TRBM"))
#' 
#' 
#' }
#'@return a list of samples identifiers filtered by type sample selected

MMRFGDC_QuerySampleTypes <- function(query,typesample){
  
  
  if(is.null(query) || is.null(typesample)){
    stop("Please provide arguments as explained in the integrated vignette.")
  }
  
  
  typesample.cn<- c("Primary Blood Derived Cancer - Bone Marrow",
                     "Recurrent Blood Derived Cancer - Bone Marrow", 
                      "Primary Blood Derived Cancer - Peripheral Blood",
                      "Recurrent Blood Derived Cancer - Peripheral Blood")
  
  typesample.sn <-  c("TBM","TRBM","TB","TRB")
  
  
  df.typesample <- data.frame(typesample.sn, typesample.cn)
  
  
  results<- getResults(query,cols=c("sample_type","cases"))
  
  identifier.all<-NULL
    
     for (typesample.i in 1:length(typesample)) {
       if(typesample[typesample.i] %in% df.typesample$typesample.sn ){
         temp<-df.typesample[df.typesample$typesample.sn==typesample[typesample.i],]
         filt<-results[results$sample_type==temp$typesample.cn,]
         identifier.all<-unique(substring(filt$cases,1,9))
      } else {
         return("Error message: one or more sample types do not exist")
        }
   } 
 
  return(identifier.all)
  
}


#----------------------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX------------------------
#' @title MMRFGDC_GetIdentifierByTherapy
#' @description
#' Search patient clinical information filtered by therapy name
#' @param therapyname Therapy name
#' @param clin.mm is a data.frame using function 'clinic' with information
#' related to identifier / samples such as bcr_patient_barcode, days_to_death ,
#' days_to_last_follow_up , vital_status, etc
#' @import dplyr
#' @export
#' @examples
#' clin.mm<-MMRFGDC_QueryClinic(type = "clinical")
#' identifier.dexa<-MMRFGDC_GetIdentifierByTherapy("Dexamethasone",clin.mm)
#' @return a character vector with samples identifiers




MMRFGDC_GetIdentifierByTherapy<- function(therapyname,clin.mm){  
  
  if(is.null(therapyname) || is.null(clin.mm)){
    stop("Please provide arguments as explained in the integrated vignette.")
  }
  
 
  
  treat.tab<-MMRFGetGDC_Treatments(clin.mm)
  
  for (i in 1:length(therapyname)) {
    
    df<-na.omit(treat.tab[treat.tab$treat==therapyname[i],])
    
  }
  
  return(unique(df$barcode))
}










#...................................new...............................


get_IDall<- function(query,typesample, clin.mm,therapyname){
  
  id1<-MMRFGDC_GetIdentifierByTherapy(therapyname,clin.mm)
 
  id2<-MMRFGDC_QuerySampleTypes(query,typesample)
  id<-intersect(id1,id2)
  return(id)
  
}










#' @title MMRFGDC_QuerySamples  
#' @description
#'   Retrieve samples identifiers filtered by the selected type sample (case 2) or
#'   Retrieve samples identifiers filtered by the selected therapy (case 3) or
#'   Retrieve samples identifiers from both case 1 and case 2
#' @param query the resuting dataframe of GDCquery
#' @param typesample a character vector indicating tissue type to query.
#' Example:
#' \tabular{ll}{
#'TBM \tab  Primary Blood Derived Cancer-Bone Marrow \cr
#'TRBM \tab Recurrent Blood Derived Cancer-Bone Marrow \cr
#'TB \tab   Primary Blood Derived Cancer-Peripheral Blood \cr
#'TRB \tab Recurrent Blood Derived Cancer - Peripheral Blood 	\cr
#'}
#' @export
#' @examples
#' \dontrun{
#' therapy<-"Dexamethasone" 
#' 
#' tsample<-"TRBM"
#' clin<-MMRFGDC_QueryClinic(type = "clinical")
#' 
#' query.mm <- GDCquery(project = "MMRF-COMMPASS", 
#'                              data.category = "Transcriptome Profiling",
#'                              data.type = "Gene Expression Quantification",
#'                              experimental.strategy = "RNA-Seq",
#'                              workflow.type="HTSeq - FPKM")
#'
#'
#'
#' IDs<-MMRFGDC_QuerySamples(query=query.mm,typesample=tsample, clin.mm=clin,therapyname=therapy) #case 1
#' IDs<-MMRFGDC_QuerySamples(query=query.mm,typesample=tsample) #case 2 
#' IDs<-MMRFGDC_QuerySamples(clin.mm=clin,therapyname=therapy) #case 3
#' 
#' }
#'@return a list of samples identifiers 

MMRFGDC_QuerySamples <- function(...){

  args<-list(...)
  mc <- match.call(expand.dots = TRUE )
 
  print(names(mc))


  
  
  if(length(args)== 0 || (length(args) %% 2)!=0 ){
    stop("Please provide arguments as explained in the integrated vignette.")
  }
  
  id.all<-NULL
 
  if(length(args)== 4){
    
    id.all<-get_IDall(...)
    print("Return sample ID filtered by sample type and therapy.")
    return(id.all)
    
    }
  
  
  
 
 if("query" %in% names(mc) & "typesample" %in%  names(mc) ){
   
    id.all<-MMRFGDC_QuerySampleTypes(...)
    print("Return sample ID filtered by sample type.")
    
    }
  else{
    
    if("clin.mm" %in% names(mc) & "therapyname" %in%  names(mc)){
      
      id.all<-MMRFGDC_GetIdentifierByTherapy(...)
      print("Return sample ID filtered by therapy.")
      
    } 
    
    
  }
  
  
       
 
  
  
 return(id.all)
  
}









#' @title MMRFRG_GetBorByID 
#' @description
#' Get Best Overall response filtered by sample identifier
#' @param identifier is a vector of samples identifiers
#' @param trt.resp is a data.frame of clinical information downloaded from MMRF-Commpass Researcher Gateway 
#' and imported into R environment
#' @examples
#'    list.samp<- c("MMRF_0001","MMRF_0002",
#'                  "MMRF_0003","MMRF_0004",
#'                  "MMRF_0005","MMRF_0006",
#'                  "MMRF_0007","MMRF_0008",
#'                  "MMRF_0009","MMRF_0010")
#'                  
#'    bestOverall<-MMRFRG_GetBorByID(clinMMGateway,list.samp)              
#' @return a dataframe



MMRFRG_GetBorByID<- function(treat.resp,identifier){ 
  inter<-intersect(identifier,treat.resp$public_id)  
  
  filt<-treat.resp[treat.resp$public_id %in% inter,]
  filt<-filt[,c("public_id","trtname","trtshnm","bestresp","bestrespsh")] 
  
  
  
  return(unique(filt))
}






#' @title MMRFRG_GetIDByBor 
#' @description
#' Get clinical data filtered by Best Overall Response (BOR) type 
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
#' bestOverallType<-MMRFRG_GetIDByBor(clinMMGateway,"PR")              

#' @return a dataframe


MMRFRG_GetIDByBor<- function(treat.resp, bor){ 
  
  
  code <- c("CR","PR","VGPR","SD","PD","sCR")
  
  bestresp<-c("Complete Response","Partial Response", "Very Good Partial Response", 
              "Stable Disease", "Progressive Disease",  "Stringent Complete Response")
  
  table.bor <- data.frame(code, bestresp)
  
  if(!bor %in% table.bor$code) 
    stop("Please set a valid argument for bor parameter: CR, PR, VGPR, SD, PD, sCR")
  
  
  
  
  
  filt<-treat.resp[treat.resp$bestrespsh==bor,]
  return(filt)
  
  
  
  
}











#' @title MMRFRG_GetBorInfo 
#' @description
#' Get Best Overall Response (BOR) type filtered by sample ID (Case 1) or 
#' Get sample ID filtered by BOR type (Case 2)
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
#' @param listSamples is a vector of samples identifiers
#' @export
#' @examples
#'  # case 1
#'  
#'  IDsamples<-MMRFRG_GetBorInfo(clin.rg=clinMMGateway,bor="PR") 
#'  
#'  # case 2
#'  
#'    list.samp<- c("MMRF_0001","MMRF_0002",
#'                  "MMRF_0003","MMRF_0004",
#'                  "MMRF_0005","MMRF_0006",
#'                  "MMRF_0007","MMRF_0008",
#'                  "MMRF_0009","MMRF_0010")
#'                  
#'  
#'  bestOverall<-MMRFRG_GetBorInfo(clin.rg=clinMMGateway,listSamples=list.samp)          
#' @export
#' @return a dataframe 





MMRFRG_GetBorInfo<- function(clin.rg,...){ 
  
  

args<-list(...)
mc <- match.call(expand.dots = TRUE )

print(names(mc))

res<-NULL


if(length(args)!= 1){
  stop("Please provide arguments as explained in the integrated vignette.")
}



if("listSamples" %in% names(mc)){
 
 res<-MMRFRG_GetBorByID(clin.rg,args[[1]]) 
 print("Return BOR filtered by sample ID.")
  
}
else{
  
  if("bor" %in% names(mc) ){
 
    res<- MMRFRG_GetIDByBor(clin.rg,args[[1]])
    print("Return sample ID filtered by BOR type.")
    
     } 
  
  
  } 


return(res)

} 
