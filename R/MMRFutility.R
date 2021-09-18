
#' @title MMRFGDC_ProjectSummary 
#' @description
#'  Get information about MMRF-COMMPASS Project form GDC Data Portal 
#' @import jsonlite
#' @export 
#' @examples
#' \dontrun{
#' summary<-MMRFGDC_ProjectSummary()
#' summary$data_categories
#' summary$experimental_strategies
#' 
#' 
#' }
#' @return a list

MMRFGDC_ProjectSummary  <- function(){
  project<-"MMRF-COMMPASS"
  baseURL <- "https://api.gdc.cancer.gov/projects/"
  url <- paste0(baseURL, project,"?expand=summary,summary.experimental_strategies,summary.data_categories&pretty=true")
  return(fromJSON(url,simplifyDataFrame = TRUE)$data$summary)
}


#' @title MMRFGDC_QuerySummary
#' @description
#'  get summarized information about data (open access) obtained from GDCquery function 
#' @param query A query form GDCquery function
#' @import dplyr
#' @export 
#' @examples
#' \dontrun{
#' query.mm.fpkm <- GDCquery(project = "MMRF-COMMPASS",
#'                           data.category = "Transcriptome Profiling",
#'                           data.type = "Gene Expression Quantification",
#'                           workflow.type="HTSeq - FPKM")
#' summary<-MMRFGDC_QuerySummary(query.mm.fpkm)
#' 
#' 
#' }
#' @return a data.frame


MMRFGDC_QuerySummary <- function(query){
  
 res.query<-getResults(query)
  
  if(query$data.category=="Transcriptome Profiling"){
    
    query.group<-group_by(res.query,data_type,experimental_strategy,sample_type,analysis_workflow_type)
    summary <- summarize(query.group, n_cases = n())
    
    } else if(query$data.category=="Simple Nucleotide Variation") {
    
    query.group<-group_by(res.query,experimental_strategy,analysis_workflow_type)
    summary <- summarize(query.group, n_cases = n())
    
  } else if(query$data.category=="Sequencing Reads") {
   
   query.group<-group_by(res.query,data_type,experimental_strategy,sample_type,analysis_workflow_type)
   summary <- summarize(query.group, n_cases = n())
   
 }
 
 else {
   stop("Query format in not allowed")
 }
 
 
  return(summary)
  
}



#' @title Convert from ensembl.gene to gene.symbol annotation
#' @description
#'  get information about query obteined from GDCquery function 
#' @param query A query form GDCquery function
#' @import EnsDb.Hsapiens.v79
#' @export 
#' @examples
#' Convert from ensembl.gene to gene.symbol
#' ensembl.genes <- c("ENSG00000150676", "ENSG00000099308", "ENSG00000142676", "ENSG00000180776", "ENSG00000108848", "ENSG00000277370", "ENSG00000103811", "ENSG00000101473")
#' symb<-Convert_toGeneSymbol(ensembl.genes)
#' @return a data.frame

Convert_toGeneSymbol <- function(ensembl.genes){
  
  symbol.gene <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  
 
  return(symbol.gene)
  
}




#' @title Convert from gene.symbol to ensembl.gene annotation
#' @description
#'  get information about query obteined from GDCquery function 
#' @param ensembl.genes list of character (character)
#' @import EnsDb.Hsapiens.v79
#' @export 
#' @examples

#' Convert from gene.symbol to ensembl.gene
#' symbol.gene <-  c('DDX26B','CCDC83',  'MAST3', 'RPL11', 'ZDHHC20',  'LUC7L3',  'SNORD49A',  'CTSH', 'ACOT8')
#' ens<-Convert_toGeneEnsembl(symbol.gene)
#' @return a data.frame

Convert_toGeneEnsembl <- function(symbol.genes){
  
  ensembl.genes <- ensembldb::select(EnsDb.Hsapiens.v79, keys= symbol.gene, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
  
  
  return(ensembl.genes)
  
}





#' @title Extract information from TCGA barcodes.
#' @description
#'    MMRFget_IDs allows user to extract metadata from identifiers The dataframe returned has columns for
#'  'project', 'patient','visit', 'source','marker'
#' @param data numeric matrix, each row represents a gene, each column represents a sample
#' @export
#' @return data frame with columns 'project', 'patient','visit', 'source','marker'
MMRFget_IDs <- function(dataMM) {
  IDs <- strsplit(c(colnames(dataMM)), "_")
  IDs <- plyr::ldply(IDs, rbind)
  colnames(IDs) <- c('project', 'patient','visit', 'source','marker')
  cols <- c("project", "patient", "visit",'marker')
#  IDs$patient <- apply(IDs[,cols],1,paste,collapse = "_" )
  barcode <- colnames(dataMM)
  IDs <- cbind(IDs, barcode)
  #condition <- gsub("11+[[:alpha:]]", "normal", as.character(IDs$sample))
 # condition  <- gsub("01+[[:alpha:]]", "cancer", condition)
  #IDs$condition <- condition
  IDs$myorder  <- 1:nrow(IDs)
  return(IDs)
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



