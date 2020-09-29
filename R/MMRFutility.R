
#' @title Get MMRF-COMMPASS Project Summary
#' @description
#'  get information about MMRF-COMMPASS Project form GDC Data Portal 
#' @import jsonlite
#' @export 
#' @examples
#' \dontrun{
#' summary<-MMRFgetProjectSummary()
#' summary$data_categories
#' summary$experimental_strategies
#' 
#' 
#' }
#' @return a list

MMRFprojectGDC_Summary <- function(){
  project<-"MMRF-COMMPASS"
  baseURL <- "https://api.gdc.cancer.gov/projects/"
  url <- paste0(baseURL, project,"?expand=summary,summary.experimental_strategies,summary.data_categories&pretty=true")
  return(fromJSON(url,simplifyDataFrame = TRUE)$data$summary)
}


#' @title Get GDC Query Summary
#' @description
#'  get information about query obteined from GDCquery function 
#' @param query A query form GDCquery function
#' @import dplyr
#' @export 
#' @examples
#' \dontrun{
#' query.mm.fpkm <- GDCquery(project = "MMRF-COMMPASS",
#'                           data.category = "Transcriptome Profiling",
#'                           data.type = "Gene Expression Quantification",
#'                           workflow.type="HTSeq - FPKM")
#' summary<-MMRFqueryGDC_Summary(query.mm.fpkm)
#' 
#' 
#' }
#' @return a data.frame


MMRFqueryGDC_Summary <- function(query){
  
 res.query<-getResults(query)
  
  if(query$data.category=="Transcriptome Profiling"){
    
    query.group<-group_by(res.query,data_type,experimental_strategy,sample_type,analysis_workflow_type)
    summary <- summarize(query.group, numb = n())
    
    } else if(query$data.category=="Simple Nucleotide Variation") {
    
    query.group<-group_by(res.query,experimental_strategy,analysis_workflow_type)
    summary <- summarize(query.group, n_cases = n())
    
  } else if(query$data.category=="Sequencing Reads") {
   
   query.group<-group_by(res.query,data_type,experimental_strategy,sample_type,analysis_workflow_type)
   summary <- summarize(query.group, numb = n())
   
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
#' Convert_toGeneSymbol(ensembl.genes)
#' 
#' 
#' 

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
#' Convert_toGeneEnsembl(symbol.gene)
 

Convert_toGeneEnsembl <- function(symbol.genes){
  
  ensembl.genes <- ensembldb::select(EnsDb.Hsapiens.v79, keys= symbol.gene, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
  
  
  return(ensembl.genes)
  
}






