#' @title Prepare GDC data
#' @description
#'   Reads the data downloaded and prepare it into an R object
#' @param query A query for GDCquery function
#' @param save Save result as RData object?
#' @param save.filename Name of the file to be save if empty an automatic will be created
#' @param directory Directory/Folder where the data was downloaded. Default: GDCdata
#' @param summarizedExperiment Create a summarizedExperiment? Default TRUE (if possible)
#' @param remove.files.prepared Remove the files read? Default: FALSE
#' This argument will be considered only if save argument is set to true
#' @param add.gistic2.mut If a list of genes (gene symbol) is given, columns with gistic2 results from GDAC firehose (hg19)
#' and a column indicating if there is or not mutation in that gene (hg38)
#' (TRUE or FALSE - use the MAF file for more information)
#' will be added to the sample matrix in the summarized Experiment object.
#' @param mut.pipeline If add.gistic2.mut is not NULL this field will be taken in consideration.
#' Four separate variant calling pipelines are implemented for GDC data harmonization.
#' Options: muse, varscan2, somaticsniper, MuTect2. For more information:
#' https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/
#' @param mutant_variant_classification List of mutant_variant_classification that will be
#' consider a sample mutant or not. Default: "Frame_Shift_Del", "Frame_Shift_Ins",
#' "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "In_Frame_Del",
#' "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation"
#' @export
#' @examples
#' \dontrun{
#' query.mm.fpkm <- GDCquery(project = "MMRF-COMMPASS",
#'                           data.category = "Transcriptome Profiling",
#'                           data.type = "Gene Expression Quantification",
#'                           workflow.type="HTSeq - FPKM")
#'
#'
#' GDCdownload(query.mm.fpkm, method = "api", files.per.chunk = 100)
#' DataGDC.prep <- MMRFGDC_prepare(query.mm.fpkm,
#'                                 save = TRUE ,
#'                                 save.filename = "MMCompassFPKM.rda" ,
#'                                 directory = "GDCdata" ,
#'                                 summarizedExperiment = TRUE)
#' }
#' @return A summarizedExperiment or a data.frame
#' @importFrom  S4Vectors DataFrame
#' @importFrom SummarizedExperiment metadata<-
#' @importFrom data.table setcolorder setnames
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom TCGAbiolinks GDCprepare

MMRFGDC_prepare <- function(query,
                        save = FALSE,
                        save.filename,
                        directory = "GDCdata",
                        summarizedExperiment = TRUE,
                        remove.files.prepared = FALSE,
                        add.gistic2.mut = NULL,
                        mut.pipeline = "mutect2",
                        mutant_variant_classification = c("Frame_Shift_Del",
                                                          "Frame_Shift_Ins",
                                                          "Missense_Mutation",
                                                          "Nonsense_Mutation",
                                                          "Splice_Site",
                                                          "In_Frame_Del",
                                                          "In_Frame_Ins",
                                                          "Translation_Start_Site",
                                                          
                                                          "Nonstop_Mutation")){
  
  
  
                object<-GDCprepare(query,
                                   save,
                                   save.filename,
                                   directory,
                                   summarizedExperiment,
                                   remove.files.prepared,
                                   add.gistic2.mut,
                                   mut.pipeline,
                                   mutant_variant_classification)
  
  
  
  
  
  
  
  
  definition<-colData(object)$sample_type
  definition <- as.data.frame(definition)
  colData(object) <- cbind(colData(object), definition)
  colnames(object) <- substr(colnames(object),1,9)
  return(object)
  
}




#' @title Prepare GDC data (extended identifiers)
#' @description
#'   Reads the data downloaded and prepare it into an R object
#' @param query A query for GDCquery function
#' @param save Save result as RData object?
#' @param save.filename Name of the file to be save if empty an automatic will be created
#' @param directory Directory/Folder where the data was downloaded. Default: GDCdata
#' @param summarizedExperiment Create a summarizedExperiment? Default TRUE (if possible)
#' @param remove.files.prepared Remove the files read? Default: FALSE
#' This argument will be considered only if save argument is set to true
#' @param add.gistic2.mut If a list of genes (gene symbol) is given, columns with gistic2 results from GDAC firehose (hg19)
#' and a column indicating if there is or not mutation in that gene (hg38)
#' (TRUE or FALSE - use the MAF file for more information)
#' will be added to the sample matrix in the summarized Experiment object.
#' @param mut.pipeline If add.gistic2.mut is not NULL this field will be taken in consideration.
#' Four separate variant calling pipelines are implemented for GDC data harmonization.
#' Options: muse, varscan2, somaticsniper, MuTect2. For more information:
#' https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/
#' @param mutant_variant_classification List of mutant_variant_classification that will be
#' consider a sample mutant or not. Default: "Frame_Shift_Del", "Frame_Shift_Ins",
#' "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "In_Frame_Del",
#' "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation"
#' @export
#' @examples
#' \dontrun{
#' query.mm.fpkm <- GDCquery(project = "MMRF-COMMPASS",
#'                           data.category = "Transcriptome Profiling",
#'                           data.type = "Gene Expression Quantification",
#'                           workflow.type="HTSeq - FPKM")
#'
#'
#' GDCdownload(query.mm.fpkm, method = "api", files.per.chunk = 100)
#' DataGDC.prep <- MMRFGDC_prepare_extended(query.mm.fpkm,
#'                                          save = TRUE ,
#'                                          save.filename = "MMCompassFPKM.rda" ,
#'                                          directory = "GDCdata" ,
#'                                          summarizedExperiment = TRUE)
#' }
#' @return A summarizedExperiment or a data.frame (with complete identifiers)
#' @importFrom  S4Vectors DataFrame
#' @importFrom SummarizedExperiment metadata<-
#' @importFrom data.table setcolorder setnames
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom TCGAbiolinks GDCprepare

MMRFGDC_prepare_extended <- function(query,
                            save = FALSE,
                            save.filename,
                            directory = "GDCdata",
                            summarizedExperiment = TRUE,
                            remove.files.prepared = FALSE,
                            add.gistic2.mut = NULL,
                            mut.pipeline = "mutect2",
                            mutant_variant_classification = c("Frame_Shift_Del",
                                                              "Frame_Shift_Ins",
                                                              "Missense_Mutation",
                                                              "Nonsense_Mutation",
                                                              "Splice_Site",
                                                              "In_Frame_Del",
                                                              "In_Frame_Ins",
                                                              "Translation_Start_Site",
                                                              
                                                              "Nonstop_Mutation")){
  
  
  
  object<-GDCprepare(query,
                     save,
                     save.filename,
                     directory,
                     summarizedExperiment,
                     remove.files.prepared,
                     add.gistic2.mut,
                     mut.pipeline,
                     mutant_variant_classification)
  
  
  
  
  
  
  
  
  definition<-colData(object)$sample_type
  definition <- as.data.frame(definition)
  colData(object) <- cbind(colData(object), definition)
#  colnames(object) <- substr(colnames(object),1,9)
  return(object)
  
}





