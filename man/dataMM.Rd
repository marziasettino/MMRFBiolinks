\name{dataMM}
\alias{dataMM}
\docType{data}
\title{
dataMM is a matrix of Gene expression  from MMRFGDC_prepare 
}
\description{
dataMM is a matrix of Gene expression  from MMRFGDC_prepare 
}
\usage{dataMM}
\format{
  A data frame with genes in rows, samples in cols
}

\examples{
listSamples <- c("MMRF_2473","MMRF_2111",
                 "MMRF_2362","MMRF_1824",
                 "MMRF_1458","MRF_1361",
                 "MMRF_2203","MMRF_2762",
                 "MMRF_2680","MMRF_1797")
                 
query <- GDCquery(project = "MMRF-COMMPASS", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type="HTSeq - FPKM",
                  barcode = listSamples)


GDCdownload(query, method = "api", files.per.chunk = 10)
clin.mm<-MMRFqueryGDC_clinic(type = "clinical")

dataPrepared<-MMRFGDC_prepare(query,
                              save = TRUE,
                              save.filename = "MMCompassFPKM.rda" ,
                              directory = "GDCdata",
                              summarizedExperiment = TRUE)
                               
dataMM<-assay(dataPrepared)                         


}
\keyword{datasets}
