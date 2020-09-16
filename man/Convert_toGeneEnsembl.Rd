% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MMRFutility.R
\name{Convert_toGeneEnsembl}
\alias{Convert_toGeneEnsembl}
\title{Convert  from ensembl.gene to gene.symbol or gene.symbol to ensembl.gene}
\usage{
Convert_toGeneEnsembl(symbol.genes)
}
\arguments{
\item{ensembl.genes}{list of character (character)}
}
\description{
get information about query obteined from GDCquery function
}
\examples{
Convert from gene.symbol to ensembl.gene
symbol.gene <-  c('DDX26B','CCDC83',  'MAST3', 'RPL11', 'ZDHHC20',  'LUC7L3',  'SNORD49A',  'CTSH', 'ACOT8')
Convert_toGeneEnsembl(symbol.gene)
}