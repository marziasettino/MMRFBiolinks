
#' @title Get MMRF-COMPASS Project Summary from GDC
#' @param project is a  GDC project
#' @export
#' @examples
#' \dontrun{
#' getProjectSummary(project="MMRF-COMMPASS")
#' }
MMRF_getProjectSummary <- function(project="MMRF-COMMPASS"){
  
  baseURL <- "https://api.gdc.cancer.gov/projects/"
  url <- paste0(baseURL, project,"?expand=summary,summary.data_categories&pretty=true")
  summary<-fromJSON(url,simplifyDataFrame = TRUE)$data$summary
  return(summary)
}




