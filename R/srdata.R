#' @title protein data
#' @description Repeated measurement dataset, for each id we have four visit observations
#' @usage data(srdata)
#' @format A dataframe with 164 rows and 30 columns
#' \describe{
#' \item{ID}{ID of subjects}
#' \item{Visit}{Number of times observations recorded}
#' \item{event}{death as event 1 if died or 0 if alive}
#' \item{OS}{Duration of overall survival}
#' \item{leftcensored}{Left censoring information}
#' \item{lc}{Left censoring information}
#' \item{C6kine,.....,GFRalpha4}{These are covariates}}
#' @examples data(srdata)
"srdata"

