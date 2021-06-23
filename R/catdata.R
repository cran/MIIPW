#' @title Categorical repeated measurement data
#'
#' @description dataset of observations made at various timepoints, with variable age,group
#' ,gender. there are missing values in variables age,group also.
#' 
#' @usage data(catdata)
#' @format A \code{tibble} with 18 columns which are :
#' \describe{
#' \item{age}{age of subject}
#' \item{grp}{group}
#' \item{gen}{gender}
#' \item{0}{Observation on timepoint 0}
#' \item{1}{Observation on timepoint 1}
#' \item{2}{Observation on timepoint 2}
#' \item{3}{Observation on timepoint 3}
#' \item{4}{Observation on timepoint 4}
#' \item{2}{Observation on timepoint 2}
#' \item{7}{Observation on timepoint 7}
#' \item{14}{Observation on timepoint 14}
#' \item{agefull}{age of subject,no missing value }
#' \item{grpfull}{group,no missing value }
#' \item{genfull}{gender, no missing value}
#' \item{X1}{takes value 0,1}
#' \item{X2}{takes value 0,1}
#' \item{X3}{takes value 0,1}
#' \item{X4}{takes value 0,1}}
#'
#'
"catdata"