## usethis namespace: start
#' @import Rcpp
#' @export FisherScoring
#' @export ReferenceF
#' @export CumulativeR
#' @export AdjacentR
#' @export SequentialR
#' @useDynLib pack, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end

loadModule("exportmod", TRUE)
loadModule("fishder", TRUE)
loadModule("referencemodule", TRUE)
loadModule("cumulativemodule", TRUE)
loadModule("sequentialmodule", TRUE)
loadModule("adjacentmodule", TRUE)
