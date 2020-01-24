## usethis namespace: start
#' @import Rcpp
#' @export FisherScoring
#' @export ReferenceF
#' @useDynLib pack, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end

loadModule("exportmod", TRUE)
loadModule("exportmoddev", TRUE)
loadModule("fishder", TRUE)
loadModule("referencemodule", TRUE)
