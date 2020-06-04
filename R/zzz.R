## usethis namespace: start
#' @import Rcpp
#' @export GLMadj
#' @export GLMcum
#' @export GLMseq
#' @export GLMref
#' @export Discrete_CM
#' @export Predict_Response
#' @export summary.pcglm
#' @export ReferenceF
#' @useDynLib pack, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end

loadModule("exportmod", TRUE)
loadModule("cumulativemodule", TRUE)
loadModule("sequentialmodule", TRUE)
loadModule("adjacentmodule", TRUE)
loadModule("referencemodule", TRUE)
# loadModule("fishder", TRUE)
