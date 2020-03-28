## usethis namespace: start
#' @import Rcpp
#' @export adj_fun
#' @export CumulativeR
#' @export SequentialR
#' @useDynLib pack, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end

loadModule("exportmod", TRUE)
loadModule("cumulativemodule", TRUE)
loadModule("sequentialmodule", TRUE)
loadModule("adjacentmodule", TRUE)
# loadModule("referencemodule", TRUE)
# loadModule("fishder", TRUE)
