// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP _rcpp_module_boot_adjacentmodule();
RcppExport SEXP _rcpp_module_boot_fishder();
RcppExport SEXP _rcpp_module_boot_cumulativemodule();
RcppExport SEXP _rcpp_module_boot_exportmod();
RcppExport SEXP _rcpp_module_boot_referencemodule();
RcppExport SEXP _rcpp_module_boot_sequentialmodule();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_adjacentmodule", (DL_FUNC) &_rcpp_module_boot_adjacentmodule, 0},
    {"_rcpp_module_boot_fishder", (DL_FUNC) &_rcpp_module_boot_fishder, 0},
    {"_rcpp_module_boot_cumulativemodule", (DL_FUNC) &_rcpp_module_boot_cumulativemodule, 0},
    {"_rcpp_module_boot_exportmod", (DL_FUNC) &_rcpp_module_boot_exportmod, 0},
    {"_rcpp_module_boot_referencemodule", (DL_FUNC) &_rcpp_module_boot_referencemodule, 0},
    {"_rcpp_module_boot_sequentialmodule", (DL_FUNC) &_rcpp_module_boot_sequentialmodule, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_pack(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
