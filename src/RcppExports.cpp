// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// adj_fun
List adj_fun(std::string response, StringVector explanatory_complete, StringVector explanatory_proportional, std::string distribution, SEXP categories_order, DataFrame dataframe);
RcppExport SEXP _pack_adj_fun(SEXP responseSEXP, SEXP explanatory_completeSEXP, SEXP explanatory_proportionalSEXP, SEXP distributionSEXP, SEXP categories_orderSEXP, SEXP dataframeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type response(responseSEXP);
    Rcpp::traits::input_parameter< StringVector >::type explanatory_complete(explanatory_completeSEXP);
    Rcpp::traits::input_parameter< StringVector >::type explanatory_proportional(explanatory_proportionalSEXP);
    Rcpp::traits::input_parameter< std::string >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< SEXP >::type categories_order(categories_orderSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type dataframe(dataframeSEXP);
    rcpp_result_gen = Rcpp::wrap(adj_fun(response, explanatory_complete, explanatory_proportional, distribution, categories_order, dataframe));
    return rcpp_result_gen;
END_RCPP
}
// GLMcum
List GLMcum(std::string response, StringVector explanatory_complete, StringVector explanatory_proportional, std::string distribution, SEXP categories_order, DataFrame dataframe, StringVector beta_t, Eigen::VectorXd beta_init);
RcppExport SEXP _pack_GLMcum(SEXP responseSEXP, SEXP explanatory_completeSEXP, SEXP explanatory_proportionalSEXP, SEXP distributionSEXP, SEXP categories_orderSEXP, SEXP dataframeSEXP, SEXP beta_tSEXP, SEXP beta_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type response(responseSEXP);
    Rcpp::traits::input_parameter< StringVector >::type explanatory_complete(explanatory_completeSEXP);
    Rcpp::traits::input_parameter< StringVector >::type explanatory_proportional(explanatory_proportionalSEXP);
    Rcpp::traits::input_parameter< std::string >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< SEXP >::type categories_order(categories_orderSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type dataframe(dataframeSEXP);
    Rcpp::traits::input_parameter< StringVector >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta_init(beta_initSEXP);
    rcpp_result_gen = Rcpp::wrap(GLMcum(response, explanatory_complete, explanatory_proportional, distribution, categories_order, dataframe, beta_t, beta_init));
    return rcpp_result_gen;
END_RCPP
}
// GLMref
List GLMref(std::string response, StringVector explanatory_complete, StringVector explanatory_proportional, std::string distribution, SEXP categories_order, DataFrame dataframe, double freedom_degrees);
RcppExport SEXP _pack_GLMref(SEXP responseSEXP, SEXP explanatory_completeSEXP, SEXP explanatory_proportionalSEXP, SEXP distributionSEXP, SEXP categories_orderSEXP, SEXP dataframeSEXP, SEXP freedom_degreesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type response(responseSEXP);
    Rcpp::traits::input_parameter< StringVector >::type explanatory_complete(explanatory_completeSEXP);
    Rcpp::traits::input_parameter< StringVector >::type explanatory_proportional(explanatory_proportionalSEXP);
    Rcpp::traits::input_parameter< std::string >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< SEXP >::type categories_order(categories_orderSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type dataframe(dataframeSEXP);
    Rcpp::traits::input_parameter< double >::type freedom_degrees(freedom_degreesSEXP);
    rcpp_result_gen = Rcpp::wrap(GLMref(response, explanatory_complete, explanatory_proportional, distribution, categories_order, dataframe, freedom_degrees));
    return rcpp_result_gen;
END_RCPP
}
// GLMseq
List GLMseq(std::string response, StringVector explanatory_complete, StringVector explanatory_proportional, std::string distribution, SEXP categories_order, DataFrame dataframe);
RcppExport SEXP _pack_GLMseq(SEXP responseSEXP, SEXP explanatory_completeSEXP, SEXP explanatory_proportionalSEXP, SEXP distributionSEXP, SEXP categories_orderSEXP, SEXP dataframeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type response(responseSEXP);
    Rcpp::traits::input_parameter< StringVector >::type explanatory_complete(explanatory_completeSEXP);
    Rcpp::traits::input_parameter< StringVector >::type explanatory_proportional(explanatory_proportionalSEXP);
    Rcpp::traits::input_parameter< std::string >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< SEXP >::type categories_order(categories_orderSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type dataframe(dataframeSEXP);
    rcpp_result_gen = Rcpp::wrap(GLMseq(response, explanatory_complete, explanatory_proportional, distribution, categories_order, dataframe));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_adjacentmodule();
RcppExport SEXP _rcpp_module_boot_fishder();
RcppExport SEXP _rcpp_module_boot_cumulativemodule();
RcppExport SEXP _rcpp_module_boot_exportmod();
RcppExport SEXP _rcpp_module_boot_sequentialmodule();

static const R_CallMethodDef CallEntries[] = {
    {"_pack_adj_fun", (DL_FUNC) &_pack_adj_fun, 6},
    {"_pack_GLMcum", (DL_FUNC) &_pack_GLMcum, 8},
    {"_pack_GLMref", (DL_FUNC) &_pack_GLMref, 7},
    {"_pack_GLMseq", (DL_FUNC) &_pack_GLMseq, 6},
    {"_rcpp_module_boot_adjacentmodule", (DL_FUNC) &_rcpp_module_boot_adjacentmodule, 0},
    {"_rcpp_module_boot_fishder", (DL_FUNC) &_rcpp_module_boot_fishder, 0},
    {"_rcpp_module_boot_cumulativemodule", (DL_FUNC) &_rcpp_module_boot_cumulativemodule, 0},
    {"_rcpp_module_boot_exportmod", (DL_FUNC) &_rcpp_module_boot_exportmod, 0},
    {"_rcpp_module_boot_sequentialmodule", (DL_FUNC) &_rcpp_module_boot_sequentialmodule, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_pack(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
