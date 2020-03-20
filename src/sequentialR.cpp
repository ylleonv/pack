#include "distribution.h"
#include "sequentialR.h"
using namespace std;
using namespace Rcpp ;

// [[Rcpp::depends(RcppEigen)]]

SequentialR::SequentialR(void) {
  Rcout << "SequentialR is being created" << endl;
}

Eigen::VectorXd SequentialR::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  double product = 1;
  for(size_t j=0; j<eta.size(); ++j)
  {
    ordered_pi[j] = product * Logistic::cdf_logit( eta(j) );
    product *= ( 1 - Logistic::cdf_logit( eta(j) ) );
  }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd SequentialR::inverse_derivative_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  double product = 1.;
  for (size_t j=0; j < eta.rows(); ++j)
  {
    M(j,j) = Logistic::pdf_logit(eta(j)) * product;
    for (size_t i=0; i<j; ++i)
    { M(i,j) = - Logistic::pdf_logit(eta(i))  * std::max(1e-10, std::min(Logistic::cdf_logit(eta(j)), 1-1e-6)) * product / std::max(1e-10, std::min( 1-Logistic::cdf_logit(eta(i)), 1-1e-6)); }
    product *= std::max(1e-10, std::min( 1-Logistic::cdf_logit(eta(j)), 1-1e-6));
  }
  return M;
}

Eigen::VectorXd SequentialR::inverse_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  double product = 1;
  for(size_t j=0; j<eta.size(); ++j)
  {
    ordered_pi[j] = product * cdf_normal( eta(j) );
    product *= ( 1 - cdf_normal( eta(j) ) );
  }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd SequentialR::inverse_derivative_normal(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  double product = 1.;
  for (size_t j=0; j < eta.rows(); ++j)
  {
    M(j,j) = pdf_normal(eta(j)) * product;
    for (size_t i=0; i<j; ++i)
    { M(i,j) = - pdf_normal(eta(i))  * std::max(1e-10, std::min(cdf_normal(eta(j)), 1-1e-6)) * product / std::max(1e-10, std::min( 1-cdf_normal(eta(i)), 1-1e-6)); }
    product *= std::max(1e-10, std::min( 1-cdf_normal(eta(j)), 1-1e-6));
  }
  return M;
}

Eigen::VectorXd SequentialR::inverse_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  double product = 1;
  for(size_t j=0; j<eta.size(); ++j)
  {
    ordered_pi[j] = product * cdf_cauchit( eta(j) );
    product *= ( 1 - cdf_cauchit( eta(j) ) );
  }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd SequentialR::inverse_derivative_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  double product = 1.;
  for (size_t j=0; j < eta.rows(); ++j)
  {
    M(j,j) = pdf_cauchit(eta(j)) * product;
    for (size_t i=0; i<j; ++i)
    { M(i,j) = - pdf_cauchit(eta(i))  * std::max(1e-10, std::min(cdf_cauchit(eta(j)), 1-1e-6)) * product / std::max(1e-10, std::min( 1-cdf_cauchit(eta(i)), 1-1e-6)); }
    product *= std::max(1e-10, std::min( 1-cdf_cauchit(eta(j)), 1-1e-6));
  }
  return M;
}

Eigen::VectorXd SequentialR::inverse_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  double product = 1;
  for(size_t j=0; j<eta.size(); ++j)
  {
    ordered_pi[j] = product * cdf_gompertz( eta(j) );
    product *= ( 1 - cdf_gompertz( eta(j) ) );
  }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd SequentialR::inverse_derivative_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  double product = 1.;
  for (size_t j=0; j < eta.rows(); ++j)
  {
    M(j,j) = pdf_gompertz(eta(j)) * product;
    for (size_t i=0; i<j; ++i)
    { M(i,j) = - pdf_gompertz(eta(i))  * std::max(1e-10, std::min(cdf_gompertz(eta(j)), 1-1e-6)) * product / std::max(1e-10, std::min( 1-cdf_gompertz(eta(i)), 1-1e-6)); }
    product *= std::max(1e-10, std::min( 1-cdf_gompertz(eta(j)), 1-1e-6));
  }
  return M;
}

List SequentialR::GLMseq(std::string response,
                                    StringVector explanatory_complete,
                                    StringVector explanatory_proportional,
                                    std::string distribution,
                                    SEXP categories_order,
                                    DataFrame dataframe){

  int P_c = 0;
  if(explanatory_complete[0] != "NA"){P_c = explanatory_complete.size(); }
  int P_p = 0;
  if(explanatory_proportional[0] != "NA"){P_p = explanatory_proportional.size(); }
  int P =  P_c +  P_p ; // Number of explanatory variables without intercept

  const int N = dataframe.nrows() ; // Number of observations

  List Full_M = distribution::select_data(dataframe, response, explanatory_complete,
                                          explanatory_proportional, categories_order);

  Eigen::MatrixXd Y_init = Full_M["Y_ext"];
  Eigen::MatrixXd X_EXT = Full_M["X_EXT"];
  // CharacterVector levs1 = Full_M["levs1"];

  int Q = Y_init.cols();
  int K = Q + 1;
  // // // Beta initialization with zeros
  Eigen::MatrixXd BETA;
  BETA = Eigen::MatrixXd::Zero(X_EXT.cols(),1);
  //
  int iteration = 0;
  // double check_tutz = 1.0;
  double Stop_criteria = 1.0;
  Eigen::MatrixXd X_M_i ;
  Eigen::VectorXd Y_M_i ;
  Eigen::VectorXd eta ;
  Eigen::VectorXd pi ;
  Eigen::MatrixXd D ;
  Eigen::MatrixXd Cov_i ;
  Eigen::MatrixXd W_in ;
  Eigen::MatrixXd Score_i_2 ;
  Eigen::MatrixXd F_i_2 ;
  Eigen::VectorXd LogLikIter;
  LogLikIter = Eigen::MatrixXd::Zero(1,1) ;
  double LogLik;
  // for (int iteration=1; iteration < 18; iteration++){
  // while (check_tutz > 0.0001){
  double epsilon = 0.0001 ;
  while (Stop_criteria >( epsilon / N)){

    Eigen::MatrixXd Score_i = Eigen::MatrixXd::Zero(BETA.rows(),1);
    Eigen::MatrixXd F_i = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());
    LogLik = 0.;

    // Loop by subject
    for (int i=0; i < N; i++){
      // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_init.row(i);
      eta = X_M_i * BETA;

      // Vector pi depends on selected distribution
      if(distribution == "logistic"){
        pi = SequentialR::inverse_logistic(eta);
        D = SequentialR::inverse_derivative_logistic(eta);
      }else if(distribution == "normal"){
        pi = SequentialR::inverse_normal(eta);
        D = SequentialR::inverse_derivative_normal(eta);
      }else if(distribution == "cauchit"){
        pi = SequentialR::inverse_cauchit(eta);
        D = SequentialR::inverse_derivative_cauchit(eta);
      }else if(distribution == "gompertz"){
        pi = SequentialR::inverse_gompertz(eta);
        D = SequentialR::inverse_derivative_gompertz(eta);
      }

      Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
      W_in = D * Cov_i.inverse();
      Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
      F_i = F_i + F_i_2;
      LogLik = LogLik + (Y_M_i.transpose().eval()*Eigen::VectorXd(pi.array().log())) + ( (1 - Y_M_i.sum()) * std::log(1 - pi.sum()) );
    }

    LogLikIter.conservativeResize(iteration+2, 1);
    LogLikIter(iteration+1) = LogLik;
    Stop_criteria = (abs(LogLikIter(iteration+1) - LogLikIter(iteration))) / (epsilon + (abs(LogLikIter(iteration+1)))) ;
    Eigen::VectorXd beta_old = BETA;
    BETA = BETA + (F_i.inverse() * Score_i);
    iteration = iteration + 1;

  }

  std::vector<std::string> text=as<std::vector<std::string>>(explanatory_complete);
  std::vector<std::string> level_text=as<std::vector<std::string>>(categories_order);
  StringVector names(Q*P_c + P_p);
  if(P_c > 0){
    for(int var = 0 ; var < explanatory_complete.size() ; var++){
      for(int cat = 0 ; cat < Q ; cat++){
        names[(Q*var) + cat] = distribution::concatenate(text[var], level_text[cat]);
      }
    }
  }
  if(P_p > 0){
    for(int var_p = 0 ; var_p < explanatory_proportional.size() ; var_p++){
      names[(Q*P_c) + var_p] = explanatory_proportional[var_p];
    }
  }

  // TO NAMED THE RESULT BETAS
  NumericMatrix BETA_2 = wrap(BETA);
  rownames(BETA_2) = names;

  return List::create(Named("Nb. iterations") = iteration-1 ,
                      Named("Coefficients") = BETA_2,
                      Named("Log-likelihood") = LogLik);
}

RCPP_MODULE(sequentialmodule){
  Rcpp::class_<SequentialR>("SequentialR")
  .constructor()
  .method( "GLMseq", &SequentialR::GLMseq )
  ;
}



