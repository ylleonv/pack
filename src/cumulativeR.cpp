#include "cumulativeR.h"
using namespace std;
using namespace Rcpp ;

#include <algorithm>    // std::sort

CumulativeR::CumulativeR(void) {
  Rcout << "CumulativeR is being created" << endl;
}

Eigen::VectorXd CumulativeR::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_logit( eta(0) );
  for(size_t j=1; j<eta.size(); ++j)
  { ordered_pi[j] = Logistic::cdf_logit( eta(j) ) -  Logistic::cdf_logit( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::VectorXd CumulativeR::inverse_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_normal( eta(0) );
  for(size_t j=1; j<eta.size(); ++j)
  { ordered_pi[j] = cdf_normal( eta(j) ) - cdf_normal( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(size_t j=0; j<eta.rows(); ++j)
  {F(j,j) =  Logistic::pdf_logit(eta(j));}
  return (F * R);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_normal(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(size_t j=0; j<eta.rows(); ++j)
  { F(j,j) = pdf_normal( eta(j) ); }
  return (F * R);
}

Eigen::VectorXd CumulativeR::inverse_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_cauchit( eta(0) );
  for(size_t j=1; j<eta.size(); ++j)
  { ordered_pi[j] = cdf_cauchit( eta(j) ) - cdf_cauchit( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(size_t j=0; j<eta.rows(); ++j)
  { F(j,j) = cdf_cauchit( eta(j) ); }
  return (F * R);
}

Eigen::VectorXd CumulativeR::inverse_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_gompertz( eta(0) );
  for(size_t j=1; j<eta.size(); ++j)
  { ordered_pi[j] = cdf_gompertz( eta(j) ) - cdf_gompertz( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(size_t j=0; j<eta.rows(); ++j)
  { F(j,j) = cdf_gompertz( eta(j) ); }
  return (F * R);
}

List CumulativeR::GLMcum(std::string response,
                                    StringVector explanatory_complete,
                                    StringVector explanatory_proportional,
                                    std::string distribution,
                                    SEXP categories_order,
                                    DataFrame dataframe,
                                    StringVector beta_t,
                                    Eigen::VectorXd beta_init){

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
  Eigen::MatrixXd BETA2;
  BETA2 = Eigen::MatrixXd::Zero(X_EXT.cols(),1);

  Eigen::VectorXd BETA3 = BETA2;

  NumericVector seq_cat2;

  // // Beta initialization modification
  for (int ex_num = 0 ; ex_num < explanatory_complete.size() ; ex_num++ ){
    if (explanatory_complete[ex_num] == "intercept"){
      IntegerVector seq_cat = seq_len(Q);
      NumericVector seq_cat1 = as<NumericVector>(seq_cat);
      seq_cat2 = seq_cat1 / K;
      Rcout << seq_cat2 << std::endl;
      Eigen::Map<Eigen::VectorXd> seq_cat3 = as<Eigen::Map<Eigen::VectorXd> >(seq_cat2);
      Eigen::VectorXd beta_ini_1;

      if(distribution == "logistic"){
        beta_ini_1 = Logistic::InverseLinkQuantileFunction(seq_cat3);
      }else if (distribution == "normal"){
        beta_ini_1 = Normal::InverseLinkQuantileFunction(seq_cat3);
      }else if(distribution == "cauchit"){
        beta_ini_1 = Cauchit::InverseLinkQuantileFunction(seq_cat3);
      }else if(distribution == "gompertz"){
        beta_ini_1 = Gompertz::InverseLinkQuantileFunction(seq_cat3);
      }

      for(int el_be = 0 ; el_be < Q ; el_be++){
        BETA3[el_be] = beta_ini_1[el_be];
      }
    }
  }



  Eigen::MatrixXd BETA = BETA3;

  if(beta_t[0] == "TRUE"){
    BETA = beta_init;
  }

  Rcout << "Beta inizialization" << std::endl;

  Rcout << BETA << std::endl;

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
        pi = CumulativeR::inverse_logistic(eta);
        D = CumulativeR::inverse_derivative_logistic(eta);
      }else if(distribution == "normal"){
        pi = CumulativeR::inverse_normal(eta);
        D = CumulativeR::inverse_derivative_normal(eta);
      }else if(distribution == "cauchit"){
        pi = CumulativeR::inverse_cauchit(eta);
        D = CumulativeR::inverse_derivative_cauchit(eta);
      }else if(distribution == "gompertz"){
        pi = CumulativeR::inverse_gompertz(eta);
        D = CumulativeR::inverse_derivative_gompertz(eta);
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

RCPP_MODULE(cumulativemodule){
  Rcpp::class_<CumulativeR>("CumulativeR")
  .constructor()
  .method( "GLMcum", &CumulativeR::GLMcum )
  ;
}




