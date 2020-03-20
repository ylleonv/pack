#include "distribution.h"
#include "adjacentR.h"
using namespace std;
using namespace Rcpp ;

// [[Rcpp::depends(RcppEigen)]]

AdjacentR::AdjacentR(void) {
  cout << "Adjacent is being created" << endl;
}

Eigen::VectorXd AdjacentR::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_logit( eta(eta.size()-1) ) / ( 1-cdf_logit( eta(eta.size()-1) ) );
  double norm = 1 + pi[eta.size()-1];
  for(size_t j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_logit( eta(j-1) ) / ( 1-cdf_logit( eta(j-1) ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_logistic(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(size_t j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_logit( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_logit(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_logit(eta(j)))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}

Eigen::VectorXd AdjacentR::inverse_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_normal( eta(eta.size()-1) ) / ( 1-cdf_normal( eta(eta.size()-1) ) );
  double norm = 1 + pi[eta.size()-1];
  for(size_t j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_normal( eta(j-1) ) / ( 1-cdf_normal( eta(j-1) ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_normal(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(size_t j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_normal( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_normal(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_normal(eta(j)))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}

Eigen::VectorXd AdjacentR::inverse_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_cauchit( eta(eta.size()-1) ) / ( 1-cdf_cauchit( eta(eta.size()-1) ) );
  double norm = 1 + pi[eta.size()-1];
  for(size_t j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_cauchit( eta(j-1) ) / ( 1-cdf_cauchit( eta(j-1) ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_cauchit(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(size_t j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_cauchit( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_cauchit(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_cauchit(eta(j)))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}

Eigen::VectorXd AdjacentR::inverse_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_gompertz( eta(eta.size()-1) ) / ( 1-cdf_gompertz( eta(eta.size()-1) ) );
  double norm = 1 + pi[eta.size()-1];
  for(size_t j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_gompertz( eta(j-1) ) / ( 1-cdf_gompertz( eta(j-1) ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_gompertz(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(size_t j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_gompertz( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_gompertz(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_gompertz(eta(j)))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}

Eigen::VectorXd AdjacentR::inverse_gumbel(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_gumbel( eta(eta.size()-1) ) / ( 1-cdf_gumbel( eta(eta.size()-1) ) );
  double norm = 1 + pi[eta.size()-1];
  for(size_t j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_gumbel( eta(j-1) ) / ( 1-cdf_gumbel( eta(j-1) ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_gumbel(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_gumbel(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(size_t j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_gumbel( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_gumbel(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_gumbel(eta(j)))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}

List AdjacentR::GLMadj(std::string response,
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
        pi = AdjacentR::inverse_logistic(eta);
        D = AdjacentR::inverse_derivative_logistic(eta);
      }else if(distribution == "normal"){
        pi = AdjacentR::inverse_normal(eta);
        D = AdjacentR::inverse_derivative_normal(eta);
      }else if(distribution == "cauchit"){
        pi = AdjacentR::inverse_cauchit(eta);
        D = AdjacentR::inverse_derivative_cauchit(eta);
      }else if(distribution == "gompertz"){
        pi = AdjacentR::inverse_gompertz(eta);
        D = AdjacentR::inverse_derivative_gompertz(eta);
      }else if(distribution == "gumbel"){
        pi = AdjacentR::inverse_gumbel(eta);
        D = AdjacentR::inverse_derivative_gumbel(eta);
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

RCPP_MODULE(adjacentmodule){
  Rcpp::class_<AdjacentR>("AdjacentR")
  .constructor()
  .method( "GLMadj", &AdjacentR::GLMadj)
  ;
}

