#include "distribution.h"
#include "reference.h"
using namespace std;
using namespace Rcpp ;

// [[Rcpp::depends(RcppEigen)]]

ReferenceF::ReferenceF(void) {
  Rcout << "Ref" << endl;
}

Eigen::VectorXd ReferenceF::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(size_t j=0; j<eta.size(); ++j)
  {
    pi[j] = Logistic::cdf_logit( eta(j) ) / ( 1-Logistic::cdf_logit( eta(j) ) );
    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::VectorXd ReferenceF::inverse_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(size_t j=0; j<eta.size(); ++j)
  {
    pi[j] = cdf_normal( eta(j) ) / ( 1-cdf_normal( eta(j) ) );
    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_logistic(const Eigen::VectorXd& eta2) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_logistic(eta2);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  { D1(j,j) = pdf_logit( eta2(j) ) /
    (Logistic::cdf_logit(eta2(j)) * (1-Logistic::cdf_logit(eta2(j))));
  }
  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}

Eigen::MatrixXd ReferenceF::inverse_derivative_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = ReferenceF::inverse_normal(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  for(size_t j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_normal( eta(j) ) /
    ( std::max(1e-10, std::min(1-1e-6,cdf_normal(eta(j)))) *
      std::max(1e-10, std::min(1-1e-6, 1-cdf_normal(eta(j)))) ); }
  return D * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose().eval() );
}

Eigen::VectorXd ReferenceF::inverse_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(size_t j=0; j<eta.size(); ++j)
  {
    pi[j] = Cauchit::cdf_cauchit( eta(j) ) / ( 1-Cauchit::cdf_cauchit( eta(j) ) );
    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_cauchit(const Eigen::VectorXd& eta2) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_cauchit(eta2);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  { D1(j,j) = pdf_cauchit( eta2(j) ) /
    (Cauchit::cdf_cauchit(eta2(j)) * (1-Cauchit::cdf_cauchit(eta2(j))));
  }
  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}


List ReferenceF::GLMref(std::string response,
                        StringVector explanatory_complete,
                        StringVector explanatory_proportional,
                        std::string distribution,
                        SEXP categories_order,
                        DataFrame dataframe){

  int P_c = 0;
  if(explanatory_complete[0] != "NA"){P_c = explanatory_complete.size(); }
  int P_p = 0;
  if(explanatory_proportional[0] == "NA"){P_p = explanatory_proportional.size(); }
  int P =  P_c +  P_p ; // Number of explanatory variables without intercept
  const int N = dataframe.nrows() ; // Number of observations

  List Full_M = distribution::select_data(dataframe, response, explanatory_complete,
                                          explanatory_proportional, categories_order);

  Eigen::MatrixXd Y_init = Full_M["Y_ext"];
  Eigen::MatrixXd X_EXT = Full_M["X_EXT"];
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
        pi = ReferenceF::inverse_logistic(eta);
        D = ReferenceF::inverse_derivative_logistic(eta);
      }else if(distribution == "normal"){
        pi = ReferenceF::inverse_normal(eta);
        D = ReferenceF::inverse_derivative_normal(eta);
      }else if(distribution == "cauchit"){
        pi = ReferenceF::inverse_cauchit(eta);
        D = ReferenceF::inverse_derivative_cauchit(eta);
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
    // check_tutz = ((BETA - beta_old).norm())/(beta_old.norm()+check_tutz);
    iteration = iteration + 1;
    // Rcout << "Beta" << endl;
    // Rcout << BETA << endl;
  }
  Rcout << "Number of iterations" << endl;
  Rcout << iteration-1 << endl;
  Rcout << "Log Likelihood" << endl;
  Rcout << LogLik << endl;
  // Rcout << "Beta" << endl;
  // Rcout << BETA << endl;
  return List::create(Named("Nb. iterations") = iteration-1 , Named("Coefficients") = BETA,
                      Named("Log-likelihood") = LogLik);


}

List ReferenceF::GLMref_ec(std::string response, std::string actual_response,
                           std::string individuals,
                           StringVector explanatory_complete,
                           StringVector depend_y,
                           std::string distribution,
                           SEXP categories_order,
                           DataFrame dataframe){
  int P_c = 0;
  if(explanatory_complete[0] != "NA"){P_c = explanatory_complete.size(); }
  int P_p = 0;
  const int N = dataframe.nrows() ; // Number of observations

  List Full_M = distribution::select_data_nested(dataframe, response, actual_response,
                                                 individuals, explanatory_complete, depend_y , categories_order);

  Eigen::MatrixXd Y_init = Full_M["Y_ext"];
  Eigen::MatrixXd X_EXT = Full_M["X_M_dep_y"];

  int Q = Y_init.cols();
  int K = Q + 1;
  // // // Beta initialization with zeros
  Eigen::MatrixXd BETA;
  BETA = Eigen::MatrixXd::Zero(X_EXT.cols(),1);

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
    for (int i=0; i < N/K; i++){
      // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_init.row(i);
      eta = X_M_i * BETA;

      // Vector pi depends on selected distribution
      if(distribution == "logistic"){
        pi = ReferenceF::inverse_logistic(eta);
        D = ReferenceF::inverse_derivative_logistic(eta);
      }else if(distribution == "normal"){
        pi = ReferenceF::inverse_normal(eta);
        D = ReferenceF::inverse_derivative_normal(eta);
      }else if(distribution == "cauchit"){
        pi = ReferenceF::inverse_cauchit(eta);
        D = ReferenceF::inverse_derivative_cauchit(eta);
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
    // check_tutz = ((BETA - beta_old).norm())/(beta_old.norm()+check_tutz);
    iteration = iteration + 1;

  }
  Rcout << "Number of iterations" << endl;
  Rcout << iteration-1 << endl;
  Rcout << "Log Likelihood" << endl;
  Rcout << LogLik << endl;
  // Rcout << "Beta" << endl;
  // Rcout << BETA << endl;
  return List::create(Named("Nb. iterations") = iteration-1 , Named("Coefficients") = BETA,
                      Named("Log-likelihood") = LogLik);
  // return List::create(Named("Coefficients") = BETA,
  //                     Named("Y_init") = Y_init,
  //                     Named("X_EXT") = X_EXT);

}

// [[Rcpp::export]]
RCPP_MODULE(referencemodule){
  Rcpp::class_<ReferenceF>("ReferenceF")
  .constructor()
  .method( "GLMref", &ReferenceF::GLMref )
  .method( "GLMref_ec", &ReferenceF::GLMref_ec )
  ;
}
