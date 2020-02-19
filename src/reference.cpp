#include "distribution.h"
#include "reference.h"
using namespace std;
using namespace Rcpp ;

#include <algorithm>    // std::sort

ReferenceF::ReferenceF(void) {
  Rcout << "Reference is being created" << endl;
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

Eigen::VectorXd ReferenceF::inverse_probit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(size_t j=0; j<eta.size(); ++j)
  {
    pi[j] = Probit::cdf_probit( eta(j) ) / ( 1-Probit::cdf_probit( eta(j) ) );
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

Eigen::MatrixXd ReferenceF::inverse_derivative_probit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = ReferenceF::inverse_probit(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  for(size_t j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_probit( eta(j) ) /
    ( std::max(1e-10, std::min(1-1e-6,Probit::cdf_probit(eta(j)))) *
      std::max(1e-10, std::min(1-1e-6, 1-Probit::cdf_probit(eta(j)))) ); }
  return D * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose().eval() );
}


Eigen::MatrixXd ReferenceF::GLMref(std::string response,
                                   StringVector explanatory_complete,
                                   StringVector explanatory_proportional,
                                   std::string distribution,
                                   NumericVector categories_order,
                                   DataFrame dataframe){
  // Number of explanatory variables
  const int P_c = explanatory_complete.size() ;
  const int P_p = explanatory_proportional.size() ;
  const int P =  P_c +  P_p ;
  const int N = dataframe.nrows() ; // Number of observations
  // Add Intercept acording to user input
  Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(N);
  dataframe["intercept"] = Ones1;

  Eigen::MatrixXd Full_M = distribution::select_data(dataframe, response, explanatory_complete,
                                                     explanatory_proportional, categories_order);
  // Order by Y
  Eigen::MatrixXd Full_M_ordered = distribution::sorted_rows(Full_M);

  // Count number of unique categories (K) in Y
  std::vector<int> Unique_k(Full_M_ordered.col(0).data(), Full_M_ordered.data() + Full_M_ordered.rows());
  int K = std::set<int>( Unique_k.begin(), Unique_k.end() ).size();
  int Q = K-1 ;
  Eigen::VectorXd Y_ordered = Full_M_ordered.col(0);

  // Expand Y to matrix where each column is a dummy variable indicating the category at which the subject row belongs:
  // Matrix of zeros of size NxK
  Eigen::MatrixXd Y_init2 = Eigen::MatrixXd::Zero(Full_M_ordered.rows(), K);
  // Replace 0 by 1 at specific column
  for (int i_1 = 0; i_1 < Y_init2.rows(); ++i_1){
    Y_init2(i_1,Y_ordered[i_1]) = 1;
  }
  // Eliminate last column (Reference level will be always last category)
  Eigen::MatrixXd Y_init = Y_init2.leftCols(Q);

  // Create design matrix
  Eigen::MatrixXd X_init = Full_M_ordered.rightCols(P);
  Eigen::MatrixXd X_EXT(2, 2);

  Eigen::MatrixXd X_M_Complete_Ext(N*Q, P_c) ;
  X_M_Complete_Ext << X_init.block(0,0, N , P_c);
  // Expand matrix
  X_M_Complete_Ext = kroneckerProduct(X_M_Complete_Ext,Eigen::MatrixXd::Identity(Q,Q)).eval();

  Eigen::MatrixXd X_M_Poportional_Ext(N*Q, P_p) ;
  for (int x = -1; x < N-1; ++x) {
    for (int j = (x+1)*Q ; j < (x+2)*Q; ++j){
      X_M_Poportional_Ext.row(j) = (X_init.rightCols(P_p)).row(x+1);
    }
  }
  X_EXT.conservativeResize( N*Q , X_M_Complete_Ext.cols()+X_M_Poportional_Ext.cols() );
  X_EXT << X_M_Complete_Ext, X_M_Poportional_Ext;

  // // Beta initialization
  Eigen::MatrixXd BETA;
  BETA = Eigen::MatrixXd::Zero(X_EXT.cols(),1);

  int iteration = 0;
  double check_tutz = 1.0;
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
  // for (int iteration=1; iteration < 18; iteration++){
  // while (check_tutz > 0.0001){
  double epsilon = 0.0001 ;
  while (Stop_criteria >( epsilon / N)){

    Eigen::MatrixXd Score_i = Eigen::MatrixXd::Zero(BETA.rows(),1);
    Eigen::MatrixXd F_i = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());

    double LogLik = 0.;

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
      }else if(distribution == "probit"){
        pi = ReferenceF::inverse_probit(eta);
        D = ReferenceF::inverse_derivative_probit(eta);
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
    Rcout << "Log Likelihood" << endl;
    Rcout << LogLik << endl;
    BETA = BETA + (F_i.inverse() * Score_i);
    // check_tutz = ((BETA - beta_old).norm())/(beta_old.norm()+check_tutz);
    iteration = iteration + 1;
    // Rcout << "Beta" << endl;
    // Rcout << BETA << endl;
  }
  Rcout << "Number of iterations" << endl;
  Rcout << iteration-1 << endl;
  return BETA;
}

// [[Rcpp::export]]
RCPP_MODULE(referencemodule){
  Rcpp::class_<ReferenceF>("ReferenceF")
  .constructor()
  .method( "GLMref", &ReferenceF::GLMref )
  ;
}
