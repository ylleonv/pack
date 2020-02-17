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

Eigen::VectorXd sort_vector(Eigen::VectorXd x1) {
  NumericVector V(x1.data(), x1.data() + x1.size());
  int x=0;
  std::iota(V.begin(),V.end(),x++);
  sort( V.begin(),V.end(), [&](int i,int j){return x1[i]<x1[j];} );
  Eigen::Map<Eigen::VectorXd> XS(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(V));
  return XS;
}

Eigen::MatrixXd sorted_rows(Eigen::MatrixXd A)
{
  Eigen::VectorXd vec1 = sort_vector(A.col(0));
  Eigen::MatrixXd B = A.row(vec1(0));
  for (int i = 1; i < A.rows(); ++i) {
    B.conservativeResize(B.rows()+1, B.cols());
    B.row(B.rows()-1) = A.row(vec1(i));
  }
  return B;
}

Eigen::MatrixXd ReferenceF::GLMref(Eigen::MatrixXd X_M, Eigen::MatrixXd Y_V, std::string link, std::string design){

  // Number of explanatory variables
  const int P = X_M.cols() ;

  // Number of observations
  const int N = X_M.rows() ;

  // Create full matrix: Y,Xs
  Eigen::MatrixXd Full_M(X_M.rows(), 1 + X_M.cols());
  Full_M << Y_V, X_M;

  // Order by Y
  Eigen::MatrixXd Full_M_ordered = sorted_rows(Full_M);

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
  Eigen::MatrixXd BETA;
  Eigen::MatrixXd X_EXT(2, 2);

  if(design == "complete"){
    // Intercept
    Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(X_init.rows());
    X_EXT.conservativeResize(Ones1.rows(), 1 + X_init.cols());
    X_EXT << Ones1, X_init;
    // Expand matrix
    X_EXT = kroneckerProduct(X_EXT,Eigen::MatrixXd::Identity(Q,Q)).eval();
    // Beta initialization
    BETA = Eigen::MatrixXd::Zero((P+1)*Q,1);
  }else if(design == "proportional"){
    Eigen::MatrixXd Ones2 = Eigen::MatrixXd::Ones(X_init.rows(),1);
    // Eigen::MatrixXd X_EXT1 = kroneckerProduct(Ones2, X_init);
    Eigen::MatrixXd Intercept_proportinal = kroneckerProduct(Ones2, Eigen::MatrixXd::Identity(Q,Q)).eval();
    Eigen::MatrixXd X_EXT2((X_init.rows())*Q, P) ;
    for (int x = -1; x < ((X_init.rows()))-1; ++x) {
      for (int j = (x+1)*Q ; j < (x+2)*Q; ++j){
        X_EXT2.row(j) = X_init.row(x+1);
      }
    }
    X_EXT.conservativeResize(((X_init.rows())*Q), Q+P);
    X_EXT << Intercept_proportinal, X_EXT2;
    // Beta initialization
    BETA = Eigen::MatrixXd::Zero((P+Q),1);
  }

  int iteration = 1;
  double check_tutz = 1.0;
  Eigen::MatrixXd X_M_i ;
  Eigen::VectorXd Y_M_i ;
  Eigen::VectorXd eta ;
  Eigen::VectorXd pi ;
  Eigen::MatrixXd D ;
  Eigen::MatrixXd Cov_i ;
  Eigen::MatrixXd W_in ;
  Eigen::MatrixXd Score_i_2 ;
  Eigen::MatrixXd F_i_2 ;

  // for (int iteration=1; iteration < 40; iteration++){
  while (check_tutz > 1e-6){

    // // Elemnts for cumulative sums
    // Eigen::MatrixXd Score_i = Eigen::MatrixXd::Zero((P+1)*Q,1);
    // Eigen::MatrixXd F_i = Eigen::MatrixXd::Zero((P+1)*Q,(P+1)*Q);

    Eigen::MatrixXd Score_i = Eigen::MatrixXd::Zero(BETA.rows(),1);
    Eigen::MatrixXd F_i = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());

    double LogLik = 0.;

    // Loop by subject

    for (int i=0; i < N; i++){
      // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());

      // X_M_i = X_EXT3.block(i*Q , 0 , Q , Q+P);;

      Y_M_i = Y_init.row(i);
      eta = X_M_i * BETA;

      // Vector pi depends on selected link

      if(link == "logistic"){
        pi = ReferenceF::inverse_logistic(eta);
        D = ReferenceF::inverse_derivative_logistic(eta);
      }else if(link == "probit"){
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
    // Stop criteria Tutz
    Eigen::VectorXd beta_old = BETA;
    BETA = BETA + (F_i.inverse() * Score_i);
    check_tutz = ((BETA - beta_old).norm())/(beta_old.norm());
    Rcout << "Log Likelihood" << endl;
    Rcout << LogLik << endl;
    iteration = iteration + 1;
  }
  Rcout << "Number of iterations" << endl;
  Rcout << iteration << endl;
  return BETA;
}

RCPP_MODULE(referencemodule){
  Rcpp::class_<ReferenceF>("ReferenceF")
  .constructor()
  .method( "GLMref", &ReferenceF::GLMref )
  ;
}
