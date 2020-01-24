#include <iostream>
#include "distribution.h"
#include "sequentialR.h"
using namespace std;
using namespace Rcpp ;
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

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

Eigen::VectorXd SequentialR::inverse_probit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  double product = 1;
  for(size_t j=0; j<eta.size(); ++j)
  {
    ordered_pi[j] = product * Probit::cdf_probit( eta(j) );
    product *= ( 1 - Probit::cdf_probit( eta(j) ) );
  }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd SequentialR::inverse_derivative_probit(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  double product = 1.;
  for (size_t j=0; j < eta.rows(); ++j)
  {
    M(j,j) = Probit::pdf_probit(eta(j)) * product;
    for (size_t i=0; i<j; ++i)
    { M(i,j) = - Probit::pdf_probit(eta(i))  * std::max(1e-10, std::min(Probit::cdf_probit(eta(j)), 1-1e-6)) * product / std::max(1e-10, std::min( 1-Probit::cdf_probit(eta(i)), 1-1e-6)); }
    product *= std::max(1e-10, std::min( 1-Probit::cdf_probit(eta(j)), 1-1e-6));
  }
  return M;
}

// Eigen::VectorXd exponential(Eigen::VectorXd vector){
//   for (size_t i = 0; i<=vector.size(); i++)
//     vector[i] = exp(vector[i]);
//   return vector;
// }
//
// Eigen::VectorXd lnatural(Eigen::VectorXd vector){
//   for (size_t i = 0; i<=vector.size(); i++)
//     vector[i] = log(vector[i]);
//   return vector;
// }

Eigen::MatrixXd SequentialR::GLMseq(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q, std::string link){
  // const int Q = K-1 ;
  // const int P = X_M.cols() -1 ;
  // const int N = X_M.n_cols ;
  Eigen::MatrixXd BETA = Eigen::MatrixXd::Zero((P+1)*Q,1);
  // Eigen::MatrixXd F_i;
  // double llikelihood1 = 0;

  for (int iteration=1; iteration < 40; iteration++){
    int i = 0;
    // Eigen::VectorXd llikelihood;
    // double llikelihood2;

    Eigen::MatrixXd X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
    Eigen::VectorXd Y_M_i = Y_EXT.segment(i*Q , Q);
    Eigen::VectorXd eta = X_M_i * BETA;

    Eigen::VectorXd pi ;
    if(link == "logistic"){
      pi = SequentialR::inverse_logistic(eta);
    }else if(link == "probit"){
      pi = SequentialR::inverse_probit(eta);
    }

    Eigen::MatrixXd Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());

    Eigen::MatrixXd D ;

    if(link == "logistic"){
      D = SequentialR::inverse_derivative_logistic(eta);
    }else if(link == "probit"){
      D = SequentialR::inverse_derivative_probit(eta);
    }

    Eigen::MatrixXd W_in = D * Cov_i.inverse();
    Eigen::MatrixXd Score_i = X_M_i.transpose() * W_in * (Y_M_i - pi);
    Eigen::MatrixXd F_i = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);

    for (i=1; i < N; i++){
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_EXT.segment(i*Q ,Q);
      eta = X_M_i * BETA;

      // pi = SequentialR::inverse_logistic(eta);

      if(link == "logistic"){
        pi = SequentialR::inverse_logistic(eta);
      }else if(link == "probit"){
        pi = SequentialR::inverse_probit(eta);
      }


      Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());

      // D = SequentialR::inverse_derivative_logistic(eta);

      if(link == "logistic"){
        D = SequentialR::inverse_derivative_logistic(eta);
      }else if(link == "probit"){
        D = SequentialR::inverse_derivative_probit(eta);
      }

      W_in = D * Cov_i.inverse();
      Eigen::MatrixXd Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      Eigen::MatrixXd F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
      F_i = F_i + F_i_2;
    }

    BETA = BETA + (F_i.inverse() * Score_i);

    //     Rcout << "ll" << endl;
    //     Rcout << llikelihood1 << endl;

    // Rcout << "Beta" << endl;
    // Rcout << BETA << endl;
  }
  // Rcout << "ll" << endl;
  // Rcout << llikelihood1 << endl;
  return BETA;
}





