#include <iostream>
#include "distribution.h"
#include "reference.h"
using namespace std;
using namespace Rcpp ;
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

ReferenceF::ReferenceF(void) {
  cout << "Reference is being created" << endl;
}

Eigen::VectorXd ReferenceF::evaluate(Eigen::VectorXd pi)
{
  Eigen::VectorXd eta(pi.size());
  for(size_t j=0; j<pi.size(); ++j)
  {
    eta[j] = (pi[j]/(pi[j]+1-pi.sum()));
  }
  eta = Logistic::Quantile(eta);
  return eta;
}

Eigen::VectorXd ReferenceF::inverse(Eigen::VectorXd eta)
{
  Eigen::VectorXd pi(eta.size());
  pi = Logistic::InverseLinkCumulativeFunction(eta);
  double norm = 1.;
  for(int j=0; j<eta.size(); ++j)
  {
    pi[j] = pi[j] / ( 1- pi[j] );
    norm += pi[j];
  }
  return pi/norm;
}

  Eigen::MatrixXd ReferenceF::inverse_derivate(Eigen::VectorXd eta)
  {
    Eigen::VectorXd pi = inverse(eta);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
    Eigen::VectorXd di = Logistic::InverseLinkDensityFunction(eta);
    Eigen::VectorXd pi_p = Logistic::InverseLinkCumulativeFunction(eta);
    for(int j=0; j<eta.size(); ++j)
    {
      D(j,j) = di[j] / (pi_p[j] *( 1 - pi_p[j]));
    }
    return D * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose().eval() );
  }


Eigen::MatrixXd ReferenceF::GLMref(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q){
  // const int Q = K-1 ;
  // const int P = X_M.cols() -1 ;
  // const int N = X_M.n_cols ;
  Eigen::MatrixXd BETA = Eigen::MatrixXd::Zero((P+1)*Q,1);
  // Eigen::MatrixXd F_i;
  for (int iteration=1; iteration < 10; iteration++){
    int i = 0;
    Eigen::MatrixXd X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
    Eigen::VectorXd Y_M_i = Y_EXT.segment(i*Q , Q);
    Eigen::VectorXd eta = X_M_i * BETA;
    Eigen::VectorXd pi = inverse(eta);
    Eigen::MatrixXd Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
    Eigen::MatrixXd D = inverse_derivate(eta);
    Eigen::MatrixXd W_in = D * Cov_i.inverse();
    Eigen::MatrixXd Score_i = X_M_i.transpose() * W_in * (Y_M_i - pi);
    Eigen::MatrixXd F_i = X_M_i.transpose() * (W_in*D) * (X_M_i);
    for (i=1; i < N; i++){
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_EXT.segment(i*Q ,Q);
      eta = X_M_i * BETA;
      pi = inverse(eta);
      Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
      D = inverse_derivate(eta);
      W_in = D * Cov_i.inverse();
      Eigen::MatrixXd Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      Eigen::MatrixXd F_i_2 = X_M_i.transpose() * (W_in*D) * (X_M_i);
      F_i = F_i + F_i_2;
    }
    BETA = BETA + (F_i.inverse() * Score_i);
  }
  return BETA;
}








