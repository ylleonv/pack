#include <iostream>
#include "distribution.h"
#include "referenceJean.h"
using namespace std;
using namespace Rcpp ;
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

ReferenceJean::ReferenceJean(void) {
  cout << "ReferenceJean is being created" << endl;
}

Eigen::VectorXd ReferenceJean::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi1( eta.size() );
  double norm1 = 1.;
  for(int j=0; j<eta.size(); ++j)
  {
    double number ;
    number = cdf_logit(eta(j));
    pi1[j] = number / ( 1-number );
    norm1 += pi1[j];
  }
  return (pi1/norm1);
}

Eigen::MatrixXd ReferenceJean::inverse_derivative_logistic(Eigen::VectorXd& eta2)
{
  Eigen::VectorXd pi1 = ReferenceJean::inverse_logistic(eta2);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
    // { D(j,j) = Logistic::pdf1( eta(j) ) /
    //   ( std::max(1e-10, std::min(1-1e-6,Logistic::cdf1(eta(j)))) *
    //     std::max(1e-10, std::min(1-1e-6, 1-Logistic::cdf1(eta(j)))) );
    // }
  { D1(j,j) = pdf_logit( eta2(j) ) /
    (cdf_logit(eta2(j)) * (1-cdf_logit(eta2(j))));
  }
  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}

Eigen::VectorXd ReferenceJean::inverse_probit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(size_t j=0; j<eta.size(); ++j)
  {
    pi[j] = cdf_probit( eta(j) ) / ( 1-cdf_probit( eta(j) ) );
    norm1 += pi[j];
  }
  return in_open_corner(pi/norm1);
}

Eigen::MatrixXd ReferenceJean::inverse_derivative_probit( Eigen::VectorXd& eta)
{
  Eigen::VectorXd pi = ReferenceJean::inverse_probit(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  for(size_t j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_probit( eta(j) ) /
    ( std::max(1e-10, std::min(1-1e-6,cdf_probit(eta(j)))) *
      std::max(1e-10, std::min(1-1e-6, 1-cdf_probit(eta(j)))) ); }
  return D * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose().eval() );
}

Eigen::VectorXd exponential(Eigen::VectorXd vector){
  for (size_t i = 0; i<=vector.size(); i++)
    vector[i] = exp(vector[i]);
  return vector;
}

Eigen::VectorXd lnatural(Eigen::VectorXd vector){
  for (size_t i = 0; i<=vector.size(); i++)
    vector[i] = log(vector[i]);
  return vector;
}

Eigen::MatrixXd ReferenceJean::GLMrefJean(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q, std::string link){
  // const int Q = K-1 ;
  // const int P = X_M.cols() -1 ;
  // const int N = X_M.n_cols ;
  Eigen::MatrixXd BETA = Eigen::MatrixXd::Zero((P+1)*Q,1);
  // Eigen::MatrixXd F_i;
  double llikelihood1 = 0;

  for (int iteration=1; iteration < 10; iteration++){
    int i = 0;
    Eigen::VectorXd llikelihood;
    double llikelihood2;

    Eigen::MatrixXd X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
    Eigen::VectorXd Y_M_i = Y_EXT.segment(i*Q , Q);
    Eigen::VectorXd eta = X_M_i * BETA;

    Eigen::VectorXd pi ;
    if(link == "logistic"){
      pi = inverse_logistic(eta);
    }else if(link == "probit"){
      pi = inverse_probit(eta);
    }

    // pi = inverse(eta);

    Eigen::MatrixXd Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());

    Eigen::MatrixXd D;

    if(link == "logistic"){
      D = inverse_derivative_logistic(eta);
    }else if(link == "probit"){
      D = inverse_derivative_probit(eta);
    }

    Eigen::MatrixXd W_in = D * Cov_i.inverse();
    Eigen::MatrixXd Score_i = X_M_i.transpose() * W_in * (Y_M_i - pi);
    Eigen::MatrixXd F_i = X_M_i.transpose() * (W_in*D) * (X_M_i);

    for (i=1; i < N; i++){
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_EXT.segment(i*Q ,Q);
      eta = X_M_i * BETA;

      if(link == "logistic"){
        pi = inverse_logistic(eta);
      }else if(link == "probit"){
        pi = inverse_probit(eta);
      }

      // pi = inverse(eta);

      Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());

      if(link == "logistic"){
        D = inverse_derivative_logistic(eta);
      }else if(link == "probit"){
        D = inverse_derivative_probit(eta);
      }

      W_in = D * Cov_i.inverse();
      Eigen::MatrixXd Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      Eigen::MatrixXd F_i_2 = X_M_i.transpose() * (W_in*D) * (X_M_i);
      F_i = F_i + F_i_2;

      // llikelihood = exponential(pi);
      // llikelihood2 = llikelihood.sum();
      //
      // llikelihood = llikelihood/llikelihood2;
      // llikelihood = lnatural(llikelihood);
      //
      // llikelihood2 = Y_M_i.transpose() * llikelihood;
      //
      // llikelihood1 = llikelihood1 + llikelihood2;

    }
    cout << "Beta" << endl;
    cout << BETA << endl;
    BETA = BETA + (F_i.inverse() * Score_i);
    //
    //     cout << "ll" << endl;
    //     cout << llikelihood1 << endl;
  }

  cout << "ll" << endl;
  cout << llikelihood1 << endl;
  return BETA;
}






