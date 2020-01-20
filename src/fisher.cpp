#include <iostream>
#include "distribution.h"
#include "fisher.h"
using namespace std;
using namespace Rcpp ;
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

FisherScoring::FisherScoring(void) {
  cout << "FisherScoring is being created" << endl;
}

arma::mat FisherScoring::GLMm(arma::mat X_M, arma::mat Y_M, std::string link){

  //    Create initial beta
  const int N = X_M.n_rows ;
  const int K = X_M.n_cols ;
  arma::mat beta(K, 1);
  beta.fill(0.0) ; // initialize betas to 0

  // logistic = false;
  arma::vec Mu;
  arma::mat D_M;
  arma::mat Cov;
  Eigen::MatrixXd eigen_B;
  double check_tutz = 1.0;
  double tol = 0.001;
  int n_iter = -1;

  arma::vec Deviance;

  //    algorithm
  while (check_tutz > tol){
    // Vector of probabilities:
    arma::mat eta = X_M * beta;
    if(link == "logistic"){
      Mu = Logistic::InverseLinkCumulativeFunction(eta);
      D_M = Logistic::InverseLinkDensityFunction(eta);
    }else if(link == "probit"){
      Mu = Probit::InverseLinkCumulativeFunction(eta);
      D_M = Probit::InverseLinkDensityFunction(eta);
    }else if(link == "cauchit"){
      Mu = Cauchit::InverseLinkCumulativeFunction(eta);
      D_M = Cauchit::InverseLinkDensityFunction(eta);
    }else if(link == "student"){
      Mu = Student::InverseLinkCumulativeFunction(eta);
      D_M = Student::InverseLinkDensityFunction(eta);
    }else if(link == "gumbel"){
      Mu = Gumbel::InverseLinkCumulativeFunction(eta);
      D_M = Gumbel::InverseLinkDensityFunction(eta);
    }else if(link == "gompertz"){
      Mu = Gompertz::InverseLinkCumulativeFunction(eta);
      D_M = Gompertz::InverseLinkDensityFunction(eta);
    }


    //  D
    D_M = arma::diagmat(D_M);

    //  Covariance
    Cov = diagvec(Mu*(1-Mu).t());
    Cov = diagmat(Cov);

    Eigen::MatrixXd CovEig = Eigen::Map<Eigen::MatrixXd>(Cov.memptr(),
                                      Cov.n_rows,
                                      Cov.n_cols);
    CovEig = CovEig.inverse();
    arma::mat Covinv = arma::mat(CovEig.data(), CovEig.rows(), CovEig.cols(),
                     false, false);

    // First derivative - ScoreFunction:
    arma::mat Score = X_M.t() * (D_M * Covinv) *  (Y_M-Mu);

    //  W
    arma::mat W_M = (D_M * Covinv) * D_M;

    // Why this line does not work?
    // arma::mat W_M = D_M * inv(diagmat(Cov)) * D_M;

    // //   Second derivate - FisherInformation
    arma::mat dffm;
    dffm = -X_M.t() * (W_M * X_M) ;

    eigen_B = Eigen::Map<Eigen::MatrixXd>(dffm.memptr(),
                                          dffm.n_rows,
                                          dffm.n_cols);
    eigen_B = eigen_B.inverse();
    dffm = arma::mat(eigen_B.data(), eigen_B.rows(), eigen_B.cols(),
                     false, false);

    // Stop criteria Tutz
    arma::vec beta_old = beta;
    arma::vec beta_new = beta - (dffm * Score);
    check_tutz = (arma::norm(beta_new - beta_old))/arma::norm(beta_old);
    // cout << n_iter << endl;
    // cout << check_tutz << endl;

    // Deviance for ungrouped -> bernulli
    Deviance = -2*(Y_M.t()*log(Mu)) + ((1-Y_M.t())*log(1-Mu));

    beta = beta_new;
    // beta.print();

    n_iter = n_iter + 1;
  }

  Rcout << "Number of iterations" << endl;
  Rcout << n_iter << endl;

  Rcout << "Deviance" << endl;
  Rcout << Deviance << endl;

  return beta;



}
