#include "binomial.h"
using namespace std;
using namespace Rcpp ;

FisherScoring::FisherScoring(void) {
  Rcpp::Rcout << "FisherScoring is being created" << std::endl;
}

Eigen::MatrixXd FisherScoring::GLMm(Eigen::MatrixXd X_M, Eigen::VectorXd Y_M, std::string link){
   //    Create initial beta
  const int N = X_M.rows() ;
  const int K = X_M.cols() ;
  // initialize betas to 0
  Eigen::MatrixXd beta = Eigen::MatrixXd::Zero(K,1);
  Eigen::VectorXd Mu;
  Eigen::VectorXd D_M;
  Eigen::VectorXd Deviance;
  double LogLik;

  double check_tutz = 1.0;
  double tol = 0.001;
  int n_iter = -1;

  //    algorithm
  while (check_tutz > tol){
    // Vector of probabilities:
    Eigen::MatrixXd eta = X_M * beta;
    if(link == "logistic"){
      Mu = Logistic::InverseLinkCumulativeFunction(eta);
      D_M = Logistic::InverseLinkDensityFunction(eta);
    }else if(link == "probit"){
      Mu = Normal::InverseLinkCumulativeFunction(eta);
      D_M = Normal::InverseLinkDensityFunction(eta);
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
    Eigen::MatrixXd D_M1 = Eigen::MatrixXd(D_M.asDiagonal());

    // Covariance
    Eigen::VectorXd Ones = Eigen::VectorXd::Ones(Mu.rows());
    Eigen::MatrixXd Covinv = ((Eigen::VectorXd(Mu.array()*(Ones-Mu).array())).asDiagonal()).inverse();

    Eigen::MatrixXd Score = X_M.transpose() * (D_M1 * Covinv) *  (Y_M-Mu);

    // W
    Eigen::MatrixXd W_M = (D_M1 * Covinv) * D_M1;
    // Second derivate - FisherInformation
    Eigen::MatrixXd dffm = (-X_M.transpose() * (W_M * X_M)).inverse() ;

    // Stop criteria Tutz
    Eigen::VectorXd beta_old = beta;
    Eigen::VectorXd beta_new = beta - (dffm * Score);
    check_tutz = ((beta_new - beta_old).norm())/(beta_old.norm());
    // // Deviance for ungrouped -> bernulli
    // // Deviance = -2*(Y_M.transpose()*log(Mu)) + ((1-Y_M.transpose())*log(1-Mu));

    // LogLik
    LogLik = (Y_M.transpose()*Eigen::VectorXd(Mu.array().log())) + ( ((Ones - Y_M).array() * ( (Ones - Mu).array()).log()).sum() );
    beta = beta_new;
    n_iter = n_iter + 1;
  }

  Rcpp::Rcout << "Number of iterations" << std::endl;
  Rcpp::Rcout << n_iter << std::endl;

  // Rcpp::Rcout << "Deviance" << std::endl;
  // Rcpp::Rcout << Deviance << std::endl;

  Rcpp::Rcout << "LogLik" << std::endl;
  Rcpp::Rcout << LogLik << std::endl;
  return beta;
}

RCPP_MODULE(fishder){
  Rcpp::class_<FisherScoring>("FisherScoring")
  .constructor()
  .method( "GLMm", &FisherScoring::GLMm )
  ;
}
