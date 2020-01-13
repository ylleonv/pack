// #include <iostream>
// #include "distribution.h"
// #include "cumulativeJean.h"
// using namespace std;
// using namespace Rcpp ;
// #include <eigen3/Eigen/Dense>
// #include <eigen3/Eigen/Core>
// #include <RcppArmadillo.h>
// #include <RcppEigen.h>
//
// // [[Rcpp::depends(RcppEigen)]]
// // [[Rcpp::depends(RcppArmadillo)]]
//
// CumulativeJean::CumulativeJean(void) {
//   cout << "cumulativeJean is being created" << endl;
// }
//
// Eigen::VectorXd CumulativeJean::inverse(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd ordered_pi( eta.size() );
//   ordered_pi[0] = Logistic::cdf1( eta(0) );
//   for(size_t j=1; j<eta.size(); ++j)
//   { ordered_pi[j] = Logistic::cdf1( eta(j) ) - Logistic::cdf1( eta(j-1) ); }
//   return in_open_corner(ordered_pi);
// }
//
// Eigen::MatrixXd CumulativeJean::inverse_derivative(const Eigen::VectorXd& eta) const
// {
//   Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
//   R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
//   Eigen::MatrixXd F = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
//   for(size_t j=0; j<eta.rows(); ++j)
//   { F(j,j) = Logistic::pdf1( eta(j) ); }
//   return (F * R);
// }
//
// Eigen::MatrixXd CumulativeJean::GLMcumJean(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q){
//   // const int Q = K-1 ;
//   // const int P = X_M.cols() -1 ;
//   // const int N = X_M.n_cols ;
//   Eigen::MatrixXd BETA = Eigen::MatrixXd::Zero((P+1)*Q,1);
//   // Eigen::MatrixXd F_i;
//   for (int iteration=1; iteration < 100; iteration++){
//     int i = 0;
//     Eigen::MatrixXd X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
//     Eigen::VectorXd Y_M_i = Y_EXT.segment(i*Q , Q);
//     Eigen::VectorXd eta = X_M_i * BETA;
//     Eigen::VectorXd pi = inverse(eta);
//     Eigen::MatrixXd Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
//     Eigen::MatrixXd D = inverse_derivative(eta);
//     Eigen::MatrixXd W_in = D * Cov_i.inverse();
//     Eigen::MatrixXd Score_i = X_M_i.transpose() * W_in * (Y_M_i - pi);
//     Eigen::MatrixXd F_i = X_M_i.transpose() * (W_in*D) * (X_M_i);
//     for (i=1; i < N; i++){
//       X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
//       Y_M_i = Y_EXT.segment(i*Q ,Q);
//       eta = X_M_i * BETA;
//       pi = inverse(eta);
//       Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
//       D = inverse_derivative(eta);
//       W_in = D * Cov_i.inverse();
//       Eigen::MatrixXd Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
//       Score_i = Score_i + Score_i_2;
//       Eigen::MatrixXd F_i_2 = X_M_i.transpose() * (W_in*D) * (X_M_i);
//       F_i = F_i + F_i_2;
//     }
//     BETA = BETA + (F_i.inverse() * Score_i);
//   }
//   return BETA;
// }
//
//
//
//
//
//
