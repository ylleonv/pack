// #include <iostream>
// #include "distribution.h"
// #include "adjacentJean.h"
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
// AdjacentJean::AdjacentJean(void) {
//   cout << "AdjacentJean is being created" << endl;
// }
//
// Eigen::VectorXd AdjacentJean::inverse(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi( eta.size() );
//   pi[eta.size()-1] = Logistic::cdf1( eta(eta.size()-1) ) / ( 1-Logistic::cdf1( eta(eta.size()-1) ) );
//   double norm = 1 + pi[eta.size()-1];
//   for(size_t j=(eta.size()-1); j>0; --j)
//   {
//     pi[j-1] = pi[j] * Logistic::cdf1( eta(j-1) ) / ( 1-Logistic::cdf1( eta(j-1) ) );
//     norm += pi[j-1];
//   }
//   return in_open_corner(pi/norm);
// }
//
// Eigen::MatrixXd AdjacentJean::inverse_derivative(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi = inverse(eta);
//   Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
//   Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
//   for(size_t j=0; j<pi.rows(); ++j)
//   { D(j,j) = Logistic::pdf1( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, Logistic::cdf1(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-Logistic::cdf1(eta(j)))) ); }
//
//   return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );
//
// }
//
// Eigen::MatrixXd AdjacentJean::GLMadjJean(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q){
//   // const int Q = K-1 ;
//   // const int P = X_M.cols() -1 ;
//   // const int N = X_M.n_cols ;
//   Eigen::MatrixXd BETA = Eigen::MatrixXd::Zero((P+1)*Q,1);
//   // Eigen::MatrixXd F_i;
//   for (int iteration=1; iteration < 10; iteration++){
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
