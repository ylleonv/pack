// #include <iostream>
// #include "distribution.h"
// #include "adjacentR.h"
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
// AdjacentR::AdjacentR(void) {
//   cout << "Reference is being created" << endl;
// }
//
// Eigen::VectorXd AdjacentR::inverse_logistic(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi( eta.size() );
//   pi[eta.size()-1] = cdf_logit( eta(eta.size()-1) ) / ( 1-cdf_logit( eta(eta.size()-1) ) );
//   double norm = 1 + pi[eta.size()-1];
//   for(size_t j=(eta.size()-1); j>0; --j)
//   {
//     pi[j-1] = pi[j] * cdf_logit( eta(j-1) ) / ( 1-cdf_logit( eta(j-1) ) );
//     norm += pi[j-1];
//   }
//   return in_open_corner(pi/norm);
// }
//
// Eigen::MatrixXd AdjacentR::inverse_derivative_logistic(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi = AdjacentR::inverse_logistic(eta);
//   Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
//   Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
//   for(size_t j=0; j<pi.rows(); ++j)
//   { D(j,j) = pdf_logit( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_logit(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_logit(eta(j)))) ); }
//
//   return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );
//
// }
//
// Eigen::VectorXd AdjacentR::inverse_probit(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi( eta.size() );
//   pi[eta.size()-1] = cdf_probit( eta(eta.size()-1) ) / ( 1-cdf_probit( eta(eta.size()-1) ) );
//   double norm = 1 + pi[eta.size()-1];
//   for(size_t j=(eta.size()-1); j>0; --j)
//   {
//     pi[j-1] = pi[j] * cdf_probit( eta(j-1) ) / ( 1-cdf_probit( eta(j-1) ) );
//     norm += pi[j-1];
//   }
//   return in_open_corner(pi/norm);
// }
//
// Eigen::MatrixXd AdjacentR::inverse_derivative_probit(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi = AdjacentR::inverse_probit(eta);
//   Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
//   Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
//   for(size_t j=0; j<pi.rows(); ++j)
//   { D(j,j) = pdf_probit( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_probit(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_probit(eta(j)))) ); }
//
//   return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );
//
// }
// // Eigen::VectorXd exponential(Eigen::VectorXd vector){
// //   for (size_t i = 0; i<=vector.size(); i++)
// //     vector[i] = exp(vector[i]);
// //   return vector;
// // }
// //
// // Eigen::VectorXd lnatural(Eigen::VectorXd vector){
// //   for (size_t i = 0; i<=vector.size(); i++)
// //     vector[i] = log(vector[i]);
// //   return vector;
// // }
//
// Eigen::MatrixXd AdjacentR::GLMadj(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q, std::string link){
//   // const int Q = K-1 ;
//   // const int P = X_M.cols() -1 ;
//   // const int N = X_M.n_cols ;
//   Eigen::MatrixXd BETA = Eigen::MatrixXd::Zero((P+1)*Q,1);
//   // Eigen::MatrixXd F_i;
//   double llikelihood1 = 0;
//
//   for (int iteration=1; iteration < 10; iteration++){
//     int i = 0;
//     Eigen::VectorXd llikelihood;
//     double llikelihood2;
//
//     Eigen::MatrixXd X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
//     Eigen::VectorXd Y_M_i = Y_EXT.segment(i*Q , Q);
//     Eigen::VectorXd eta = X_M_i * BETA;
//
//     Eigen::VectorXd pi ;
//     if(link == "logistic"){
//       pi = AdjacentR::inverse_logistic(eta);
//     }else if(link == "probit"){
//       pi = AdjacentR::inverse_probit(eta);
//     }
//
//     // pi = inverse(eta);
//
//     Eigen::MatrixXd Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
//
//     Eigen::MatrixXd D;
//
//     if(link == "logistic"){
//       D = AdjacentR::inverse_derivative_logistic(eta);
//     }else if(link == "probit"){
//       D = AdjacentR::inverse_derivative_probit(eta);
//     }
//
//     Eigen::MatrixXd W_in = D * Cov_i.inverse();
//     Eigen::MatrixXd Score_i = X_M_i.transpose() * W_in * (Y_M_i - pi);
//     Eigen::MatrixXd F_i = X_M_i.transpose() * (W_in*D) * (X_M_i);
//
//     for (i=1; i < N; i++){
//       X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
//       Y_M_i = Y_EXT.segment(i*Q ,Q);
//       eta = X_M_i * BETA;
//
//       if(link == "logistic"){
//         pi = AdjacentR::inverse_logistic(eta);
//       }else if(link == "probit"){
//         pi = AdjacentR::inverse_probit(eta);
//       }
//
//       // pi = inverse(eta);
//
//       Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
//
//       if(link == "logistic"){
//         D = AdjacentR::inverse_derivative_logistic(eta);
//       }else if(link == "probit"){
//         D = AdjacentR::inverse_derivative_probit(eta);
//       }
//
//       W_in = D * Cov_i.inverse();
//       Eigen::MatrixXd Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
//       Score_i = Score_i + Score_i_2;
//       Eigen::MatrixXd F_i_2 = X_M_i.transpose() * (W_in*D) * (X_M_i);
//       F_i = F_i + F_i_2;
//
//       // llikelihood = exponential(pi);
//       // llikelihood2 = llikelihood.sum();
//       //
//       // llikelihood = llikelihood/llikelihood2;
//       // llikelihood = lnatural(llikelihood);
//       //
//       // llikelihood2 = Y_M_i.transpose() * llikelihood;
//       //
//       // llikelihood1 = llikelihood1 + llikelihood2;
//
//     }
//     cout << "Beta" << endl;
//     cout << BETA << endl;
//     BETA = BETA + (F_i.inverse() * Score_i);
//     //
//     //     cout << "ll" << endl;
//     //     cout << llikelihood1 << endl;
//   }
//
//   cout << "ll" << endl;
//   cout << llikelihood1 << endl;
//   return BETA;
// }
//
//
//
//
//
//
