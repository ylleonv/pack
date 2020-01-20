// #ifndef ADJACENTR_H_
// #define ADJACENTR_H_
// #include <RcppArmadillo.h>
// #include "distribution.h"
//
// class AdjacentR : public Logistic, Probit, Cauchit, Student, Gumbel, Gompertz{
// public:
//   AdjacentR();
//   const arma::mat X_M;
//   const arma::mat Y_M;
//
//   virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta) const;
//
//   virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;
//
//   virtual Eigen::VectorXd inverse_probit(const Eigen::VectorXd& eta) const;
//   virtual Eigen::MatrixXd inverse_derivative_probit(const Eigen::VectorXd& eta) const ;
//
//   Eigen::MatrixXd GLMadj(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q, std::string link);
// };
//
// #endif
//
//
