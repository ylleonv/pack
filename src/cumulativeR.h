// #ifndef CUMULATIVER_H_
// #define CUMULATIVER_H_
// #include "distribution.h"
//
// class CumulativeR : virtual public distribution, Logistic, Normal, Cauchit{
// public:
//   CumulativeR();
//
//   virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta1) const;
//   virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;
//
//   virtual Eigen::VectorXd inverse_normal(const Eigen::VectorXd& eta) const;
//   virtual Eigen::MatrixXd inverse_derivative_normal(const Eigen::VectorXd& eta) const;
//
//   virtual Eigen::VectorXd inverse_cauchit(const Eigen::VectorXd& eta) const;
//   virtual Eigen::MatrixXd inverse_derivative_cauchit(const Eigen::VectorXd& eta) const;
//
//   Eigen::MatrixXd GLMcum(Eigen::MatrixXd X_M, Eigen::MatrixXd Y_V, std::string link,  std::string design);
//
//
// };
//
// #endif
//
//
