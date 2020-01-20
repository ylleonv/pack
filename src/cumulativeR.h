#ifndef CUMULATIVER_H_
#define CUMULATIVER_H_
#include <RcppArmadillo.h>
#include "distribution.h"
#include <boost/math/distributions/logistic.hpp>

class CumulativeR : public Logistic, Probit, Cauchit{
public:
  CumulativeR();
  const arma::mat X_M;
  const arma::mat Y_M;

  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta1) const;
  virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_probit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_probit(const Eigen::VectorXd& eta) const;

  virtual Eigen::VectorXd inverse_cauchit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_cauchit(const Eigen::VectorXd& eta) const;

  Eigen::MatrixXd GLMcum(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q, std::string link);


};

#endif


