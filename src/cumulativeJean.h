#ifndef CUMULATIVEJEAN_H_
#define CUMULATIVEJEAN_H_
#include <RcppArmadillo.h>
#include "distribution.h"

class CumulativeJean : public Logistic, Probit, Cauchit, Student, Gumbel, Gompertz{
public:
  CumulativeJean();
  const arma::mat X_M;
  const arma::mat Y_M;

  virtual Eigen::VectorXd inverse(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative(const Eigen::VectorXd& eta) const;
  Eigen::MatrixXd GLMcumJean(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q);

};

#endif


