#ifndef ADJACENTJEAN_H_
#define ADJACENTJEAN_H_
#include <RcppArmadillo.h>
#include "distribution.h"

class AdjacentJean : public Logistic{
public:
  AdjacentJean();
  const arma::mat X_M;
  const arma::mat Y_M;

  virtual Eigen::VectorXd inverse(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative(const Eigen::VectorXd& eta) const;
  Eigen::MatrixXd GLMadjJean(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q);

};

#endif


