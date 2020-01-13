#ifndef CUMULATIVE_H_
#define CUMULATIVE_H_
#include <RcppArmadillo.h>
#include "distribution.h"

class CumulativeF : public Logistic, Probit, Cauchit, Student, Gumbel, Gompertz{
public:
  CumulativeF();
  const arma::mat X_M;
  const arma::mat Y_M;
  Eigen::VectorXd inverse(Eigen::VectorXd eta);
  Eigen::MatrixXd inverse_derivate(Eigen::VectorXd eta);
  Eigen::MatrixXd GLMcum(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q);
};

#endif


