#ifndef ADJACENT_H_
#define ADJACENT_H_
#include <RcppArmadillo.h>
#include "distribution.h"

class AdjacentF : public Logistic, Probit, Cauchit, Student, Gumbel, Gompertz{
public:
  AdjacentF();
  const arma::mat X_M;
  const arma::mat Y_M;

  Eigen::VectorXd inverse(Eigen::VectorXd eta);
  Eigen::MatrixXd inverse_derivate(Eigen::VectorXd eta);
  Eigen::MatrixXd GLMadj(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q);
};

#endif


