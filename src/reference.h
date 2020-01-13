#ifndef REFERENCE_H_
#define REFERENCE_H_
#include <RcppArmadillo.h>
#include "distribution.h"

class ReferenceF : public Logistic, Probit, Cauchit, Student, Gumbel, Gompertz{
public:
  ReferenceF();
  const arma::mat X_M;
  const arma::mat Y_M;
  Eigen::VectorXd evaluate(Eigen::VectorXd pi);
  Eigen::VectorXd inverse(Eigen::VectorXd eta);
  Eigen::MatrixXd inverse_derivate(Eigen::VectorXd eta);
  Eigen::MatrixXd GLMref(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q);
};

#endif


