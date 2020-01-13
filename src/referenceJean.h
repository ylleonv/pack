#ifndef REFERENCEJEAN_H_
#define REFERENCEJEAN_H_
#include <RcppArmadillo.h>
#include "distribution.h"

class ReferenceJean : public Logistic, Probit, Cauchit, Student, Gumbel, Gompertz{
public:
  ReferenceJean();
  const arma::mat X_M;
  const arma::mat Y_M;

  // virtual Eigen::VectorXd inverse( Eigen::VectorXd& eta, std::string link) ;

  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta1) const;
  virtual Eigen::MatrixXd inverse_derivative_logistic( Eigen::VectorXd& eta) ;

  virtual Eigen::VectorXd inverse_probit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_probit( Eigen::VectorXd& eta) ;

  Eigen::MatrixXd GLMrefJean(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q, std::string link);

};

#endif


