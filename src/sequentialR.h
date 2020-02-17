#ifndef SEQUENTIALR_H_
#define SEQUENTIALR_H_
#include "distribution.h"

class SequentialR : public Logistic, Probit{
public:
  SequentialR();

  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta1) const;
  virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_probit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_probit(const Eigen::VectorXd& eta) const;

  Eigen::MatrixXd GLMseq(Eigen::MatrixXd X_M, Eigen::MatrixXd Y_V, std::string link,  std::string design);


};

#endif

