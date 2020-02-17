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

  Eigen::MatrixXd GLMseq(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q, std::string link);


};

#endif

