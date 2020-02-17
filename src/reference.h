#ifndef REFERENCEF_H_
#define REFERENCEF_H_
#include "distribution.h"

class ReferenceF : public virtual Logistic, Probit, Cauchit, Student, Gumbel, Gompertz{
public:
  ReferenceF();

  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta1) const;
  virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_probit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_probit(const Eigen::VectorXd& eta) const ;

  Eigen::MatrixXd GLMref(Eigen::MatrixXd X_M, Eigen::MatrixXd Y_V, std::string link,  std::string design);

};

#endif
