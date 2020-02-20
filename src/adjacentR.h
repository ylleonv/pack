#ifndef ADJACENTR_H_
#define ADJACENTR_H_
#include "distribution.h"

class AdjacentR : virtual public Logistic, Probit, Cauchit, Student, Gumbel, Gompertz{
public:
  AdjacentR();


  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_probit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_probit(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_cauchit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_cauchit(const Eigen::VectorXd& eta) const ;

  Eigen::MatrixXd GLMadj(Eigen::MatrixXd X_M, Eigen::MatrixXd Y_V, std::string link,  std::string design);
};

#endif


