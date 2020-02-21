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

  List GLMadj(std::string response,
              StringVector explanatory_complete,
              StringVector explanatory_proportional,
              std::string distribution,
              NumericVector categories_order,
              DataFrame dataframe);
};

#endif


