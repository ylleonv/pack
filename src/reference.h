#ifndef REFERENCEF_H_
#define REFERENCEF_H_
#include "distribution.h"

class ReferenceF : public virtual Logistic, Normal, Cauchit, Student, Gumbel, Gompertz{
public:

  ReferenceF();

  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta1) const;
  virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_normal(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_normal(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_cauchit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_cauchit(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_student(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_student(const Eigen::VectorXd& eta) const ;

  List GLMref(std::string response,
              StringVector explanatory_complete,
              StringVector explanatory_proportional,
              std::string distribution,
              SEXP categories_order,
              DataFrame dataframe);

  List GLMref_ec(std::string response, std::string actual_response,
                 std::string individuals,
                 StringVector explanatory_complete,
                 StringVector depend_y,
                 std::string distribution,
                 SEXP categories_order,
                 DataFrame dataframe,
                 std::string design);

};

#endif
