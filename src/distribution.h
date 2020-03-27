#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <RcppEigen.h>
using namespace std;
using namespace Rcpp;

class distribution{
public:
  double _epsilon_0 = 1e-10;
  double _epsilon_1 = 1e-6;

  std::string concatenate(std::string x, std::string level);

  List select_data(DataFrame x1, std::string response,
                              StringVector explanatory_complete,
                              StringVector explanatory_proportional,
                              SEXP order) ;

  List select_data_nested(DataFrame x1, std::string response,std::string actual_response,
                                        std::string individuals,
                                        StringVector explanatory_complete,
                                        StringVector depend_y,
                                        SEXP order);

  // distribution(int x);
  distribution();
};

class Logistic : virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  virtual Eigen::VectorXd in_open_corner(const Eigen::VectorXd& p) const;

  virtual double cdf_logit(const double& value) const;
  virtual double pdf_logit(const double& value) const;

  // Logistic(int x):distribution()   {
  //   cout<<"Logistic::Logistic(int ) called"<< endl;
  // }

  Logistic();
};

class Normal : virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  virtual double cdf_normal(const double& value) const;
  virtual double pdf_normal(const double& value) const;

  // Normal(int x):distribution()   {
  //   cout<<"normal::normal(int ) called"<< endl;
  // }

  Normal();
};

class Cauchit : virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  virtual double cdf_cauchit(const double& value) const;
  virtual double pdf_cauchit(const double& value) const;

  // Cauchit(int x):distribution()   {
  //   cout<<"Cauchit::Cauchit(int ) called"<< endl;
  // }

  Cauchit();
};

class Student :  virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);

  virtual double cdf_student(const double& value, double df_student) const;
  virtual double pdf_student(const double& value, double df_student) const;

  Student();
};

class Gumbel :  virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);

  virtual double cdf_gumbel(const double& value) const;
  virtual double pdf_gumbel(const double& value) const;

  Gumbel();
};

class Gompertz : virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  virtual double cdf_gompertz(const double& value) const;
  virtual double pdf_gompertz(const double& value) const;

  Gompertz();
};


#endif
