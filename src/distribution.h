#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <RcppEigen.h>
using namespace std;
using namespace Rcpp;

class distribution{
public:
  double _epsilon_0 = 1e-10;
  double _epsilon_1 = 1e-6;
  DataFrame sort_by_user(DataFrame A, NumericVector order);
  Eigen::VectorXd sort_vector(Eigen::VectorXd x1) ;
  Eigen::MatrixXd sorted_rows(Eigen::MatrixXd A) ;
  Eigen::MatrixXd select_data(DataFrame x1, std::string response,
                              StringVector explanatory_complete,
                              StringVector explanatory_proportional,
                              NumericVector order) ;
  distribution(int x)  { cout << "distribution::distribution(int ) called" << endl;   }
  distribution();
};

class Logistic : virtual public distribution{
public:
  Eigen::VectorXd Quantile(Eigen::VectorXd vectordis1);
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  virtual Eigen::VectorXd in_open_corner(const Eigen::VectorXd& p) const;

  virtual double cdf_logit(const double& value) const;
  virtual double pdf_logit(const double& value) const;

  Logistic(int x):distribution(x)   {
    cout<<"Logistic::Logistic(int ) called"<< endl;
  }

  Logistic();
};

class Probit : virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);

  virtual double cdf_probit(const double& value) const;
  virtual double pdf_probit(const double& value) const;

  Probit(int x):distribution(x)   {
    cout<<"Probit::Probit(int ) called"<< endl;
  }


  Probit();
};

class Cauchit : virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);

  virtual double cdf_cauchit(const double& value) const;
  virtual double pdf_cauchit(const double& value) const;

  Cauchit(int x):distribution(x)   {
    cout<<"Cauchit::Cauchit(int ) called"<< endl;
  }

  Cauchit();
};

class Student :  virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  Student();
};

class Gumbel :  virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  Gumbel();
};

class Gompertz : virtual public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);

  virtual double cdf_gompertz(const double& value) const;
  virtual double pdf_gompertz(const double& value) const;

  Gompertz();
};


#endif
