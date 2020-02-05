#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>

class distribution{
public:
  double _epsilon_0 = 1e-10;
  double _epsilon_1 = 1e-6;
  distribution();
};

class Logistic : virtual distribution{
public:
  Eigen::VectorXd Quantile(Eigen::VectorXd vectordis1);
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  virtual Eigen::VectorXd in_open_corner(const Eigen::VectorXd& p) const;

  virtual double cdf_logit(const double& value) const;
  virtual double pdf_logit(const double& value) const;

  Logistic();
};

class Probit : public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);

  virtual double cdf_probit(const double& value) const;
  virtual double pdf_probit(const double& value) const;

  Probit();
};

class Cauchit : public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);

  virtual double cdf_cauchit(const double& value) const;

  Cauchit();
};

class Student : public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  Student();
};

class Gumbel : public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  Gumbel();
};

class Gompertz : public distribution{
public:
  Eigen::VectorXd InverseLinkCumulativeFunction(Eigen::VectorXd vectordis);
  Eigen::VectorXd InverseLinkDensityFunction(Eigen::VectorXd vectordis);
  Gompertz();
};

#endif
