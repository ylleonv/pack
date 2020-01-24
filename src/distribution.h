#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

<<<<<<< HEAD
=======

>>>>>>> a2c45f04bfe1c4bc8eb9be0778af64684ad35010
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
  arma::vec InverseLinkCumulativeFunction(arma::vec vectordis);
  arma::vec InverseLinkDensityFunction(arma::vec vectordis);
  virtual Eigen::VectorXd in_open_corner(const Eigen::VectorXd& p) const;

  virtual double cdf_logit(const double& value) const;
  virtual double pdf_logit(const double& value) const;

  Logistic();
};

class Probit : public distribution{
public:
  arma::vec InverseLinkCumulativeFunction(arma::vec vectordis);
  arma::vec InverseLinkDensityFunction(arma::vec vectordis);

  virtual double cdf_probit(const double& value) const;
  virtual double pdf_probit(const double& value) const;

  Probit();
};

class Cauchit : public distribution{
public:
  arma::vec InverseLinkCumulativeFunction(arma::vec vectordis);
  arma::vec InverseLinkDensityFunction(arma::vec vectordis);

  virtual double cdf_cauchit(const double& value) const;

  Cauchit();
};

class Student : public distribution{
public:
  arma::vec InverseLinkCumulativeFunction(arma::vec vectordis);
  arma::vec InverseLinkDensityFunction(arma::vec vectordis);
  Student();
};

class Gumbel : public distribution{
public:
  arma::vec InverseLinkCumulativeFunction(arma::vec vectordis);
  arma::vec InverseLinkDensityFunction(arma::vec vectordis);
  Gumbel();
};

class Gompertz : public distribution{
public:
  arma::vec InverseLinkCumulativeFunction(arma::vec vectordis);
  arma::vec InverseLinkDensityFunction(arma::vec vectordis);
  Gompertz();
};

#endif
