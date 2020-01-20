// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "distribution.h"
#include <iostream>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/students_t.hpp>
using namespace boost::math;
using namespace std;
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace Rcpp ;

distribution::distribution(void) {
  cout << "Distribution is being created" << endl;
}

Eigen::VectorXd Logistic::in_open_corner(const Eigen::VectorXd& p) const
{
  Eigen::VectorXd pi = p;
  int J = pi.size() + 1;
  for(int j=0; j<J-1; ++j)
  { pi[j] = std::max(_epsilon_0, std::min(pi[j], 1-_epsilon_1)); }
  double sum = pi.sum();
  if(sum > 1-_epsilon_1)
  {
    for(int j=0; j<J-1; ++j)
    { pi[j] *= (1.-_epsilon_1)/sum;  }
  }
  return pi;
}

Logistic::Logistic(void) {
  cout << "Logistic is being created" << endl;
}

Probit::Probit(void) {
  cout << "Probit is being created" << endl;
}

// Methods
// Eigen::VectorXd Logistic::Quantile(Eigen::VectorXd vector){
//   logistic_distribution<> myLogistic1(0., 1.);
//   for (int i = 0; i<=vector.size()-1; i++)
//     vector[i] = quantile(myLogistic1, vector[i]);
//   return vector;
// }

arma::vec Logistic::InverseLinkCumulativeFunction(arma::vec vector){
  logistic_distribution<> myLogistic1(0., 1.);
  for (int i = 0; i<=vector.size()-1; i++)
    vector[i] = cdf(myLogistic1, vector[i]);
  return vector;
}

arma::vec Logistic::InverseLinkDensityFunction(arma::vec vector){
  logistic_distribution<> myLogistic1(0., 1.);
  for (int i = 0; i<=vector.size()-1; i++)
    vector[i] = pdf(myLogistic1, vector[i]);
  return vector;
}

// For constant values

double Logistic::cdf_logit(const double& value) const
{
  logistic_distribution<> myLogistic1(0., 1.);
  return cdf(myLogistic1, value);
}
double Logistic::pdf_logit(const double& value) const
{
  // logistic_distribution<> myLogistic1(0., 1.);
  // return pdf(myLogistic1, value);
  boost::math::logistic dist(0., 1.);
  return boost::math::pdf(dist, value);
}

double Probit::cdf_probit(const double& value) const
{
  boost::math::normal norm;
  return cdf(norm, value);
}
double Probit::pdf_probit(const double& value) const
{
  boost::math::normal norm;
  return cdf(norm, value);
}


arma::vec Probit::InverseLinkCumulativeFunction(arma::vec vector ){
  boost::math::normal norm;
  for (int i = 0; i<=vector.n_elem-1; i++)
    vector(i) = cdf(norm, vector(i));
  return vector;
}
arma::vec Probit::InverseLinkDensityFunction(arma::vec vector ){
  boost::math::normal norm;
  for (int i = 0; i<=vector.n_elem-1; i++)
    vector(i) = pdf(norm, vector(i));
  return vector;
}

Cauchit::Cauchit(void) {
  cout << "Cauchit is being created" << endl;
}

double Cauchit::cdf_cauchit(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  return cdf(extreme_value, value);
}

arma::vec Cauchit::InverseLinkCumulativeFunction(arma::vec vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.n_elem-1; i++)
    vector(i) = cdf(extreme_value, vector(i));
  return vector;
}
arma::vec Cauchit::InverseLinkDensityFunction(arma::vec vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.n_elem-1; i++)
    vector(i) = pdf(extreme_value, vector(i));
  return vector;
}

Student::Student(void) {
  cout << "Student is being created" << endl;
}
arma::vec Student::InverseLinkCumulativeFunction(arma::vec vector ){
  double _degrees = 1.0;
  boost::math::students_t_distribution<> student(_degrees);
  for (int i = 0; i<=vector.n_elem-1; i++)
    vector(i) = cdf(student, vector(i));
  return vector;
}
arma::vec Student::InverseLinkDensityFunction(arma::vec vector ){
  double _degrees = 1.0;
  boost::math::students_t_distribution<> student(_degrees);
  for (int i = 0; i<=vector.n_elem-1; i++)
    vector(i) = pdf(student, vector(i));
  return vector;
}

Gumbel::Gumbel(void) {
  cout << "Gumbel is being created" << endl;
}
arma::vec Gumbel::InverseLinkCumulativeFunction(arma::vec vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.n_elem-1; i++)
    vector(i) = cdf(extreme_value, vector(i));
  return vector;
}
arma::vec Gumbel::InverseLinkDensityFunction(arma::vec vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.n_elem-1; i++)
    vector(i) = pdf(extreme_value, vector(i));
  return vector;
}

Gompertz::Gompertz(void) {
  cout << "Gompertz is being created" << endl;
}
arma::vec Gompertz::InverseLinkCumulativeFunction(arma::vec vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.n_elem-1; i++)
    vector(i) = 1-cdf(extreme_value, -vector(i));
  return vector;
}
arma::vec Gompertz::InverseLinkDensityFunction(arma::vec vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.n_elem-1; i++)
    vector(i) = pdf(extreme_value, -vector(i));
  return vector;
}

