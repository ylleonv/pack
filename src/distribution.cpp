// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

#include "distribution.h"
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/students_t.hpp>

using namespace boost::math;
using namespace std;
using namespace Rcpp ;

distribution::distribution(void) {
  Rcout << "Distribution is being created" << endl;
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
  Rcout << "Logistic is being created" << endl;
}

Probit::Probit(void) {
  Rcout << "Probit is being created" << endl;
}

Eigen::VectorXd Logistic::InverseLinkCumulativeFunction(Eigen::VectorXd vector){
  boost::math::logistic dist(0., 1.);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = boost::math::cdf(dist, vector(i));
  return vector;
}

Eigen::VectorXd Logistic::InverseLinkDensityFunction(Eigen::VectorXd vector){
  boost::math::logistic dist(0., 1.);
  for (int i = 0; i<=vector.size()-1; i++)
    vector(i) = boost::math::pdf(dist, vector(i));
  return vector;
}

// For constant values
double Logistic::cdf_logit(const double& value) const
{
  boost::math::logistic dist(0., 1.);
  return boost::math::cdf(dist, value);
}
double Logistic::pdf_logit(const double& value) const
{
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


Eigen::VectorXd Probit::InverseLinkCumulativeFunction(Eigen::VectorXd vector ){
  boost::math::normal norm;
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = cdf(norm, vector(i));
  return vector;
}
Eigen::VectorXd Probit::InverseLinkDensityFunction(Eigen::VectorXd vector ){
  boost::math::normal norm;
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = pdf(norm, vector(i));
  return vector;
}

Cauchit::Cauchit(void) {
  Rcout << "Cauchit is being created" << endl;
}

double Cauchit::cdf_cauchit(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  return cdf(extreme_value, value);
}
Eigen::VectorXd Cauchit::InverseLinkCumulativeFunction(Eigen::VectorXd vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = cdf(extreme_value, vector(i));
  return vector;
}
Eigen::VectorXd Cauchit::InverseLinkDensityFunction(Eigen::VectorXd vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = pdf(extreme_value, vector(i));
  return vector;
}

Student::Student(void) {
  Rcout << "Student is being created" << endl;
}
Eigen::VectorXd Student::InverseLinkCumulativeFunction(Eigen::VectorXd vector ){
  double _degrees = 1.0;
  boost::math::students_t_distribution<> student(_degrees);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = cdf(student, vector(i));
  return vector;
}
Eigen::VectorXd Student::InverseLinkDensityFunction(Eigen::VectorXd vector ){
  double _degrees = 1.0;
  boost::math::students_t_distribution<> student(_degrees);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = pdf(student, vector(i));
  return vector;
}

Gumbel::Gumbel(void) {
  Rcout << "Gumbel is being created" << endl;
}
Eigen::VectorXd Gumbel::InverseLinkCumulativeFunction(Eigen::VectorXd vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = cdf(extreme_value, vector(i));
  return vector;
}
Eigen::VectorXd Gumbel::InverseLinkDensityFunction(Eigen::VectorXd vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = pdf(extreme_value, vector(i));
  return vector;
}

Gompertz::Gompertz(void) {
  Rcout << "Gompertz is being created" << endl;
}
Eigen::VectorXd Gompertz::InverseLinkCumulativeFunction(Eigen::VectorXd vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = 1-cdf(extreme_value, -vector(i));
  return vector;
}
Eigen::VectorXd Gompertz::InverseLinkDensityFunction(Eigen::VectorXd vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = pdf(extreme_value, -vector(i));
  return vector;
}

RCPP_MODULE(exportmod){
  using namespace Rcpp ;
  class_<distribution>("distribution")
    .constructor()
  ;
}

// RCPP_MODULE(exportmoddev){
//   using namespace Rcpp ;
//   class_<distribution>("distribution")
//     .constructor()
//   ;
//   class_<Logistic>("Logistic")
//     .derives<distribution>("distribution")
//     .constructor()
//     .method( "InverseLinkCumulativeFunction", &Logistic::InverseLinkCumulativeFunction )
//   ;
// }

