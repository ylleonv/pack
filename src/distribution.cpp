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

Eigen::VectorXd distribution::sort_vector(Eigen::VectorXd x1) {
  NumericVector V(x1.data(), x1.data() + x1.size());
  int x=0;
  std::iota(V.begin(),V.end(),x++);
  sort( V.begin(),V.end(), [&](int i,int j){return x1[i]<x1[j];} );
  Eigen::Map<Eigen::VectorXd> XS(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(V));
  return XS;
}

DataFrame distribution::sort_by_user(DataFrame A, NumericVector order)
{
  NumericVector y_1 = A[0];
  NumericVector y_n(y_1.length());
  for (int element_order = 0 ; element_order < order.length(); element_order++){
    LogicalVector v0 = (y_1 == order[element_order]);
    y_n[v0] = element_order;
  }
  A[0] = y_n;
  DataFrame B = A;
  return B;
}

Eigen::MatrixXd distribution::sorted_rows(Eigen::MatrixXd A)
{
  Eigen::VectorXd vec1 = distribution::sort_vector(A.col(0));
  Eigen::MatrixXd B = A.row(vec1(0));
  for (int i = 1; i < A.rows(); ++i) {
    B.conservativeResize(B.rows()+1, B.cols());
    B.row(B.rows()-1) = A.row(vec1(i));
  }
  return B;
}

Eigen::MatrixXd distribution::select_data(DataFrame x1, std::string response,
                            StringVector explanatory_complete,
                            StringVector explanatory_proportional,
                            NumericVector order) {
  // Zero initialization
  NumericVector a1(1);
  a1[0] = x1.findName(response);
  for (int element = 0 ; element < explanatory_complete.size() ; element++ ){
    String element_1 = explanatory_complete[element];
    a1.push_back(x1.findName(element_1));
  }
  for (int element = 0 ; element < explanatory_proportional.size() ; element++ ){
    String element_1 = explanatory_proportional[element];
    a1.push_back(x1.findName(element_1));
  }
  x1 = distribution::sort_by_user(x1[a1], order);
  NumericMatrix x2 = internal::convert_using_rfunction(x1, "as.matrix");
  Eigen::Map<Eigen::MatrixXd> P = as<Eigen::Map<Eigen::MatrixXd> >(x2);
  return P;
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
  return boost::math::cdf(norm, value);
}


double Probit::pdf_probit(const double& value) const
{
  boost::math::normal norm;
  return boost::math::pdf(norm, value);
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
double Cauchit::pdf_cauchit(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  return pdf(extreme_value, value);
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
  double _degrees = 2.0;
  boost::math::students_t_distribution<> student(_degrees);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = cdf(student, vector(i));
  return vector;
}
Eigen::VectorXd Student::InverseLinkDensityFunction(Eigen::VectorXd vector ){
  double _degrees = 2.0;
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

Eigen::VectorXd sort_vector(Eigen::VectorXd x1) {
  NumericVector V(x1.data(), x1.data() + x1.size());
  int x=0;
  std::iota(V.begin(),V.end(),x++);
  sort( V.begin(),V.end(), [&](int i,int j){return x1[i]<x1[j];} );
  Eigen::Map<Eigen::VectorXd> XS(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(V));
  return XS;
}

Eigen::MatrixXd sorted_rows(Eigen::MatrixXd A)
{
  Eigen::VectorXd vec1 = sort_vector(A.col(0));
  Eigen::MatrixXd B = A.row(vec1(0));
  for (int i = 1; i < A.rows(); ++i) {
    B.conservativeResize(B.rows()+1, B.cols());
    B.row(B.rows()-1) = A.row(vec1(i));
  }
  return B;
}

RCPP_MODULE(exportmod){
  using namespace Rcpp ;
  class_<distribution>("distribution")
    .constructor()
  ;
}

RCPP_MODULE(exportmoddev){
  using namespace Rcpp ;
  class_<distribution>("distribution")
    .constructor()
  ;
  class_<Logistic>("Logistic")
    .derives<distribution>("distribution")
    .constructor()
    .method( "InverseLinkCumulativeFunction", &Logistic::InverseLinkCumulativeFunction )
  ;
}

