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
  // Rcout << "Distribution is being created" << endl;
}

LogicalVector is_character(DataFrame A) {
  LogicalVector res(A.cols());
  for (int column = 0 ; column < A.cols() ; column++){
    // bool a89 = (TYPEOF(A[column]) == STRSXP);
    // bool a90 = ( TYPEOF(A[column]) == STRSXP || TYPEOF(A[column]) ==  INTSXP );
    bool a90 = Rf_isFactor(A[column]) || (TYPEOF(A[column]) == STRSXP);
    res[column] = a90;
  }
  return res;
}

template <int RTYPE>
IntegerVector fast_factor_template( const Vector<RTYPE>& x ) {
  Vector<RTYPE> levs = sort_unique(x);
  IntegerVector out = match(x, levs);
  out.attr("levels") = as<CharacterVector>(levs);
  out.attr("class") = "integer";
  return out;
}

IntegerVector fast_factor( SEXP x ) {
  switch( TYPEOF(x) ) {
  case INTSXP: return fast_factor_template<INTSXP>(x);
  case REALSXP: return fast_factor_template<REALSXP>(x);
  case STRSXP: return fast_factor_template<STRSXP>(x);
  }
  return R_NilValue;
}

Eigen::VectorXd sort_vector_getindex(Eigen::VectorXd x1) {
  NumericVector V(x1.data(), x1.data() + x1.size());
  int x=0;
  std::iota(V.begin(),V.end(),x++);
  sort( V.begin(),V.end(), [&](int i,int j){return x1[i]<x1[j];} );
  Eigen::Map<Eigen::VectorXd> XS(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(V));
  return XS;
}

Eigen::MatrixXd sorted_rows(Eigen::MatrixXd A)
{
  Eigen::VectorXd vec1 = sort_vector_getindex(A.col(0));
  Eigen::MatrixXd B = A.row(vec1(0));
  for (int i = 1; i < A.rows(); ++i) {
    B.conservativeResize(B.rows()+1, B.cols());
    B.row(B.rows()-1) = A.row(vec1(i));
  }
  return B;
}

NumericMatrix to_dummy(SEXP A)
{
  IntegerVector cha_to_fact = fast_factor(A);
  CharacterVector levs1 = cha_to_fact.attr("levels");
  int var_lev = levs1.length();
  NumericMatrix B_Ma(cha_to_fact.length(), var_lev);
  for (int i_1 = 0; i_1 < cha_to_fact.length(); ++i_1){
    int col_ind = cha_to_fact[i_1] - 1;
    B_Ma(i_1, col_ind) = 1;
  }
  if (var_lev != 2){
    B_Ma = B_Ma( _ , Range(0,var_lev-2) );
  } else {
    B_Ma = B_Ma( _ , Range(1,var_lev-1) );
  }
  return B_Ma;
}

DataFrame sort_by_user(DataFrame A, SEXP order2)
{
  IntegerVector y_1 = (A[0]);
  StringVector y_n(A.rows());
  IntegerVector order = (fast_factor(order2));
  CharacterVector levs1 = order.attr("levels");
  CharacterVector levs2 = y_1.attr("levels");
  if (is_true(all(levs1 == levs2))) {
    for (int element_order = 0 ; element_order <= order.length(); element_order++){
      LogicalVector v0 = (y_1 == order[element_order]);
      y_n[v0] = element_order;
    }
    DataFrame B = A;
    CharacterVector y_2 = as<CharacterVector>(y_n);
    B[0] = y_2;
    return B;
  }stop("The response categories do not match the proposed order of entry");
}

std::string distribution::concatenate(std::string x, std::string level)
{
  return (x + " " +level);
}


List distribution::select_data(DataFrame x1, std::string response,
                               StringVector explanatory_complete,
                               StringVector explanatory_proportional,
                               SEXP order) {

  int P_c = explanatory_complete.size();
  if(explanatory_complete[0] == "NA"){P_c = 0; }
  // Rcout << P_c << std::endl;
  int P_p = explanatory_proportional.size();
  if(explanatory_proportional[0] == "NA"){P_p = 0; }
  // Rcout << P_p << std::endl;
  const int N = x1.nrows() ; // Number of observations

  // ADD INTERCEPT
  Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(x1.rows());
  x1["intercept"] = Ones1;

  // Zero initialization
  NumericVector a1(1);
  a1[0] = x1.findName(response);
  int n_com_cat = 0;
  for (int element = 0 ; element < explanatory_complete.size() ; element++ ){
    if(explanatory_complete[0] != "NA"){
      String element_1 = explanatory_complete[element];
      a1.push_back(x1.findName(element_1));
    }else {}
  }

  // SOLO CONTEO
  DataFrame x21 = x1[a1];
  LogicalVector n_com_cat1 = is_character(x21);
  n_com_cat = sum(n_com_cat1)-1;

  // CONTINUA PARA AÑADUR PROPORTIONAL
  for (int element_p = 0 ; element_p < explanatory_proportional.size() ; element_p++ ){
    if(explanatory_proportional[0] != "NA"){
      String element_2 = explanatory_proportional[element_p];
      a1.push_back(x1.findName(element_2));
    }else {}
  }

  // SE CREA LA MATRIZ COMPLETA DONDE Y ES LA PRIMERA COLUMNA
  DataFrame x23 = x1[a1];
  // Just assign factor vector to categorical order proposed
  DataFrame x2 = sort_by_user(x23, order); // solo reemplaza por nuevo orden

  LogicalVector character_var = is_character(x2);
  NumericMatrix X_com_cat_int( x2.nrows() , 1 ) ;

  // Solo para conteo de categoricas completas
  if(explanatory_complete[0] != "NA"){
    for (int column_char = 1 ; column_char < P_c+1; column_char++){
      if (character_var(column_char)){
        NumericMatrix a2 = to_dummy(x2[column_char]);
        X_com_cat_int = cbind(X_com_cat_int, a2);
      }else{
        NumericVector a3 = x2[column_char];
        X_com_cat_int = cbind(X_com_cat_int, a3);
      }
    }
  }
  int col_com = X_com_cat_int.cols() - 1;
  // Rcout << col_com << std::endl;

  NumericMatrix result_m( x2.nrows() , 1 ) ;
  for (int column_char1 = 0 ; column_char1 < x2.cols(); column_char1++){
    if(character_var[column_char1] == TRUE){
      NumericMatrix a8 = to_dummy(x2[column_char1]);
      result_m = cbind(result_m, a8);
    }
    else{
      NumericVector a9 = x2[column_char1];
      result_m = cbind(result_m, a9);
    }
  }

  IntegerVector a4 = fast_factor(x2[0]);
  CharacterVector levs1 = a4.attr("levels");
  int K = levs1.length();
  int Q = K-1;
  NumericVector a5 = as<NumericVector>(a4);
  result_m = result_m(_, Range(1,result_m.cols()-1));
  result_m = cbind(a5,result_m);
  Eigen::Map<Eigen::MatrixXd> P = as<Eigen::Map<Eigen::MatrixXd> >(result_m);
  P = sorted_rows(P);
  Eigen::MatrixXd P1 = P.rightCols(P.cols()-1);
  Eigen::MatrixXd Y_ext = P1.leftCols(K-1);
  Eigen::MatrixXd Com_ext = P1.block(0 , K-1 , P1.rows() , col_com );
  Eigen::MatrixXd Pro_ext = P1.rightCols(P1.cols() - Com_ext.cols() - Y_ext.cols()) ;

  Eigen::MatrixXd X_EXT(2, 2);
  Eigen::MatrixXd X_M_Complete_Ext;
  if(P_c > 0){
    X_M_Complete_Ext = Eigen::kroneckerProduct(Com_ext,Eigen::MatrixXd::Identity(Q,Q)).eval();
  }
  Eigen::MatrixXd X_M_Poportional_Ext(N*Q, Pro_ext.cols()) ;
  if(P_p > 0){
    for (int x = -1; x < N-1; ++x) {
      for (int j = (x+1)*Q ; j < (x+2)*Q; ++j){
        X_M_Poportional_Ext.row(j) = (Pro_ext).row(x+1);
      }
    }
  }
  X_EXT.conservativeResize( N*Q , X_M_Complete_Ext.cols()+X_M_Poportional_Ext.cols() );
  X_EXT << X_M_Complete_Ext, X_M_Poportional_Ext;

  return List::create(_["Y_ext"] = Y_ext,
                      _["X_EXT"] = X_EXT,
                      _["levs1"] = levs1);
}

List distribution::select_data_nested(DataFrame x1, std::string response, std::string actual_response,
                                      std::string individuals,
                                      StringVector explanatory_complete,
                                      StringVector depend_y,
                                      SEXP order) {

  const int N = x1.nrows() ;

  // ADD INTERCEPT
  Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(N);
  x1["intercept"] = Ones1;

  NumericVector a1(1);
  a1[0] = x1.findName(response);
  a1.push_back(x1.findName(individuals));
  a1.push_back(x1.findName(actual_response));

  for (int element = 0 ; element < explanatory_complete.size() ; element++ ){
    if(explanatory_complete[0] != "NA"){
      String element_1 = explanatory_complete[element];
      a1.push_back(x1.findName(element_1));
    }else {}
  }
  for (int element = 0 ; element < depend_y.size() ; element++ ){
    if(depend_y[0] != "NA"){
      String element_1 = depend_y[element];
      a1.push_back(x1.findName(element_1));
    }else {}
  }

  DataFrame x2 = sort_by_user(x1[a1], order); // solo reemplaza por nuevo orden
  // SOLO PARA CONTEO DE NIVELES RESPUESTA
  IntegerVector a4 = fast_factor(x2[0]);
  CharacterVector levs1 = a4.attr("levels");
  int K = levs1.length();
  int Q = K-1;
  NumericVector ref_cat(N) ;

  for (int i = 0 ; i < N; i++){
    if(a4[i] == K){
      ref_cat[i] = 1;
    }else {}
  }

  CharacterVector response_to_reemplace = x2[0];
  x2[0] = x2[1];
  x2[1] = response_to_reemplace; //ahora tienen los nombres opuestos

  LogicalVector choice_ind  = x2[2];

  NumericVector ind_response(N/K);
  int iter = 0;
  for (int i = 0; i < choice_ind.length(); i++) {
    if (choice_ind[i]) {
      ind_response[iter] = a4[i];
      iter++;
    } else {
    }
  }

  NumericMatrix Y_ext = to_dummy(ind_response);

  NumericVector a88(1);
  if(explanatory_complete[0] != "NA"){
    String element_1 = explanatory_complete[0];
    a88[0] = x2.findName(element_1);
    for (int element_p = 1 ; element_p < explanatory_complete.size() ; element_p++ ){
      String element_2 = explanatory_complete[element_p];
      a88.push_back(x2.findName(element_2));
    }
  }else {}
  NumericMatrix complete_var = internal::convert_using_rfunction(x2[a88], "as.matrix");

  NumericVector a89(1);
  if(depend_y[0] != "NA"){
    String element_1 = depend_y[0];
    a89[0] = x2.findName(element_1);
    for (int element_p = 1 ; element_p < depend_y.size() ; element_p++ ){
      String element_2 = depend_y[element_p];
      a89.push_back(x2.findName(element_2));
    }
  }else {}

  NumericMatrix to_change44 = internal::convert_using_rfunction(x2[a89], "as.matrix");
  to_change44 = cbind(to_change44, ref_cat);

  NumericVector to_change_sub(N);
  NumericMatrix vec;
  NumericVector vec1, vec_ref , ref5;

  // AHORA COMENZAMOS EL BUQLE POR INDIVIDUOS
  for(int vector = 0 ; vector < depend_y.length(); vector++){
    if(depend_y[0] != "NA"){
      for(int indi = 1 ; indi <= (N/K) ; indi++)
      {
        vec = to_change44( Range(K*(indi-1),(indi*K)-1) , _);
        vec1 = vec( _ , vector);
        vec_ref = vec( _ , depend_y.size());
        ref5 = vec1[vec_ref == 1];
        for (int y = 1; y <= K; ++y)
        {
          to_change_sub[K*(indi-1)+y-1] = vec1[y-1] - ref5[0];
        }
      }
      to_change44 = cbind(to_change44, to_change_sub);
    }else {}
  }
  to_change44 = to_change44(_, Range(to_change44.cols()-depend_y.length(), to_change44.cols() - 1));
  to_change44 = cbind(complete_var, to_change44);

  // ELIMINO LA FILA DE REFERENCIA
  NumericMatrix x52(Dimension(N- N/K, to_change44.ncol()));
  NumericMatrix x_to_ext(Dimension(N/K, to_change44.ncol()-depend_y.length()));
  int iter2 = 0;
  int def = 0;
  for (int i = 0; i < to_change44.nrow(); i++) {
    if (ref_cat[i] != 1) {
      x52.row(iter2) = to_change44.row(i);
      iter2++;
    } else {
      x_to_ext.row(def) = to_change44.row(i);
      def++;
    }
  }
  Eigen::Map<Eigen::MatrixXd> P = as<Eigen::Map<Eigen::MatrixXd> >(x_to_ext);
  Eigen::Map<Eigen::MatrixXd> ext_dep_y = as<Eigen::Map<Eigen::MatrixXd> >(x52);

  Eigen::MatrixXd X_M_Complete_Ext = Eigen::kroneckerProduct(P, Eigen::MatrixXd::Identity(Q,Q)).eval();




  // Eigen::MatrixXd X1 = X_M_Complete_Ext.block(0,0,X_M_Complete_Ext.rows(),Q+1) ;
  //
  // Eigen::MatrixXd X2 = X_M_Complete_Ext.block(0,(2*Q),X_M_Complete_Ext.rows(),(1)) ;



  Eigen::MatrixXd X1 = X_M_Complete_Ext.block(0,0,X_M_Complete_Ext.rows(),Q) ;
  Eigen::MatrixXd X2(X_M_Complete_Ext.rows(),1);
  // Eigen::MatrixXd X_mod = X1;
  //
  //
  if(explanatory_complete.size() >= 2){

    for (int element_p1 = 0 ; element_p1 < explanatory_complete.size() - 1; element_p1++ ){

      X2 = X_M_Complete_Ext.block(0,((element_p1+1)*Q),X_M_Complete_Ext.rows(),1) ;


      X1.conservativeResize(X1.rows(), X1.cols() + 1);

      X1.col(X1.cols()-1) = X2;

      // X_mod << X_mod, X2;
    }

  }else{}

  Eigen::MatrixXd X_M_dep_y(X_M_Complete_Ext.rows(),X_M_Complete_Ext.cols()+ depend_y.length());

  X_M_dep_y << X_M_Complete_Ext, ext_dep_y.rightCols(depend_y.length());

  Eigen::MatrixXd X_M_dep_y_alt(X1.rows(),X1.cols()+ depend_y.length());

  X_M_dep_y_alt << X1, ext_dep_y.rightCols(depend_y.length());

  return List::create( _["X1"] = X1,
                       _["X_M_dep_y"] = X_M_dep_y,
                       _["X_M_dep_y_alt"] = X_M_dep_y_alt,
                       _["X_M_Complete_Ext"] = X_M_Complete_Ext,
                       _["Y_ext"] = Y_ext,
                       _["levs1"] = levs1
  );
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
  // Rcout << "Logistic is being created" << endl;
}

Normal::Normal(void) {
  // Rcout << "normal is being created" << endl;
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

double Normal::cdf_normal(const double& value) const
{
  boost::math::normal norm;
  return boost::math::cdf(norm, value);
}

double Normal::pdf_normal(const double& value) const
{
  boost::math::normal norm;
  return boost::math::pdf(norm, value);
}


Eigen::VectorXd Normal::InverseLinkCumulativeFunction(Eigen::VectorXd vector ){
  boost::math::normal norm;
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = cdf(norm, vector(i));
  return vector;
}
Eigen::VectorXd Normal::InverseLinkDensityFunction(Eigen::VectorXd vector ){
  boost::math::normal norm;
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = pdf(norm, vector(i));
  return vector;
}

Cauchit::Cauchit(void) {
  // Rcout << "Cauchit is being created" << endl;
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
  // Rcout << "Student is being created" << endl;
}
double Student::cdf_student(const double& value) const
{
  double _degrees = 1.35;
  boost::math::students_t_distribution<> student(_degrees);
  return cdf(student, value);
}
double Student::pdf_student(const double& value) const
{
  double _degrees = 1.35;
  boost::math::students_t_distribution<> student(_degrees);
  return pdf(student, value);
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
  // Rcout << "Gumbel is being created" << endl;
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
  // Rcout << "Gompertz is being created" << endl;
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

double Gompertz::cdf_gompertz(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  return cdf(extreme_value, value);
}
double Gompertz::pdf_gompertz(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  return pdf(extreme_value, value);
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

