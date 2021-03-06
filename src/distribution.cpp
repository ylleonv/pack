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

std::string distribution::concatenate(std::string x, std::string level)
{
  return (x + " " +level);
}

// ORDINAL DATA SELECTION

DataFrame my_transpose(DataFrame dat_in){
  Environment base_base("package:base");
  Function my_transpose1 = base_base["t"];
  DataFrame data_out = my_transpose1(dat_in);
  return data_out;
}

NumericMatrix to_dummy1(NumericVector A, CharacterVector levels)
{
  NumericVector cha_to_fact = A;
  CharacterVector levs1 = levels;
  int var_lev = levs1.length();

  // Rcout << var_lev << std::endl;

  NumericMatrix B_Ma(cha_to_fact.length(), var_lev);
  for (int i_1 = 0; i_1 < cha_to_fact.length(); ++i_1){
    int col_ind = cha_to_fact[i_1] - 1;
    B_Ma(i_1, col_ind) = 1;
  }

  // print(B_Ma);

  B_Ma = B_Ma( _ , Range(0,var_lev-2) );

  // if (var_lev != 2){
  //   B_Ma = B_Ma( _ , Range(0,var_lev-2) );
  // } else {
  //   B_Ma = B_Ma( _ , Range(1,var_lev-1) );
  // }
  return B_Ma;
}

List Model_Matrix_or(DataFrame data, Formula formula) {
  Environment stats_env("package:stats");
  Function model_frame = stats_env["model.frame"];
  Function model_matrix = stats_env["model.matrix"];
  DataFrame df_new1 = model_frame(Rcpp::_["formula"] = formula, _["data"] = data);
  SEXP Response = df_new1[0];
  NumericMatrix df_new = model_matrix(df_new1, _["data"] = data);

  return List::create(
    Named("df_new") = df_new,
    Named("Response") = Response
  );

}


List Cat_ref_order(CharacterVector categories_order, SEXP response_categories){
  // Environment base_env("package:base");
  // Function my_unique = base_env["unique"];

  CharacterVector response_categories1 = response_categories;
  // CharacterVector levels = my_unique(response_categories1);

  IntegerVector num_categories_order = seq_len(categories_order.length());

  DataFrame response_neworder = DataFrame::create(  _["num_categories_order"] = num_categories_order );

  response_neworder = my_transpose(response_neworder);
  response_neworder.names() = categories_order;

  String a1;
  CharacterVector a2;
  for(int i = 0; i <response_categories1.length(); i++){
    a1 = response_categories1[i];
    a2 = response_neworder[a1];
    response_categories1[i] = a2[0];
  }

  // print(response_categories1);

  return List::create(
    Named("response_neworder") = response_neworder,
    Named("response_categories2") = response_categories1,
    Named("levels") = categories_order
  );
}



List distribution::All_pre_data_or(Formula formula, DataFrame input_data,
                                   CharacterVector categories_order,
                                   CharacterVector proportional_effect,
                                   std::string threshold){

  Environment base_env("package:base");
  Function my_asnumeric = base_env["as.numeric"];
  Function my_cbind = base_env["cbind"];
  Function my_order = base_env["order"];

  List M_matrix = Model_Matrix_or(input_data, formula);
  List Cat_ref_or_L = Cat_ref_order(categories_order, M_matrix["Response"]);
  NumericMatrix Design = M_matrix["df_new"];

  CharacterVector a1 = Cat_ref_or_L["response_categories2"];

  NumericVector Num_res = my_asnumeric(a1);

  Design = my_cbind(Num_res, Design);

  // Now order dataset with respect to the repsonse variables in the given order
  DataFrame df_tans = my_transpose(Design);
  DataFrame df_tans_2 = Design ;
  NumericVector order_var_sel = my_order(df_tans_2[0]);
  order_var_sel = order_var_sel - 1 ;
  df_tans = df_tans[order_var_sel];

  df_tans_2 = my_transpose(df_tans);

  CharacterVector Levels = Cat_ref_or_L["levels"];


  int N_cats = Levels.length();

  LogicalVector any_alternative_specific = !is_na(proportional_effect); // TRUE IF THERE ARE
  CharacterVector colnames_final_m = df_tans_2.names();

  NumericVector x(colnames_final_m.length());
  if (any_alternative_specific[0]) {
    // x(colnames_final_m.length());
    for(int indi = 0 ; indi < proportional_effect.length(); indi++){
      String var_1 = proportional_effect[indi];
      int indi_var = df_tans_2.findName(var_1);
      x[indi_var] = indi_var;
    }
    colnames_final_m = colnames_final_m[x==0]; // Case where there no proportional effects
  }

  // Now extend
  // Y EXTEND
  NumericVector Response = df_tans_2[0];
  NumericMatrix Response_EXT = to_dummy1(Response, Levels);

  // X EXTEND

  // X COMPLETE

  DataFrame DF_complete_effect = df_tans_2[colnames_final_m];
  NumericMatrix Pre_Design1 = internal::convert_using_rfunction(DF_complete_effect, "as.matrix");
  Eigen::Map<Eigen::MatrixXd> Pre_Design = as<Eigen::Map<Eigen::MatrixXd> >(Pre_Design1);
  Eigen::MatrixXd Iden_Q1 = Eigen::MatrixXd::Identity(N_cats-1,N_cats-1);

  Eigen::MatrixXd X_EXT_COMPLETE;

  if (threshold == "equidistant"){ // eRASE THE LAST TWO COLUMNS ONE CORRESPONDING TO DR AND OTHER TO THE INTERCEPT
    X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 2), Iden_Q1).eval();
    colnames_final_m.erase(0);
    colnames_final_m.erase(0);
  } else {
    X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 1), Iden_Q1).eval();
    colnames_final_m.erase(0);
  }

  Eigen::MatrixXd Design_Matrix;

  if (any_alternative_specific[0]) {

    // PONER ESA PARTE ACA

    DataFrame DF_proportional_effect = df_tans_2[proportional_effect];
    NumericMatrix Pre_DF_proportional1 = internal::convert_using_rfunction(DF_proportional_effect, "as.matrix");
    Eigen::Map<Eigen::MatrixXd> Pre_DF_proportional2 = as<Eigen::Map<Eigen::MatrixXd> >(Pre_DF_proportional1);

    Eigen::MatrixXd Pre_DF_proportional = Pre_DF_proportional2;

    Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N_cats-1);
    Eigen::MatrixXd X_EXT_PROPORTIONAL = Eigen::kroneckerProduct(Pre_DF_proportional, Ones).eval();

    // TENGO QUE PONER EL IF ACA

    if (threshold == "equidistant"){
      NumericMatrix tJac = my_cbind(1, seq_len(categories_order.length() -1 )-1 );
      Eigen::Map<Eigen::MatrixXd> tJac2 = as<Eigen::Map<Eigen::MatrixXd> >(tJac);
      Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(Response_EXT.rows());
      Eigen::MatrixXd Threshold_M = Eigen::kroneckerProduct(Ones1, tJac2).eval();
      X_EXT_PROPORTIONAL.conservativeResize(X_EXT_PROPORTIONAL.rows(),X_EXT_PROPORTIONAL.cols()+2);
      X_EXT_PROPORTIONAL.block(0,X_EXT_PROPORTIONAL.cols()-2, X_EXT_PROPORTIONAL.rows(),2) = Threshold_M;
    }


    Design_Matrix.conservativeResize(X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()+X_EXT_PROPORTIONAL.cols());
    Design_Matrix.block(0,0,X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()) = X_EXT_COMPLETE;
    Design_Matrix.block(0,X_EXT_COMPLETE.cols(),X_EXT_COMPLETE.rows(),X_EXT_PROPORTIONAL.cols()) = X_EXT_PROPORTIONAL;

  }else{Design_Matrix = X_EXT_COMPLETE;}



  return List::create(
    Named("Design_Matrix") = Design_Matrix,
    Named("Response_EXT") = Response_EXT,
    Named("Levels") = Levels,
    Named("Complete_effects") = colnames_final_m,
    Named("N_cats") = N_cats
  );
}




// [[Rcpp::export]]
List All_pre_data_or2(Formula formula, DataFrame input_data,
                      CharacterVector categories_order,
                      CharacterVector proportional_effect){

  Environment base_env("package:base");
  Function my_asnumeric = base_env["as.numeric"];
  Function my_cbind = base_env["cbind"];
  Function my_order = base_env["order"];

  List M_matrix = Model_Matrix_or(input_data, formula);
  List Cat_ref_or_L = Cat_ref_order(categories_order, M_matrix["Response"]);
  NumericMatrix Design = M_matrix["df_new"];

  CharacterVector a1 = Cat_ref_or_L["response_categories2"];
  NumericVector Num_res = my_asnumeric(a1);

  Design = my_cbind(Num_res, Design);

  // Now order dataset with respect to the repsonse variables in the given order
  DataFrame df_tans = my_transpose(Design);
  DataFrame df_tans_2 = Design ;
  NumericVector order_var_sel = my_order(df_tans_2[0]);
  order_var_sel = order_var_sel - 1 ;
  df_tans = df_tans[order_var_sel];

  df_tans_2 = my_transpose(df_tans);

  CharacterVector Levels = Cat_ref_or_L["levels"];
  int N_cats = Levels.length();

  LogicalVector any_alternative_specific = !is_na(proportional_effect); // TRUE IF THERE ARE
  CharacterVector colnames_final_m = df_tans_2.names();

  NumericVector x(colnames_final_m.length());
  if (any_alternative_specific[0]) {
    // x(colnames_final_m.length());
    for(int indi = 0 ; indi < proportional_effect.length(); indi++){
      String var_1 = proportional_effect[indi];
      int indi_var = df_tans_2.findName(var_1);
      x[indi_var] = indi_var;
    }
    colnames_final_m = colnames_final_m[x==0]; // Case where there no proportional effects
  }

  // Now extend

  // X COMPLETE
  DataFrame DF_complete_effect = df_tans_2[colnames_final_m];
  NumericMatrix Pre_Design1 = internal::convert_using_rfunction(DF_complete_effect, "as.matrix");
  Eigen::Map<Eigen::MatrixXd> Pre_Design = as<Eigen::Map<Eigen::MatrixXd> >(Pre_Design1);
  Eigen::MatrixXd Iden_Q1 = Eigen::MatrixXd::Identity(N_cats-1,N_cats-1);
  Eigen::MatrixXd X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 1), Iden_Q1).eval();

  // X PROPOTIONAL

  Eigen::MatrixXd Design_Matrix;

  // Y EXTEND
  NumericVector Response = df_tans_2[0];
  NumericMatrix Response_EXT = to_dummy1(Response, Levels);

  NumericMatrix tJac = my_cbind(1, seq_len(categories_order.length() -1 )-1 );
  Rcout << tJac << endl;

  NumericMatrix Response_EXT1 = internal::convert_using_rfunction(Response_EXT, "as.matrix");
  NumericMatrix tJac1 = internal::convert_using_rfunction(tJac, "as.matrix");

  Eigen::Map<Eigen::MatrixXd> Response_EXT2 = as<Eigen::Map<Eigen::MatrixXd> >(Response_EXT1);
  Eigen::Map<Eigen::MatrixXd> tJac2 = as<Eigen::Map<Eigen::MatrixXd> >(tJac1);


  // Eigen::MatrixXd MAR = Response_EXT2 * tJac2;

  Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(Response_EXT2.rows());
  Eigen::MatrixXd MAR = Eigen::kroneckerProduct(Ones1, tJac2).eval();


  if (any_alternative_specific[0]) {

    DataFrame DF_proportional_effect = df_tans_2[proportional_effect];
    NumericMatrix Pre_DF_proportional1 = internal::convert_using_rfunction(DF_proportional_effect, "as.matrix");
    Eigen::Map<Eigen::MatrixXd> Pre_DF_proportional2 = as<Eigen::Map<Eigen::MatrixXd> >(Pre_DF_proportional1);

    Eigen::MatrixXd Pre_DF_proportional = Pre_DF_proportional2;

    // Pre_DF_proportional.conservativeResize(Pre_DF_proportional.rows(),Pre_DF_proportional.cols()+2);
    // Pre_DF_proportional.block(0,Pre_DF_proportional2.cols(),Pre_DF_proportional.rows(),2) = MAR;


    Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N_cats-1);
    Eigen::MatrixXd X_EXT_PROPORTIONAL = Eigen::kroneckerProduct(Pre_DF_proportional, Ones).eval();

    X_EXT_PROPORTIONAL.conservativeResize(X_EXT_PROPORTIONAL.rows(),X_EXT_PROPORTIONAL.cols()+2);
    X_EXT_PROPORTIONAL.block(0,X_EXT_PROPORTIONAL.cols()-2,X_EXT_PROPORTIONAL.rows(),2) = MAR;

    Design_Matrix.conservativeResize(X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()+X_EXT_PROPORTIONAL.cols());
    Design_Matrix.block(0,0,X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()) = X_EXT_COMPLETE;
    Design_Matrix.block(0,X_EXT_COMPLETE.cols(),X_EXT_COMPLETE.rows(),X_EXT_PROPORTIONAL.cols()) = X_EXT_PROPORTIONAL;

  }else{Design_Matrix = X_EXT_COMPLETE;}



  colnames_final_m.erase(0);

  return List::create(
    Named("Design_Matrix") = Design_Matrix,
    Named("Response_EXT") = Response_EXT,
    Named("Levels") = Levels,
    Named("Complete_effects") = colnames_final_m,
    Named("MAR2") = MAR,
    Named("N_cats") = N_cats
  );
}


List distribution::All_pre_data_NEWDATA(Formula formula,
                                        DataFrame NEWDATA,
                                        CharacterVector categories_order,
                                        CharacterVector proportional_effect,
                                        int N_cats){

  Environment base_env("package:base");
  Function my_asnumeric = base_env["as.numeric"];
  Function my_cbind = base_env["cbind"];
  Function my_order = base_env["order"];

  List M_matrix = Model_Matrix_or(NEWDATA, formula);
  // List Cat_ref_or_L = Cat_ref_order(categories_order, M_matrix["Response"]);
  NumericMatrix Design = M_matrix["df_new"];

  // CharacterVector a1 = Cat_ref_or_L["response_categories2"];
  // NumericVector Num_res = my_asnumeric(a1);
  //
  // Design = my_cbind(Num_res, Design);

  // Now order dataset with respect to the repsonse variables in the given order
  // DataFrame df_tans = my_transpose(Design);
  DataFrame df_tans_2 = Design ;


  // NumericVector order_var_sel = my_order(df_tans_2[0]);
  // order_var_sel = order_var_sel - 1 ;
  // df_tans = df_tans[order_var_sel];

  // df_tans_2 = my_transpose(df_tans);

  // CharacterVector Levels = Cat_ref_or_L["levels"];
  // int N_cats = Levels.length();

  LogicalVector any_alternative_specific = !is_na(proportional_effect); // TRUE IF THERE ARE
  CharacterVector colnames_final_m = df_tans_2.names();


  NumericVector x(colnames_final_m.length());
  if (any_alternative_specific[0]) {
    // x(colnames_final_m.length());
    for(int indi = 0 ; indi < proportional_effect.length(); indi++){
      String var_1 = proportional_effect[indi];
      int indi_var = df_tans_2.findName(var_1);
      x[indi_var] = indi_var;
    }
    colnames_final_m = colnames_final_m[x==0]; // Case where there no proportional effects
  }


  // X EXTEND

  // X COMPLETE

  DataFrame DF_complete_effect = df_tans_2[colnames_final_m];

  NumericMatrix Pre_Design1 = internal::convert_using_rfunction(DF_complete_effect, "as.matrix");
  Eigen::Map<Eigen::MatrixXd> Pre_Design = as<Eigen::Map<Eigen::MatrixXd> >(Pre_Design1);
  Eigen::MatrixXd Iden_Q1 = Eigen::MatrixXd::Identity(N_cats-1,N_cats-1);

  Eigen::MatrixXd X_EXT_COMPLETE;

  // if (threshold == "equidistant"){ // eRASE THE LAST TWO COLUMNS ONE CORRESPONDING TO DR AND OTHER TO THE INTERCEPT
  //   X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 2), Iden_Q1).eval();
  //   colnames_final_m.erase(0);
  //   colnames_final_m.erase(0);
  // } else {
  X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design, Iden_Q1).eval();


  // colnames_final_m.erase(0);
  // }

  Eigen::MatrixXd Design_Matrix;

  if (any_alternative_specific[0]) {

    // PONER ESA PARTE ACA

    DataFrame DF_proportional_effect = df_tans_2[proportional_effect];
    NumericMatrix Pre_DF_proportional1 = internal::convert_using_rfunction(DF_proportional_effect, "as.matrix");
    Eigen::Map<Eigen::MatrixXd> Pre_DF_proportional2 = as<Eigen::Map<Eigen::MatrixXd> >(Pre_DF_proportional1);

    Eigen::MatrixXd Pre_DF_proportional = Pre_DF_proportional2;

    Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N_cats-1);
    Eigen::MatrixXd X_EXT_PROPORTIONAL = Eigen::kroneckerProduct(Pre_DF_proportional, Ones).eval();

    // TENGO QUE PONER EL IF ACA
    //
    // if (threshold == "equidistant"){
    //   NumericMatrix tJac = my_cbind(1, seq_len(categories_order.length() -1 )-1 );
    //   Eigen::Map<Eigen::MatrixXd> tJac2 = as<Eigen::Map<Eigen::MatrixXd> >(tJac);
    //   Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(Response_EXT.rows());
    //   Eigen::MatrixXd Threshold_M = Eigen::kroneckerProduct(Ones1, tJac2).eval();
    //   X_EXT_PROPORTIONAL.conservativeResize(X_EXT_PROPORTIONAL.rows(),X_EXT_PROPORTIONAL.cols()+2);
    //   X_EXT_PROPORTIONAL.block(0,X_EXT_PROPORTIONAL.cols()-2, X_EXT_PROPORTIONAL.rows(),2) = Threshold_M;
    // }


    Design_Matrix.conservativeResize(X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()+X_EXT_PROPORTIONAL.cols());
    Design_Matrix.block(0,0,X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()) = X_EXT_COMPLETE;
    Design_Matrix.block(0,X_EXT_COMPLETE.cols(),X_EXT_COMPLETE.rows(),X_EXT_PROPORTIONAL.cols()) = X_EXT_PROPORTIONAL;

  }else{Design_Matrix = X_EXT_COMPLETE;}

  // Now extend


  // // X COMPLETE
  // DataFrame DF_complete_effect = df_tans_2[colnames_final_m];
  // NumericMatrix Pre_Design1 = internal::convert_using_rfunction(DF_complete_effect, "as.matrix");
  // Eigen::Map<Eigen::MatrixXd> Pre_Design = as<Eigen::Map<Eigen::MatrixXd> >(Pre_Design1);
  // // print(colnames_final_m);
  // Eigen::MatrixXd Iden_Q1 = Eigen::MatrixXd::Identity(N_cats-1,N_cats-1);
  //
  // Eigen::MatrixXd X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 1), Iden_Q1).eval();
  //
  // // X PROPOTIONAL
  //
  // Eigen::MatrixXd Design_Matrix;
  //
  // if (any_alternative_specific[0]) {
  //
  //   DataFrame DF_proportional_effect = df_tans_2[proportional_effect];
  //   NumericMatrix Pre_DF_proportional1 = internal::convert_using_rfunction(DF_proportional_effect, "as.matrix");
  //   Eigen::Map<Eigen::MatrixXd> Pre_DF_proportional = as<Eigen::Map<Eigen::MatrixXd> >(Pre_DF_proportional1);
  //   Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N_cats-1);
  //   Eigen::MatrixXd X_EXT_PROPORTIONAL = Eigen::kroneckerProduct(Pre_DF_proportional, Ones).eval();
  //
  //   Design_Matrix.conservativeResize(X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()+X_EXT_PROPORTIONAL.cols());
  //   Design_Matrix.block(0,0,X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()) = X_EXT_COMPLETE;
  //   Design_Matrix.block(0,X_EXT_COMPLETE.cols(),X_EXT_COMPLETE.rows(),X_EXT_PROPORTIONAL.cols()) = X_EXT_PROPORTIONAL;
  //
  // }else{Design_Matrix = X_EXT_COMPLETE;}

  return List::create(
    Named("Design") = Design,
    Named("Design_Matrix") = Design_Matrix
  );
}


CharacterVector Var_Not_In(DataFrame final_matrix, CharacterVector alternative_specific){

  LogicalVector any_alternative_specific = !is_na(alternative_specific);

  CharacterVector colnames_final_m = final_matrix.names();

  if (any_alternative_specific[0]) {
    NumericVector x(colnames_final_m.length());
    for(int indi = 0 ; indi < alternative_specific.length(); indi++){
      String var_1 = alternative_specific[indi];
      int indi_var = final_matrix.findName(var_1);
      x[indi_var] = indi_var;
    }
    colnames_final_m = colnames_final_m[x==0]; // Case where there no alternative specific variables
  }
  colnames_final_m.erase(0, 4); // Eliminate the first 4 information variables

  return colnames_final_m;
}

List formula_entry(Formula formula1){
  Environment base_base("package:base");
  Function my_strsplit = base_base["strsplit"];
  Function my_format = base_base["format"];
  Function my_paste = base_base["paste"];
  Function my_sub = base_base["sub"];
  Function my_rev = base_base["rev"];
  Function my_trimws = base_base["trimws"];

  Environment base_stats("package:stats");
  Function my_as_formula = base_stats["as.formula"];

  CharacterVector st1 = my_format(formula1);
  String str_for = my_paste(st1,  _["collapse"] = "");
  List list1 = (my_strsplit(str_for, "~"));
  CharacterVector vars = list1[0];
  String vars_string = vars[1];
  String res1 = vars[0];
  String res = my_trimws(res1);
  List list_vars = my_strsplit(vars_string, "[+]");
  CharacterVector char_vars = list_vars(0);

  String vars_for = my_paste(my_sub("\\[[^()]*\\]", "", char_vars),  _["collapse"] = " + ");
  CharacterVector form = my_paste(res, vars_for, _["sep"] = "~");
  Formula formula2 = my_as_formula(form);

  String firs_col;
  List list_no_spaces, list_var_rev, list_cat_rev;
  StringVector Vars(char_vars.length()), Alternatives(char_vars.length());

  for (int i = 0; i < char_vars.length() ; i++) {
    firs_col = char_vars[i];
    String no_spaces = my_trimws(firs_col);
    list_no_spaces = my_strsplit(no_spaces, "");
    String rev_string = my_paste(my_rev(list_no_spaces[0]), _["collapse"] = "");

    String var = my_sub(".*\\[", "", rev_string);
    list_var_rev = my_strsplit(var, "");
    String var1 = my_paste(my_rev(list_var_rev[0]), _["collapse"] = "");

    String cat = my_sub(".*\\]", "", my_sub("\\[.*", "", rev_string));
    list_cat_rev = my_strsplit(cat, "");
    String cat1 = my_paste(my_rev(list_cat_rev[0]), _["collapse"] = "");

    if (var1 != cat1){ Alternatives[i] = cat1; }
    if (var1 == "1"){ var1 = "(Intercept)"; }
    Vars[i] = var1;
  }

  DataFrame Var_alt = DataFrame::create(Named("Alternatives") = Alternatives);
  Var_alt = my_transpose(Var_alt);
  Var_alt.names() = Vars;

  DataFrame Var_alt1 = Var_alt;

  return List::create(
    Named("Var_alt") = Var_alt,
    Named("Response") = res,
    Named("formula_model") = formula2
  );
}

NumericVector Cat_ref(CharacterVector alternatives, String Ref_cat){
  NumericVector Ref_vec(alternatives.length());
  for (int i = 0; i < alternatives.length(); i++){
    if(alternatives(i) == Ref_cat){
      Ref_vec(i) = 1;
    }else{Ref_vec(i) = 0;}
  }
  return Ref_vec;
}

DataFrame Sort_DataFrame(DataFrame ModelMatrix,
                         DataFrame InputData,
                         CharacterVector names,
                         String Choice_vec,
                         String Ref_cat) {
  // Order DATAFRAME ACCORDING TO VARIABLES GIVEN IN VECTOR NAMES
  // CBIND OF DATA SETS AND THEN ORDER ACCORDING TO VARIABLES
  Environment base_env("package:base");
  Function my_order = base_env["order"];
  Function my_cbind = base_env["cbind"];
  Function my_asnumeric = base_env["as.numeric"];

  String alt = names[0];
  String id_case_0 = names[2];

  NumericVector Cat_ref_vec = Cat_ref(InputData[alt], Ref_cat);

  DataFrame A2 = my_cbind( _["alternatives"] = InputData[alt],
                           _["Cat_ref_vec"] = Cat_ref_vec,
                           _["id_case"] = InputData[id_case_0],
                                                   _["choice"] = my_asnumeric(InputData[Choice_vec]),
                                                   ModelMatrix);

  CharacterVector names1 = {"alternatives", "Cat_ref_vec", "id_case"};

  DataFrame df_tans = my_transpose(A2);
  DataFrame df_tans_2 = A2 ;

  for (int i = 0; i < names1.length() ; i++) {
    String var = names1(i);
    NumericVector order_var_sel = my_order(df_tans_2[var]);
    order_var_sel = order_var_sel - 1 ;
    df_tans = df_tans[order_var_sel];
    df_tans_2 = my_transpose(df_tans);
  }
  return  df_tans_2;
}


DataFrame my_AsNumericMatrix(DataFrame dat_in){
  Environment base_env("package:base");
  Function my_asnumeric = base_env["as.numeric"];
  Function my_ascharacter = base_env["as.character"];
  Function my_cbind = base_env["cbind"];
  DataFrame data_out = dat_in ;

  for (int i = 4; i < dat_in.length() ; i++) {
    NumericVector vec = my_asnumeric(my_ascharacter(dat_in[i]));
    data_out[i] = vec;
  }
  return data_out;
}

NumericMatrix Model_Matrix(DataFrame data, Formula formula) {
  Environment stats_env("package:stats");
  Function model_frame = stats_env["model.frame"];
  Function model_matrix = stats_env["model.matrix"];
  DataFrame df_new1 = model_frame(Rcpp::_["formula"] = formula, _["data"] = data);
  NumericMatrix df_new = model_matrix(df_new1, _["data"] = data);
  return  df_new;
}

DataFrame All_pre_data(Formula formula, DataFrame input_data, CharacterVector var_informatives, String choice, String Ref_cat){
  DataFrame data_output = my_AsNumericMatrix(Sort_DataFrame(
    Model_Matrix(input_data,
                 formula_entry(formula)["formula_model"]),
                 input_data,
                 var_informatives,
                 choice,
                 Ref_cat));
  return data_output;
}

Eigen::MatrixXd Extend_alt_specific(DataFrame alt_specific, int N_cats, int N_ind, CharacterVector var_alt_specific){
  // cat_index = 1;

  Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(N_cats-1);
  Eigen::MatrixXd Iden_Q = Eigen::MatrixXd::Identity(N_cats,N_cats);
  Iden_Q.conservativeResize(Iden_Q.rows() - 1, Iden_Q.cols());
  Iden_Q.col(N_cats-1) = -Ones1;
  NumericMatrix alt_specific_num = internal::convert_using_rfunction(alt_specific[var_alt_specific], "as.matrix");
  Eigen::Map<Eigen::MatrixXd> M_alt_specific = as<Eigen::Map<Eigen::MatrixXd> >(alt_specific_num);
  Eigen::MatrixXd Matrix_trans((N_cats-1)*N_ind,var_alt_specific.length());
  for(int indi = 1 ; indi <= N_ind ; indi++)
  {
    Eigen::MatrixXd Block_ind =  M_alt_specific.block((indi-1) * N_cats, 0, N_cats, var_alt_specific.length());
    Eigen::MatrixXd Block_RES = Block_ind.transpose() * Iden_Q.transpose();
    Matrix_trans.block((indi-1) * (N_cats-1), 0, N_cats-1, var_alt_specific.length()) = Block_RES.transpose();
  }

  return Matrix_trans;
}

Eigen::MatrixXd Extend_case_specific(DataFrame case_specific, int N_cats, int N_ind,
                                     CharacterVector var_alt_specific,  DataFrame Var_alt){

  CharacterVector var_case_specific = Var_Not_In(case_specific, var_alt_specific);
  DataFrame effect_specific_for1 = Var_alt[var_case_specific];
  effect_specific_for1 = my_transpose(effect_specific_for1);
  CharacterVector effect_specific_for2 = effect_specific_for1[0];

  Eigen::MatrixXd Iden_Q1 = Eigen::MatrixXd::Identity(N_cats-1,N_cats-1);
  NumericMatrix case_specific_num = internal::convert_using_rfunction(case_specific[var_case_specific], "as.matrix");
  Eigen::Map<Eigen::MatrixXd> M_case_specific = as<Eigen::Map<Eigen::MatrixXd> >(case_specific_num);

  Eigen::MatrixXd Matrix_trans((N_cats-1)*N_ind,(N_cats-1)*var_case_specific.length());

  for(int indi = 1 ; indi <= N_ind ; indi++)
  {
    Eigen::MatrixXd Block_ind =  M_case_specific.block((indi-1) * N_cats, 0, 1, var_case_specific.length());
    Eigen::MatrixXd Block_RES = Eigen::kroneckerProduct(Block_ind, Iden_Q1).eval();
    Matrix_trans.block((indi-1) * (N_cats-1), 0, N_cats-1, (N_cats-1)*var_case_specific.length()) = Block_RES;
  }
  // Rcout << Matrix_trans << std::endl;
  // Create vector relation between category and order in dataset
  CharacterVector ordered_cat = case_specific["alternatives"];
  String cat_1 = ordered_cat[0];
  NumericVector cat_index = NumericVector::create(Named(cat_1 , 0));

  for(int j_1 = 1; j_1 < N_cats ; j_1++){
    cat_1 = ordered_cat[j_1];
    cat_index.push_back(j_1, cat_1);
  }
  int count_re_var = 0;
  for(int ind_b_var = 0; ind_b_var < var_case_specific.length() ; ind_b_var++){

    Eigen::MatrixXd Block_cat = Matrix_trans.block(0,ind_b_var*(N_cats-1),Matrix_trans.rows(), N_cats-1);

    if( effect_specific_for2[ind_b_var] != "" ){
      String cat_loop = effect_specific_for2[ind_b_var];
      int Var_to_keep = cat_index[cat_loop];
      Matrix_trans.block(0,count_re_var,Matrix_trans.rows(), 1) = Block_cat.block(0, Var_to_keep, Matrix_trans.rows(),1);
      count_re_var = count_re_var+1;

    }else{

      Matrix_trans.block(0,count_re_var,Matrix_trans.rows(), N_cats-1) = Block_cat;
      count_re_var = count_re_var+N_cats-1;

    }
  }
  Eigen::MatrixXd Matrix_trans1 = Matrix_trans.block(0,0,Matrix_trans.rows(), count_re_var);
  return Matrix_trans1;
}

Eigen::MatrixXd Extend_All_design(DataFrame Final_mat, DataFrame Var_alt, CharacterVector var_alt_specific,
                                  int N_ind, int N_cats){

  LogicalVector any_alternative_specific = !is_na(var_alt_specific); // TRUE IF ANY

  CharacterVector var_case_specific = Var_Not_In(Final_mat, var_alt_specific);
  Eigen::MatrixXd Ex_case_M = Extend_case_specific(Final_mat, N_cats, N_ind, var_alt_specific, Var_alt );

  Eigen::MatrixXd Design_Matrix;

  if (any_alternative_specific[0]) {
    Eigen::MatrixXd Ex_alt_M = Extend_alt_specific(Final_mat, N_cats, N_ind, var_alt_specific);
    int rows_t = Ex_case_M.rows();
    Design_Matrix.conservativeResize(rows_t,Ex_case_M.cols()+Ex_alt_M.cols());
    Design_Matrix.block(0,0,rows_t,Ex_case_M.cols()) = Ex_case_M;
    Design_Matrix.block(0,Ex_case_M.cols(),rows_t,Ex_alt_M.cols()) = Ex_alt_M;
  }else{Design_Matrix = Ex_case_M;}

  return Design_Matrix;
}

List Extend_Response(DataFrame Final_mat ){

  Environment base_env("package:base");
  Function my_unique = base_env["unique"];
  Function my_length = base_env["length"];
  Function my_matrix = base_env["matrix"];
  Function my_asnumeric = base_env["as.numeric"];

  SEXP N_cat_1 = my_length(my_unique(Final_mat["alternatives"]));
  int N_cat = Rcpp::as<int>(N_cat_1);

  NumericVector y11 = my_asnumeric(Final_mat["choice"]) ;
  DataFrame Y_Ext = my_transpose(my_matrix(y11 ,  _["nrow"] = N_cat));
  NumericMatrix Y_Ext1 = internal::convert_using_rfunction(Y_Ext, "as.matrix");

  Eigen::MatrixXd Y_n2 = as<Eigen::Map<Eigen::MatrixXd> >(Y_Ext1-1);
  Y_n2.conservativeResize(Y_n2.rows(), Y_n2.cols() - 1);

  return List::create(
    Named("N_cat") = N_cat,
    Named("Y_Ext") = Y_n2
  );
}

List distribution::select_data_nested(Formula formula,
                                      String individuals,
                                      String Alternatives,
                                      SEXP ref_cat,
                                      CharacterVector var_alt_specific,
                                      DataFrame input_data
) {

  List Formula_l = formula_entry(formula);
  SEXP Var_spe_alt1 = Formula_l["Var_alt"];
  String Response = Formula_l["Response"];
  DataFrame Var_spe_alt = Rcpp::as<DataFrame>(Var_spe_alt1);

  CharacterVector var_informatives = {"Alternatives", "Cat_ref_vec", "individuals"};
  var_informatives[0] = Alternatives;
  var_informatives[2] = individuals;

  DataFrame Final_mat1 = All_pre_data(Formula_l["formula_model"],
                                      input_data,
                                      var_informatives,
                                      Response,
                                      ref_cat);

  List Response_L = Extend_Response(Final_mat1);
  Eigen::MatrixXd Response_M = Response_L["Y_Ext"];

  Eigen::MatrixXd Design_Matrix = Extend_All_design(Final_mat1,
                                                    Var_spe_alt,
                                                    var_alt_specific,
                                                    Response_M.rows(),
                                                    Response_L["N_cat"]);


  return List::create(
    _["Design_Matrix"] = Design_Matrix,
    _["Response_M"] = Response_M
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
}
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
Eigen::VectorXd Logistic::InverseLinkQuantileFunction(Eigen::VectorXd vector ){
  boost::math::logistic dist(0., 1.);
  for (int i = 0; i<=vector.size()-1; i++)
    vector(i) = quantile(dist, vector(i));
  return vector;
}


Normal::Normal(void) {
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
Eigen::VectorXd Normal::InverseLinkQuantileFunction(Eigen::VectorXd vector ){
  boost::math::normal norm;
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = quantile(norm, vector(i));
  return vector;
}

Cauchit::Cauchit(void) {
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
Eigen::VectorXd Cauchit::InverseLinkQuantileFunction(Eigen::VectorXd vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = quantile(extreme_value, vector(i));
  return vector;
}

Student::Student(void) {
}

double Student::cdf_student(const double& value, const double& freedom_degrees) const
{
  double z;
  if(freedom_degrees < 2 * pow(value, 2) )
  { z = boost::math::ibeta(freedom_degrees * 0.5, 0.5, freedom_degrees / (freedom_degrees + pow(value, 2))) * 0.5; }
  else
  { z = boost::math::ibetac(0.5, freedom_degrees * 0.5, pow(value, 2) / (freedom_degrees + pow(value, 2))) * 0.5; }
  if(value > 0)
  { return 1-z; }
  else
  {return z; }
}

double Student::pdf_student(const double& value, const double& freedom_degrees) const
{ return pow( freedom_degrees/(freedom_degrees + pow(value, 2)) , (1+freedom_degrees) * 0.5 ) / ( pow(freedom_degrees,0.5) * boost::math::beta(freedom_degrees*0.5, 0.5) ); }

Gumbel::Gumbel(void) {
  // Rcout << "Gumbel is being created" << endl;
}
double Gumbel::cdf_gumbel(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  return cdf(extreme_value, value);
}
double Gumbel::pdf_gumbel(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  return pdf(extreme_value, value);
}

Gompertz::Gompertz(void) {
  // Rcout << "Gompertz is being created" << endl;
}

double Gompertz::pdf_gompertz(const double& value) const
{ double _mu = 0.0;
  double _sigma = 1.0;

  return (exp((value - _mu)/ _sigma) *  exp( - exp ((value - _mu)/ _sigma) ) ) / _sigma ; }

double Gompertz::cdf_gompertz(const double& value) const
{ double _mu = 0.0;
  double _sigma = 1.0;
  return  1 - exp( - exp((value - _mu) / _sigma) ); }


RCPP_MODULE(exportmod){
  using namespace Rcpp ;
  class_<distribution>("distribution")
    .constructor()
  ;
  class_<Student>("Student")
    .derives<distribution>("distribution")
    .constructor()
    .method( "cdf_student", &Student::cdf_student )
    .method( "pdf_student", &Student::pdf_student )
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

