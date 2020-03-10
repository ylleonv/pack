// #include "distribution.h"
// #include "adjacentR.h"
// using namespace std;
// using namespace Rcpp ;
//
// // [[Rcpp::depends(RcppEigen)]]
//
// AdjacentR::AdjacentR(void) {
//   cout << "Adjacent is being created" << endl;
// }
//
// Eigen::VectorXd AdjacentR::inverse_logistic(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi( eta.size() );
//   pi[eta.size()-1] = cdf_logit( eta(eta.size()-1) ) / ( 1-cdf_logit( eta(eta.size()-1) ) );
//   double norm = 1 + pi[eta.size()-1];
//   for(size_t j=(eta.size()-1); j>0; --j)
//   {
//     pi[j-1] = pi[j] * cdf_logit( eta(j-1) ) / ( 1-cdf_logit( eta(j-1) ) );
//     norm += pi[j-1];
//   }
//   return in_open_corner(pi/norm);
// }
//
// Eigen::MatrixXd AdjacentR::inverse_derivative_logistic(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi = AdjacentR::inverse_logistic(eta);
//   Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
//   Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
//   for(size_t j=0; j<pi.rows(); ++j)
//   { D(j,j) = pdf_logit( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_logit(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_logit(eta(j)))) ); }
//
//   return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );
//
// }
//
// Eigen::VectorXd AdjacentR::inverse_normal(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi( eta.size() );
//   pi[eta.size()-1] = cdf_normal( eta(eta.size()-1) ) / ( 1-cdf_normal( eta(eta.size()-1) ) );
//   double norm = 1 + pi[eta.size()-1];
//   for(size_t j=(eta.size()-1); j>0; --j)
//   {
//     pi[j-1] = pi[j] * cdf_normal( eta(j-1) ) / ( 1-cdf_normal( eta(j-1) ) );
//     norm += pi[j-1];
//   }
//   return in_open_corner(pi/norm);
// }
//
// Eigen::MatrixXd AdjacentR::inverse_derivative_normal(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi = AdjacentR::inverse_normal(eta);
//   Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
//   Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
//   for(size_t j=0; j<pi.rows(); ++j)
//   { D(j,j) = pdf_normal( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_normal(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_normal(eta(j)))) ); }
//
//   return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );
//
// }
//
// Eigen::VectorXd AdjacentR::inverse_cauchit(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi( eta.size() );
//   pi[eta.size()-1] = cdf_cauchit( eta(eta.size()-1) ) / ( 1-cdf_cauchit( eta(eta.size()-1) ) );
//   double norm = 1 + pi[eta.size()-1];
//   for(size_t j=(eta.size()-1); j>0; --j)
//   {
//     pi[j-1] = pi[j] * cdf_cauchit( eta(j-1) ) / ( 1-cdf_cauchit( eta(j-1) ) );
//     norm += pi[j-1];
//   }
//   return in_open_corner(pi/norm);
// }
//
// Eigen::MatrixXd AdjacentR::inverse_derivative_cauchit(const Eigen::VectorXd& eta) const
// {
//   Eigen::VectorXd pi = AdjacentR::inverse_cauchit(eta);
//   Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
//   Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
//   for(size_t j=0; j<pi.rows(); ++j)
//   { D(j,j) = pdf_cauchit( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_cauchit(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_cauchit(eta(j)))) ); }
//
//   return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );
//
// }
//
// List AdjacentR::GLMadj(std::string response,
//                         StringVector explanatory_complete,
//                         StringVector explanatory_proportional,
//                         std::string distribution,
//                         NumericVector categories_order,
//                         DataFrame dataframe){
//   // Number of explanatory variables
//   int P_c = explanatory_complete.size();
//   if(explanatory_complete[0] == "NA"){P_c = 0; }
//
//   int P_p = explanatory_proportional.size();
//   if(explanatory_proportional[0] == "NA"){P_p = 0; }
//   // const int P_c = explanatory_complete.size() ;
//   // const int P_p = explanatory_proportional.size() ;
//   int P =  P_c +  P_p ;
//   const int N = dataframe.nrows() ; // Number of observations
//   // Add Intercept acording to user input
//   Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(N);
//   dataframe["intercept"] = Ones1;
//
//   Eigen::MatrixXd Full_M = distribution::select_data(dataframe, response, explanatory_complete,
//                                                      explanatory_proportional, categories_order);
//   // Order by Y
//   Eigen::MatrixXd Full_M_ordered = distribution::sorted_rows(Full_M);
//
//   // Count number of unique categories (K) in Y
//   std::vector<int> Unique_k(Full_M_ordered.col(0).data(), Full_M_ordered.data() + Full_M_ordered.rows());
//   int K = std::set<int>( Unique_k.begin(), Unique_k.end() ).size();
//   int Q = K-1 ;
//   Eigen::VectorXd Y_ordered = Full_M_ordered.col(0);
//
//   // Expand Y to matrix where each column is a dummy variable indicating the category at which the subject row belongs:
//   // Matrix of zeros of size NxK
//   Eigen::MatrixXd Y_init2 = Eigen::MatrixXd::Zero(Full_M_ordered.rows(), K);
//   // Replace 0 by 1 at specific column
//   for (int i_1 = 0; i_1 < Y_init2.rows(); ++i_1){
//     Y_init2(i_1,Y_ordered[i_1]) = 1;
//   }
//   // Eliminate last column (Adjacent level will be always last category)
//   Eigen::MatrixXd Y_init = Y_init2.leftCols(Q);
//
//   // Create design matrix
//   Eigen::MatrixXd X_init = Full_M_ordered.rightCols(P);
//   Eigen::MatrixXd X_EXT(2, 2);
//
//   Eigen::MatrixXd X_M_Complete_Ext;
//   if(P_c > 0){
//     // X_M_Complete_Ext << X_init.block(0,0, N , P_c);
//     X_M_Complete_Ext = kroneckerProduct(X_init.block(0,0, N , P_c),Eigen::MatrixXd::Identity(Q,Q)).eval();
//   }
//
//
//   Eigen::MatrixXd X_M_Poportional_Ext(N*Q, P_p) ;
//   if(P_p > 0){
//     for (int x = -1; x < N-1; ++x) {
//       for (int j = (x+1)*Q ; j < (x+2)*Q; ++j){
//         X_M_Poportional_Ext.row(j) = (X_init.rightCols(P_p)).row(x+1);
//       }
//     }
//   }
//
//   X_EXT.conservativeResize( N*Q , X_M_Complete_Ext.cols()+X_M_Poportional_Ext.cols() );
//   X_EXT << X_M_Complete_Ext, X_M_Poportional_Ext;
//
//   // // Beta initialization
//   Eigen::MatrixXd BETA;
//   BETA = Eigen::MatrixXd::Zero(X_EXT.cols(),1);
//
//   int iteration = 0;
//   double check_tutz = 1.0;
//   double Stop_criteria = 1.0;
//   Eigen::MatrixXd X_M_i ;
//   Eigen::VectorXd Y_M_i ;
//   Eigen::VectorXd eta ;
//   Eigen::VectorXd pi ;
//   Eigen::MatrixXd D ;
//   Eigen::MatrixXd Cov_i ;
//   Eigen::MatrixXd W_in ;
//   Eigen::MatrixXd Score_i_2 ;
//   Eigen::MatrixXd F_i_2 ;
//   Eigen::VectorXd LogLikIter;
//   LogLikIter = Eigen::MatrixXd::Zero(1,1) ;
//   double LogLik;
//   // for (int iteration=1; iteration < 18; iteration++){
//   // while (check_tutz > 0.0001){
//   double epsilon = 0.0001 ;
//   while (Stop_criteria >( epsilon / N)){
//
//     Eigen::MatrixXd Score_i = Eigen::MatrixXd::Zero(BETA.rows(),1);
//     Eigen::MatrixXd F_i = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());
//
//     LogLik = 0.;
//
//     // Loop by subject
//
//     for (int i=0; i < N; i++){
//       // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
//       X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
//       Y_M_i = Y_init.row(i);
//       eta = X_M_i * BETA;
//
//       // Vector pi depends on selected distribution
//       if(distribution == "logistic"){
//         pi = AdjacentR::inverse_logistic(eta);
//         D = AdjacentR::inverse_derivative_logistic(eta);
//       }else if(distribution == "normal"){
//         pi = AdjacentR::inverse_normal(eta);
//         D = AdjacentR::inverse_derivative_normal(eta);
//       }else if(distribution == "cauchit"){
//         pi = AdjacentR::inverse_cauchit(eta);
//         D = AdjacentR::inverse_derivative_cauchit(eta);
//       }
//
//       Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
//       W_in = D * Cov_i.inverse();
//       Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
//       Score_i = Score_i + Score_i_2;
//       F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
//       F_i = F_i + F_i_2;
//       LogLik = LogLik + (Y_M_i.transpose().eval()*Eigen::VectorXd(pi.array().log())) + ( (1 - Y_M_i.sum()) * std::log(1 - pi.sum()) );
//     }
//
//     LogLikIter.conservativeResize(iteration+2, 1);
//     LogLikIter(iteration+1) = LogLik;
//     Stop_criteria = (abs(LogLikIter(iteration+1) - LogLikIter(iteration))) / (epsilon + (abs(LogLikIter(iteration+1)))) ;
//     Eigen::VectorXd beta_old = BETA;
//     BETA = BETA + (F_i.inverse() * Score_i);
//     // check_tutz = ((BETA - beta_old).norm())/(beta_old.norm()+check_tutz);
//     iteration = iteration + 1;
//   }
//   Rcout << "Number of iterations" << endl;
//   Rcout << iteration-1 << endl;
//   Rcout << "Log Likelihood" << endl;
//   Rcout << LogLik << endl;
//   Rcout << "Beta" << endl;
//   Rcout << BETA << endl;
//   return List::create(Named("Nb. iterations") = iteration-1 , Named("Coefficients") = BETA,
//                       Named("Log-likelihood") = LogLik);
// }
//
//
// RCPP_MODULE(adjacentmodule){
//   Rcpp::class_<AdjacentR>("AdjacentR")
//   .constructor()
//   .method( "GLMadj", &AdjacentR::GLMadj)
//   ;
// }
//
