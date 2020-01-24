#include "distribution.h"
#include "reference.h"
using namespace std;
using namespace Rcpp ;

ReferenceF::ReferenceF(void) {
  Rcout << "Reference is being created" << endl;
}

Eigen::VectorXd ReferenceF::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi1(eta.size());
  double norm1 = 1.;
  for(int j=0; j<eta.size(); ++j)
  {
    double number ;
    number = Logistic::cdf_logit(eta(j));
    pi1[j] = number / ( 1-number );
    norm1 += pi1[j];
  }
  return (pi1/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_logistic(const Eigen::VectorXd& eta2) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_logistic(eta2);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  { D1(j,j) = pdf_logit( eta2(j) ) /
    (Logistic::cdf_logit(eta2(j)) * (1-Logistic::cdf_logit(eta2(j))));
  }
  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}

Eigen::VectorXd ReferenceF::inverse_probit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(size_t j=0; j<eta.size(); ++j)
  {
    pi[j] = Probit::cdf_probit( eta(j) ) / ( 1-Probit::cdf_probit( eta(j) ) );
    norm1 += pi[j];
  }
  return in_open_corner(pi/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_probit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = ReferenceF::inverse_probit(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  for(size_t j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_probit( eta(j) ) /
    ( std::max(1e-10, std::min(1-1e-6,Probit::cdf_probit(eta(j)))) *
      std::max(1e-10, std::min(1-1e-6, 1-Probit::cdf_probit(eta(j)))) ); }
  return D * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose().eval() );
}

// Eigen::VectorXd ReferenceF::lnatural(Eigen::VectorXd vector){
//   for (size_t i = 0; i<=vector.size(); i++)
//     vector[i] = std::log(vector[i]);
//   return vector;
// }


Eigen::MatrixXd ReferenceF::GLMref(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q, std::string link){
  // const int Q = K-1 ;
  // const int P = X_M.cols() -1 ;
  // const int N = X_M.n_cols ;
  Eigen::MatrixXd BETA = Eigen::MatrixXd::Zero((P+1)*Q,1);
  int iteration = 1;
  double check_tutz = 1.0;

  // for (int iteration=1; iteration < 40; iteration++){
  while (check_tutz > 1e-6){
    Eigen::VectorXd ll_vector;
    int i = 0;
    Eigen::MatrixXd X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
    Eigen::VectorXd Y_M_i = Y_EXT.segment(i*Q , Q);
    Eigen::VectorXd eta = X_M_i * BETA;
    Eigen::VectorXd pi ;
    if(link == "logistic"){
      pi = ReferenceF::inverse_logistic(eta);
    }else if(link == "probit"){
      pi = ReferenceF::inverse_probit(eta);
    }
    Eigen::MatrixXd Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
    Eigen::MatrixXd D ;
    if(link == "logistic"){
      D = ReferenceF::inverse_derivative_logistic(eta);
    }else if(link == "probit"){
      D = ReferenceF::inverse_derivative_probit(eta);
    }
    Eigen::MatrixXd W_in = D * Cov_i.inverse();
    Eigen::MatrixXd Score_i = X_M_i.transpose() * W_in * (Y_M_i - pi);
    Eigen::MatrixXd F_i = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
    // double llikelihood_prev = (Y_M_i.transpose()*ReferenceF::lnatural(pi))  ;
    // double llikelihood3 = (1-(Y_M_i.sum()*std::log( 1 - pi.sum())));
    // llikelihood_prev =  llikelihood_prev + llikelihood3;

    for (i=1; i < N; i++){
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_EXT.segment(i*Q ,Q);
      eta = X_M_i * BETA;
      if(link == "logistic"){
        pi = ReferenceF::inverse_logistic(eta);
      }else if(link == "probit"){
        pi = ReferenceF::inverse_probit(eta);
      }
      Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
      if(link == "logistic"){
        D = ReferenceF::inverse_derivative_logistic(eta);
      }else if(link == "probit"){
        D = ReferenceF::inverse_derivative_probit(eta);
      }
      W_in = D * Cov_i.inverse();
      Eigen::MatrixXd Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      Eigen::MatrixXd F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
      F_i = F_i + F_i_2;
      // double llikelihood1 = (Y_M_i.transpose()*ReferenceF::lnatural(pi))  ;
      // llikelihood3 = (1-(Y_M_i.sum()*std::log( 1 - pi.sum())));
      // llikelihood1 =  llikelihood1 + llikelihood3;
      // llikelihood_prev = llikelihood_prev + llikelihood1;
    }
    // Stop criteria Tutz
    Eigen::VectorXd beta_old = BETA;
    BETA = BETA + (F_i.inverse() * Score_i);
    check_tutz = ((BETA - beta_old).norm())/(beta_old.norm());
    // ll_vector.resize(iteration);
    // ll_vector(iteration-1) = (double) llikelihood_prev;
    // Rcout << "Log Likelihood" << endl;
    // Rcout << ll_vector(iteration-1) << endl;
    iteration = iteration + 1;
  }
  Rcout << "Number of iterations" << endl;
  Rcout << iteration-2 << endl;
  return BETA;
}

RCPP_MODULE(referencemodule){
  Rcpp::class_<ReferenceF>("ReferenceF")
  .constructor()
  .method( "GLMref", &ReferenceF::GLMref )
  ;
}




