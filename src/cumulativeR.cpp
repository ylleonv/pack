#include "cumulativeR.h"
using namespace std;
using namespace Rcpp ;

CumulativeR::CumulativeR(void) {
  Rcout << "CumulativeR is being created" << endl;
}

Eigen::VectorXd CumulativeR::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_logit( eta(0) );
  for(size_t j=1; j<eta.size(); ++j)
  { ordered_pi[j] = Logistic::cdf_logit( eta(j) ) -  Logistic::cdf_logit( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(size_t j=0; j<eta.rows(); ++j)
  {F(j,j) =  Logistic::pdf_logit(eta(j));}
  return (F * R);
}

Eigen::VectorXd CumulativeR::inverse_probit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = Probit::cdf_probit( eta(0) );
  for(size_t j=1; j<eta.size(); ++j)
  { ordered_pi[j] = Probit::cdf_probit( eta(j) ) - Probit::cdf_probit( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_probit(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(size_t j=0; j<eta.rows(); ++j)
  { F(j,j) = Probit::cdf_probit( eta(j) ); }
  return (F * R);
}

Eigen::VectorXd CumulativeR::inverse_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = Cauchit::cdf_cauchit( eta(0) );
  for(size_t j=1; j<eta.size(); ++j)
  { ordered_pi[j] = Cauchit::cdf_cauchit( eta(j) ) - Cauchit::cdf_cauchit( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(size_t j=0; j<eta.rows(); ++j)
  { F(j,j) = Cauchit::cdf_cauchit( eta(j) ); }
  return (F * R);
}

Eigen::MatrixXd CumulativeR::GLMcum(Eigen::MatrixXd X_EXT, Eigen::VectorXd Y_EXT, int K, int P, int N, int Q, std::string link){
  // const int Q = K-1 ;
  // const int P = X_M.cols() -1 ;
  // const int N = X_M.n_cols ;
  Eigen::MatrixXd BETA = Eigen::MatrixXd::Zero((P+1)*Q,1);
  // Eigen::MatrixXd F_i;
  // double llikelihood1 = 0;

  for (int iteration=1; iteration < 40; iteration++){
    int i = 0;
    // Eigen::VectorXd llikelihood;
    // double llikelihood2;

    Eigen::MatrixXd X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
    Eigen::VectorXd Y_M_i = Y_EXT.segment(i*Q , Q);
    Eigen::VectorXd eta = X_M_i * BETA;

    Eigen::VectorXd pi ;
    if(link == "logistic"){
      pi = CumulativeR::inverse_logistic(eta);
    }else if(link == "probit"){
      pi = CumulativeR::inverse_probit(eta);
    }

    Eigen::MatrixXd Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());

    Eigen::MatrixXd D ;

    if(link == "logistic"){
      D = CumulativeR::inverse_derivative_logistic(eta);
    }else if(link == "probit"){
      D = CumulativeR::inverse_derivative_probit(eta);
    }

    Eigen::MatrixXd W_in = D * Cov_i.inverse();
    Eigen::MatrixXd Score_i = X_M_i.transpose() * W_in * (Y_M_i - pi);
    Eigen::MatrixXd F_i = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);

    for (i=1; i < N; i++){
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_EXT.segment(i*Q ,Q);
      eta = X_M_i * BETA;

      if(link == "logistic"){
        pi = CumulativeR::inverse_logistic(eta);
      }else if(link == "probit"){
        pi = CumulativeR::inverse_probit(eta);
      }


      Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());

      // D = CumulativeR::inverse_derivative_logistic(eta);

      if(link == "logistic"){
        D = CumulativeR::inverse_derivative_logistic(eta);
      }else if(link == "probit"){
        D = CumulativeR::inverse_derivative_probit(eta);
      }

      W_in = D * Cov_i.inverse();
      Eigen::MatrixXd Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      Eigen::MatrixXd F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
      F_i = F_i + F_i_2;
    }

    BETA = BETA + (F_i.inverse() * Score_i);

    //     cout << "ll" << endl;
    //     cout << llikelihood1 << endl;

    // cout << "Beta" << endl;
    // cout << BETA << endl;
  }
  // cout << "ll" << endl;
  // cout << llikelihood1 << endl;
  return BETA;
}


RCPP_MODULE(cumulativemodule){
  Rcpp::class_<CumulativeR>("CumulativeR")
  .constructor()
  .method( "GLMcum", &CumulativeR::GLMcum )
  ;
}




