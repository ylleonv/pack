#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace std;
using namespace Rcpp ;
using Eigen::MatrixXd;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

// logit_model <- function(X, b){exp(X%*%b)/(1+sum(exp(X%*%b)))}

// [[Rcpp::export]]
arma::mat LogisticMulti(arma::mat a, arma::mat b)
{
  arma::mat r;
  r = exp(a * b);
  r.each_row() /= (1+sum(r));
  return r;
}

// [[Rcpp::export]]
arma::mat MultinomialRegressionAlgo(arma::mat X_M, arma::mat Y_EXT,
                                    int K){
  const int Q = K-1 ;
  const int P = X_M.n_cols -1 ;
  const int N = X_M.n_cols ;
  arma::mat Id_Q = arma::eye(Q,Q);
  arma::mat X_EXT = kron(X_M, Id_Q );
  arma::mat BETA((P+1)*Q, 1);
  BETA.fill(0.0) ; // initialize betas to 0
  arma::mat block_i;
  arma::mat W_G_A;
  Eigen::MatrixXd block_i_l_e;
  arma::mat pi;
  arma::mat vec_i;
  arma::mat block_i_l ;


  for (int iteration=1; iteration < 10; iteration++){

    // int i=0;
    // arma::mat pi = LogisticMulti(X_EXT.rows(i*Q, (i*Q)+Q-1), BETA);
    //
    // arma::mat vec_i = pi.rows(((i*Q)), ((i*Q)+Q-1));

    // block_i = arma::diagmat(vec_i) - (vec_i * vec_i.t());
    Eigen::MatrixXf W_G = Eigen::MatrixXf::Identity(N*Q, N*Q);

    for (int i=0; i < N-1; i++){
      pi.rows(((i*Q)), ((i*Q)+Q-1)) = LogisticMulti(X_EXT.rows(i*Q, (i*Q)+Q-1), BETA);
      vec_i = pi.rows(((i*Q)), ((i*Q)+Q-1));
      block_i_l = arma::diagmat(vec_i) - (vec_i * vec_i.t());

      // block_i_l_e = Eigen::Map<Eigen::MatrixXd>(block_i_l.memptr(),
      //                                                       block_i_l.n_rows,
      //                                                       block_i_l.n_cols);

      // W_G.block<block_i_l_e.rows(),block_i_l_e.rows()>(0,0) = block_i;
      //
      // W_G_A = arma::mat(W_G.data(), W_G.rows(), W_G.cols(),
      //                              false, false);

      // block_i <- adiag(block_i, block_i_l)

    }

    //   for (i in 1:(nrow(data2)-1)){
    //     pi[((i*Q)+1):((i*Q)+Q)] = logit_model(X_EXT[((i*Q)+1):((i*Q)+Q),],BETA)
    //     vec_i <- pi[((i*Q)+1):((i*Q)+Q)]
    //     block_i_l <- diag(vec_i) - (as.matrix(vec_i) %*% t(as.matrix(vec_i)))
    //     block_i <- adiag(block_i, block_i_l)
    //   }
    //   Score <- t(X_EXT) %*% (Y_EXT - pi)
    //     F_matrix <- (t(X_EXT) %*% block_i) %*% (X_EXT)
    //
    //     BETA = BETA + (solve(F_matrix) %*% Score)
    //     print(BETA)
  }

  return block_i_l;

}
