#ifndef FISHER_H_
#define FISHER_H_
#include "distribution.h"

class FisherScoring : public Logistic, Probit, Cauchit, Student, Gumbel, Gompertz{
public:
  FisherScoring();
  const arma::mat X_M;
  const arma::mat Y_M;
  std::string link;
  arma::mat GLMm(const arma::mat X_M, const arma::mat Y_M, std::string link);
};

#endif


