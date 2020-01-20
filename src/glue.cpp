#include <RcppArmadillo.h>
#include "distribution.h"
#include "fisher.h"
#include "reference.h"
#include "cumulativeR.h"
#include "adjacentR.h"
#include "sequentialR.h"

using namespace Rcpp;

RCPP_MODULE(distribution){
  class_<distribution>("distribution")
  .constructor()
  ;
  class_<Logistic>("Logistic")
    .derives<distribution>("distribution")
    .constructor()
    .method( "InverseLinkCumulativeFunction", &Logistic::InverseLinkCumulativeFunction )
  ;
  class_<Probit>("Probit")
    .derives<distribution>("distribution")
    .constructor()
    .method( "InverseLinkCumulativeFunction", &Probit::InverseLinkCumulativeFunction )
  ;
  class_<Cauchit>("Cauchit")
    .derives<distribution>("distribution")
    .constructor()
    .method( "InverseLinkCumulativeFunction", &Cauchit::InverseLinkCumulativeFunction )
  ;
  class_<Student>("Student")
    .derives<distribution>("distribution")
    .constructor()
    .method( "InverseLinkCumulativeFunction", &Student::InverseLinkCumulativeFunction )
  ;
  class_<Gumbel>("Gumbel")
    .derives<distribution>("distribution")
    .constructor()
    .method( "InverseLinkCumulativeFunction", &Gumbel::InverseLinkCumulativeFunction )
  ;
  class_<Gompertz>("Gompertz")
    .derives<distribution>("distribution")
    .constructor()
    .method( "InverseLinkCumulativeFunction", &Gompertz::InverseLinkCumulativeFunction )
  ;
  class_<FisherScoring>("FisherScoring")
    // .derives<Logistic>("Logistic")
    .constructor()
    .method( "GLMm", &FisherScoring::GLMm )
  ;
  class_<ReferenceF>("ReferenceF")
    .constructor()
    .method( "GLMref", &ReferenceF::GLMref )
  ;
  class_<CumulativeR>("CumulativeR")
    .constructor()
    .method( "GLMcum", &CumulativeR::GLMcum )
  ;
  // class_<AdjacentR>("AdjacentR")
  //   .constructor()
  //   .method( "GLMadj", &AdjacentR::GLMadj)
  // ;
  class_<SequentialR>("SequentialR")
    .constructor()
    .method( "GLMseq", &SequentialR::GLMseq)
  ;
}

