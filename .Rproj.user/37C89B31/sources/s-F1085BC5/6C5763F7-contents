# library(devtools);library(tidyverse);library(fs);library(Rcpp);library(RcppArmadillo);library(RcppEigen)
# library(dplyr); library(tidyr); library(nnet); library(varhandle);library(catdata); library(magic)
# library(VGAM)
#
# # Initial configuration ---------------------------------------------------
#
# load_all()
# check()
# use_mit_license("Lorena LEON")
# check()
# document()
#
# usethis::use_rcpp()
# # NO funciona timestwo y he creado el paquete siguiendo los pasos del libro
# timesTwo(5)
# # Cuando hago esto por segunda vez, funciona.
# usethis::use_rcpp()
#
# use_package("forcats")
# use_package("RcppArmadillo",type = "LinkingTo")
# use_package("RcppEigen",type = "LinkingTo")
# # Problema de rcpparmadillo not include solo se cambia el orden como sigue:
# # LinkingTo:
# # Rcpp,
# # RcppArmadillo,
# # RcppEigen
# use_package("RcppEigen",type = "Imports")
# use_package("RcppArmadillo",type = "Imports")
# #error "The file 'Rcpp.h' should not be included. Please correct to include only 'RcppArmadillo.h'."
# # Borre del distribution.h el include rcpp
# getLoadedDLLs()
# # Using modules
# # We need to tell R that we want to use a C++11 compiler
# Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# pkgbuild::compile_dll()
# compileAttributes()
#
#
# # Binomial case -----------------------------------------------------------
#
# # # DATA
# {
#   # Libraries
#   library(dobson)
#
#   # Data
#   beetle
#   i=1
#   beetle_ext <- data.frame(x = as.matrix(rep(beetle[i,1],beetle[i,2])), y= c(rep(1,beetle[i,3]), rep(0,beetle[i,2]-beetle[i,3]) ))
#
#   for (i in 2:nrow(beetle)) {
#     beetle_ext <- rbind(beetle_ext,data.frame(x = as.matrix(rep(beetle[i,1],beetle[i,2])), y= c(rep(1,beetle[i,3]), rep(0,beetle[i,2]-beetle[i,3]))))
#   }
#   beetle_ext <- as.data.frame(beetle_ext)
#   names(beetle_ext) <- c("x", "y")
#   head(beetle_ext)
#
#   # Matrix and vectors
#   X = as.matrix(data.frame(x0 = as.vector(rep(1,nrow(beetle_ext))), x1 = as.vector(unlist(as.vector(beetle_ext$x)))))
#   Y = beetle_ext$y
#   Y = matrix(Y)
#   beta = as.matrix(rep(0,2))
#   mu = as.matrix(X)%*%beta
# }
#
# # Pruebas con librerìa usual
#
# {
#   y=c(6,13,18,28,52,53,61,60)
#   n=c(59,60,62,56,63,59,62,60)
#   x=c(1.6907,1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)
#   n_y=n-y
#   beetle.mat=cbind(y,n_y)
#   glm(beetle.mat~x, family = binomial(link = "logit"))
#   glm(beetle.mat~x, family = binomial(link = "probit"))
#   glm(beetle.mat~x, family = binomial(link = "cauchit"))
#   glm(beetle.mat~x, family = binomial(link = "cloglog"))
# }
#
# # Logit
# FisherScoring
# dist3 <- new(FisherScoring)
# a3 = dist3$GLMm(X_M = X,
#                 Y_M = Y,
#                 link = "logistic")
#
# # Probit
# FisherScoring
# dist3 <- new(FisherScoring)
# a3 = dist3$GLMm(X_M = X,
#                 Y_M = Y,
#                 link = "probit")
#
# # Cauchit
# FisherScoring
# dist3 <- new(FisherScoring)
# dist3$GLMm(X_M = X,
#            Y_M = Y,
#            link = "cauchit")
#
# # Student_1
# FisherScoring
# dist3 <- new(FisherScoring)
# dist3$GLMm(X_M = X,
#            Y_M = Y,
#            link = "student")
# # anadir otro parametro para los grados de libertad
#
# # Gumbel
# FisherScoring
# dist3 <- new(FisherScoring)
# dist3$GLMm(X_M = X,
#            Y_M = Y,
#            link = "gumbel")
#
# # Gompertz
# FisherScoring
# dist3 <- new(FisherScoring)
# dist3$GLMm(X_M = X,
#            Y_M = Y,
#            link = "gompertz")
#
#
# # Multicategorical response -----------------------------------------------
#
# # DATA
# {
#   data(addiction)
#   summary(addiction)
#   head(addiction)
#   # vignette("multinomial-addiction1")
#   attach(addiction)
#   ill <- as.factor(ill)
#   addiction$ill<-as.factor(addiction$ill)
#   data1 <- addiction[,c("ill","gender","university")]
#   data1 <- data1[order(data1[,c("ill")]),]
#   data1 <- data1[order(data1[,c("gender")]),]
#   data1 <- data1[order(data1[,c("university")]),]
#   data2 <- na.omit(data1)
#
#   X_1 <- data2[,c("gender","university")]
#   X <- cbind("inter" = rep(1, nrow(X_1)), X_1)
#
#   K = length(unique(data2$ill))
#   Q = K-1
#   P = ncol(X_1)
#   N = nrow(data2)
#
#   Y_EXT_1 <- to.dummy(data2$ill, "category")
#   Y_EXT <- Y_EXT_1[,-ncol(Y_EXT_1)]
#   Y_EXT <- as.vector(t((Y_EXT)))
#   X_EXT <- kronecker(as.matrix(X), diag(Q))
# }
#
# ReferenceF
# dist1 <- new(ReferenceF)
# dist1$GLMref(as.matrix(X_EXT), as.matrix(Y_EXT), K, P, N, Q, link = "logistic")
# dist1$GLMref(as.matrix(X_EXT), as.matrix(Y_EXT), K, P, N, Q, link = "probit")
#
# CumulativeR
# dist2 <- new(CumulativeR)
# dist2$GLMcum(as.matrix(X_EXT), as.matrix(Y_EXT), K, P, N, Q, link = "logistic")
# dist2$GLMcum(as.matrix(X_EXT), as.matrix(Y_EXT), K, P, N, Q, link = "probit")
#
# dist1 <- new(AdjacentR)
# dist1$GLMadj(as.matrix(X_EXT), as.matrix(Y_EXT), K, P, N, Q, link = "logistic")
# dist1$GLMadj(as.matrix(X_EXT), as.matrix(Y_EXT), K, P, N, Q, link = "probit")
#
# dist3 <- new(SequentialR)
# dist3$GLMseq(as.matrix(X_EXT), as.matrix(Y_EXT), K, P, N, Q, link = "logistic")
# dist3$GLMseq(as.matrix(X_EXT), as.matrix(Y_EXT), K, P, N, Q, link = "probit")
#
# # Pruebas con la librerìa usual
#
# data2$ill = relevel(data2$ill, ref = "2")
# (multinom0 <- multinom(ill ~ gender + university, data=data2))
#
# dim(Y_EXT)
# head(X)
#
# (vglm(Y_EXT_1 ~ . , data = as.data.frame(X_1), multinomial(refLevel = 3)))
#
# (vglm(Y_EXT_1 ~ . , data = as.data.frame(X_1), acat(link = "logitlink"),
#       trace = TRUE, maxit = 1))
# (vglm(Y_EXT_1 ~ . , data = as.data.frame(X_1), acat(link = "loglink"),
#       trace = TRUE, maxit = 1))
#
# (vglm(Y_EXT_1 ~ . , data = as.data.frame(X_1), cumulative(link = "logit")))
# (vglm(Y_EXT_1 ~ . , data = as.data.frame(X_1), cumulative(link = "probit")))
#
# (vglm(Y_EXT_1 ~ . , data = as.data.frame(X_1), sratio(link = "logit")))
# (vglm(Y_EXT_1 ~ . , data = as.data.frame(X_1), sratio(link = "probit")))
#
#
# # Others ------------------------------------------------------------------
# vignette("ordinal-knee1")
#
# data(ModeChoice, package="Ecdat")
# # Aca estan cada 4 la eleccion, asì lo necesito en mi modelo
# head(ModeChoice); summary(ModeChoice)
#
# travelmode <- matrix(ModeChoice$mode, byrow = T, ncol = 4)
# colnames(travelmode) <- c("air","train","bus","car")
# travelhinc <- matrix(ModeChoice$hinc, byrow = T, ncol = 4)
# travelhinc <- travelhinc[,1]
# travelinvt <- matrix(ModeChoice$invt, byrow = T, ncol = 4)
# colnames(travelinvt) <- c("invtair","invttrain","invtbus","invtcar")
# travelgc <- matrix(ModeChoice$gc, byrow = T, ncol = 4)
# colnames(travelgc) <- c("gcair","gctrain","gcbus","gccar")
# travelinvt <- sweep(travelinvt[,-1], 1, travelinvt[,1])
# travelgc <- sweep(travelgc[,-1], 1, travelgc[,1])
# Invt <- travelinvt[,1]
# Gc <- travelgc[,1]
# traveldat <- cbind(travelhinc, travelinvt, Invt, travelgc, Gc)
# traveldat <- as.data.frame(traveldat)
#
# head(traveldat)
# head(travelmode)
#
#
