# # Libraries ---------------------------------------------------------------
# library(devtools)
# library(tidyverse)
# library(fs)
# library(Rcpp)
# library(RcppArmadillo)
# library(RcppEigen)
# library(dplyr)
# library(tidyr)
# library(nnet)
# library(varhandle)
# library(catdata)
# library(magic)
# library(VGAM)
# library(dobson)
# library(gtools) # For permutations
#
# library(gridExtra) #For plots
# library(grid)
#
# # Initial configuration ---------------------------------------------------
#
# # load_all()
# # check()
# # use_mit_license("Lorena LEON")
# # check()
# # document()
#
# # usethis::use_rcpp()
# # NO funciona timestwo y he creado el paquete siguiendo los pasos del libro
# # timesTwo(5)
# # Cuando hago esto por segunda vez, funciona.
# # usethis::use_rcpp()
#
# # use_package("forcats")
# # use_package("RcppArmadillo",type = "LinkingTo")
# # use_package("RcppEigen",type = "LinkingTo")
# # Problema de rcpparmadillo not include solo se cambia el orden como sigue:
# # LinkingTo:
# # Rcpp,
# # RcppArmadillo,
# # RcppEigen
# # use_package("RcppEigen",type = "Imports")
# # use_package("RcppArmadillo",type = "Imports")
# # error "The file 'Rcpp.h' should not be included. Please correct to include only 'RcppArmadillo.h'."
# # Borre del distribution.h el include rcpp
# getLoadedDLLs()
# # Using modules
# # We need to tell R that we want to use a C++11 compiler
# Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
# pkgbuild::compile_dll()
# pkgbuild::build()
#
# Rcpp::compileAttributes("~/Desktop/pack_270320")
#
# GLMseq()
# adj_fun()
#
# compileAttributes()
# devtools::load_all()
# devtools::document()
#
# tools::package_native_routine_registration_skeleton(".", character_only = FALSE)
#
# # Binomial case -----------------------------------------------------------
#
# # # DATA
# {
#   beetle
#   i <- 1
#   beetle_ext <- data.frame(x = as.matrix(rep(beetle[i, 1], beetle[i, 2])), y = c(rep(1, beetle[i, 3]), rep(0, beetle[i, 2] - beetle[i, 3])))
#   for (i in 2:nrow(beetle)) {
#     beetle_ext <- rbind(beetle_ext, data.frame(x = as.matrix(rep(beetle[i, 1], beetle[i, 2])), y = c(rep(1, beetle[i, 3]), rep(0, beetle[i, 2] - beetle[i, 3]))))
#   }
#   beetle_ext <- as.data.frame(beetle_ext)
#   names(beetle_ext) <- c("x", "y")
#   head(beetle_ext)
#
#   # Matrix and vectors
#   X <- as.matrix(data.frame(x0 = as.vector(rep(1, nrow(beetle_ext))), x1 = as.vector(unlist(as.vector(beetle_ext$x)))))
#   Y <- beetle_ext$y
#   Y <- matrix(Y)
#   beta <- as.matrix(rep(0, 2))
#   mu <- as.matrix(X) %*% beta
#
#   df_beetle <- data.frame(t(data.frame(matrix(unlist(beetle_ext), nrow = length(beetle_ext), byrow = T))))
#   colnames(df_beetle) <- c("x_beetle", "y_beetle")
# }
#
# # Pruebas con librerìa usual
#
# {
#   y <- c(6, 13, 18, 28, 52, 53, 61, 60)
#   n <- c(59, 60, 62, 56, 63, 59, 62, 60)
#   x <- c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
#   n_y <- n - y
#   beetle.mat <- cbind(y, n_y)
#   glm(beetle.mat ~ x, family = binomial(link = "logit"))
#   glm(beetle.mat ~ x, family = binomial(link = "normal"))
#   glm(beetle.mat ~ x, family = binomial(link = "cauchit"))
#   glm(beetle.mat ~ x, family = binomial(link = "cloglog"))
# }
#
# # Logit
# FisherScoring
# dist3 <- new(FisherScoring)
# dist3$GLMm(
#   X_M = X,
#   Y_M = Y,
#   link = "logistic"
# )
#
#
# df_beetle$y_beetle <- as.factor(df_beetle$y_beetle)
#
# library(plyr)
# df_beetle$y_beetle <- revalue(df_beetle$y_beetle, c("1" = "one", "0" = "zero"))
#
# str((df_beetle))
# summary(df_beetle)
#
# dim(df_beetle)
#
# dist1 <- new(ReferenceF)
# dist1$GLMref(
#   response = "y_beetle",
#   explanatory_complete = c("intercept", "x_beetle"),
#   explanatory_proportional = NA,
#   distribution = "logistic",
#   categories_order = c("zero", "one"),
#   dataframe = df_beetle
# )
#
# # For same results than previous one
# dist1 <- new(ReferenceF)
# dist1$GLMref(
#   response = "y_beetle",
#   explanatory_complete = NA,
#   explanatory_proportional = c("intercept", "x_beetle"),
#   distribution = "logistic",
#   categories_order = c(0, 1),
#   dataframe = df_beetle
# )
#
# # normal
# dist3 <- new(FisherScoring)
# dist3$GLMm(
#   X_M = as.matrix(X),
#   Y_M = as.matrix(Y),
#   link = "normal"
# )
#
#
#
# # Multicategorical response -----------------------------------------------
#
# # ADICTION DATA -----------------------------------------------------------
#
# # DATA
# {
#   data(addiction)
#   summary(addiction)
#   head(addiction)
#   # vignette("multinomial-addiction1")
#   addiction$ill <- as.factor(addiction$ill)
#   data1 <- addiction[, c("ill", "gender", "university", "age")]
#   data2 <- na.omit(data1)
#   # library(plyr)
#   # data2$ill <- revalue(data2$ill, c("0"="d", "1"="o", "2"="z"))
# }
# summary(data2)
# colnames(data2)
# str(data2) # RESPONSE FACTOR. COV AS INT
#
# dist1 <- new(ReferenceF)
# (mod1 <- dist1$GLMref(
#   response = "ill",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("NA"),
#   distribution = "logistic",
#   categories_order = c("0", "1", "2"),
#   dataframe = data2
# ))
#
# mod1$deviance
# mod1$predicted
# head(mod1$fitted)
# head(mod1$residuals)
# head(mod1$dev_r)
# -2 * sum(mod1$dev_log)
#
#
#
# summary.fastLm <- function(object, ...) {
#   coef <- object$coefficients
#   se <- object$stderr
#   tval <- coef / se
#
#   object$coefficients <- cbind(
#     Estimate = coef,
#     "Std. Error" = se,
#     "t value" = tval,
#     "Pr(>|t|)" = 2 * pt(-abs(tval), df = object$df)
#   )
#
#   # cf src/stats/R/lm.R and case with no weights and an intercept
#   f <- object$fitted.values
#   r <- object$residuals
#   # mss <- sum((f - mean(f))^2)
#   mss <- if (object$intercept) sum((f - mean(f))^2) else sum(f^2)
#   rss <- sum(r^2)
#
#   object$r.squared <- mss / (mss + rss)
#   df.int <- if (object$intercept) 1L else 0L
#   n <- length(f)
#   rdf <- object$df
#   object$adj.r.squared <- 1 - (1 - object$r.squared) * ((n - df.int) / rdf)
#   class(object) <- "summary.fastLm"
#   object
# }
#
# summary.fastLm(mod1)
#
# (mod13b <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender", "university"), explanatory_proportional = c("intercept", "age"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# ))
#
# dist1 <- new(AdjacentR)
# (mod1 <- dist1$GLMadj(
#   response = "ill",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("NA"),
#   distribution = "gompetz",
#   categories_order = c("0", "1", "2"),
#   dataframe = dat
# ))
#
# head(data2)
#
# ratio_cum <- new(CumulativeR)
# (mod1 <- GLMcum(
#   response = "ill",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "age"),
#   distribution = "logistic",
#   categories_order = c("0", "1", "2"),
#   dataframe = dat
# ))
#
# ratio_seq <- new(SequentialR)
# (mod1 <- ratio_seq$GLMseq(
#   response = "ill",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "age"),
#   distribution = "gompertz",
#   categories_order = c("0", "1", "2"),
#   dataframe = dat
# ))
#
# (mod2 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c(1, 2, 0), dataframe = dat
# )) # -724.8664
#
# (mod3 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) # -724.8664
# (mod4 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# )) # -737.5639
# (mod5 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#   distribution = "logistic", categories_order = c("0", "2", "1"), dataframe = dat
# )) # -726.4392
# (mod6 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) # -746.9176
#
# (mod7 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# )) #  -730.6547
# (mod8 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "2", "1"), dataframe = dat
# )) #  -730.6547
# (mod8 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("2", "0", "1"), dataframe = dat
# )) #
# (mod9 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) #
#
# summary(data2)
#
# (mod10 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# ))
#
# (mod11 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept", "gender", "university"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# ))
#
# (mod12 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept", "gender", "university"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# )) # -698.9851
# # NO LO DEJA EL OTRO
# (mod13 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept", "gender", "university"), explanatory_proportional = c("age"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) # -698.9851
# # FUNCIONA
# (mod13a <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept", "age", "gender", "university"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) # -679.9081
# (mod13a_1 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept", "age", "gender", "university"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "2", "1"), dataframe = dat
# )) # -679.9081
# (mod13b <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender", "university"), explanatory_proportional = c("intercept", "age"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) # -679.9081
#
# # With categorical variables: Not working with just 2 categories
# data2$gender <- as.factor(data2$gender)
# data2$university <- as.factor(data2$university)
# str(data2)
#
# (mod7 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# ))
# (mod8 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "2", "1"), dataframe = dat
# ))
# (mod9 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# ))
# (mod11 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept", "gender", "university"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# ))
# (mod12 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept", "gender", "university"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# ))
#
#
# # POLITICAL VIEW DATA -----------------------------------------------------
#
# # DATA 2
# Polviews2 <- read.table("http://www.stat.ufl.edu/~aa/cat/data/Polviews2.dat", header = TRUE)
#
# str(Polviews2)
# M2 <- as.data.frame(sapply(Polviews2[, c("ideology", "party", "gender")], unclass))
# M2 <- as.data.frame(M2)
# sum(complete.cases(M2))
# M2$ideology <- as.factor(M2$ideology)
# str(M2)
# summary(M2)
#
# dist1 <- new(ReferenceF)
# (mod1_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -980.4022
# (mod2_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "5", "4"), dataframe = M2
# )) # -980.4022
# (mod3_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2
# )) # -980.4022
# (mod4_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -1043.362
# (mod5_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -977.8485
# (mod6_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "gender", "party"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -774.6079
# (mod7_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("party"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -977.8485
# (mod8_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("intercept", "party"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -996.6733
#
# # NO FUNCIONA ALLA
# (mod9_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("party"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -977.8485
#
# (mod10_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2
# )) # -977.8485
#
# # ALGO RARO TIENE PARTY
# (mod11_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2
# )) # -775.7563
#
#
# (mod6_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#   distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2
# ))
#
# # SI FUNCIONA
# (mod14_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party", "gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -774.6079
# (mod15_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -775.5648
# (mod16_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("5", "4", "3", "2", "1"), dataframe = M2
# )) # -775.5648
#
# # TIENE OTRO ORDEN
# (mod16_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2
# )) # -775.5648
# (mod16_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("5", "2", "3", "4", "1"), dataframe = M2
# )) # -775.5648
#
# # ALGO RARO CON PARTY
# (mod16_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("5", "2", "3", "4", "1"), dataframe = M2
# )) # -775.5648
#
# # GENDER WITHOUT INTERCEPT
#
# # All posible permutations
# a1 <- dist1$GLMref(
#   response = "ideology",
#   explanatory_complete = c("intercept", "party"),
#   explanatory_proportional = "gender",
#   distribution = "logistic",
#   categories_order = c(0, 1, 3, 2, 4),
#   dataframe = M2
# )
#
# ## Reference
# ### Invariance under permutations
# # Multinomial logit model is invariant under all permutations of the response categories
#
# dist1 <- new(ReferenceF)
# all_permutations <- permutations(v = c(0, 1, 2, 3, 4), repeats.allowed = F, n = 5, r = 5)
# Log_lik_Vec <- NA
# for (element in 1:nrow(all_permutations)) {
#   a1 <- dist1$GLMref(
#     response = "ideology",
#     explanatory_complete = c("intercept", "party", "gender"),
#     explanatory_proportional = NA,
#     distribution = "logistic",
#     categories_order = all_permutations[element, ],
#     dataframe = M2
#   )
#   Log_lik_Vec[element] <- a1$`Log-likelihood`
# }
# Log_lik_Vec
#
#
# # MULTINOMIAL PARTY EX 8.3 TUTZ -------------------------------------------
# # DATA
# {
#   partypref <- matrix(data = c(
#     114, 10, 53, 224, 134, 9, 42, 226, 114, 8, 23, 174, 339, 30, 13,
#     414, 42, 5, 44, 161, 88, 10, 60, 171, 90, 8, 31, 168, 413, 23, 14, 375
#   ), nrow = 8, byrow = TRUE)
#   partydat <- data.frame(
#     party = c(
#       rep("CDU", sum(partypref[, 1])), rep("SPD", sum(partypref[, 4])),
#       rep("The Liberals", sum(partypref[, 2])), rep("The Greens", sum(partypref[, 3]))
#     ),
#     sex = c(
#       rep(0, sum(partypref[1:4, 1])), rep(1, sum(partypref[5:8, 1])),
#       rep(0, sum(partypref[1:4, 4])), rep(1, sum(partypref[5:8, 4])),
#       rep(0, sum(partypref[1:4, 2])), rep(1, sum(partypref[5:8, 2])),
#       rep(0, sum(partypref[1:4, 3])), rep(1, sum(partypref[5:8, 3]))
#     ),
#     age = c(
#       rep(c(1:4, 1:4), partypref[, 1]), rep(c(1:4, 1:4), partypref[, 4]),
#       rep(c(1:4, 1:4), partypref[, 2]), rep(c(1:4, 1:4), partypref[, 3])
#     )
#   )
#   partydat$age <- as.factor(partydat$age)
# }
# head(partydat)
# str(partydat)
# summary(partydat)
#
# (m_party_1 <- dist1$GLMref(
#   response = "party", explanatory_complete = c("intercept", "sex", "age"), explanatory_proportional = NA,
#   distribution = "logistic", categories_order = c("SPD", "The Greens", "The Liberals", "CDU"), dataframe = partydat
# )) # -3521.484
#
# # MULTINOMIAL TRAVEL EX 8.3 TUTZ ------------------------------------------
# # vignette("multinomial-travel")
#
# library(mlogit)
# data(ModeChoice, package = "Ecdat")
# head(ModeChoice)
# travel.long <- mlogit.data(ModeChoice, choice = "mode", shape = "long", alt.levels = c("air", "train", "bus", "car"))
#
# head(travel.long)
# travel.kat.id <- mlogit(mode ~ invt + gc | hinc, data = travel.long)
# logLik(travel.kat.id)
#
# # Now with VGAM
# travelmode <- matrix(ModeChoice$mode, byrow = T, ncol = 4)
# colnames(travelmode) <- c("air", "train", "bus", "car")
# travelhinc <- matrix(ModeChoice$hinc, byrow = T, ncol = 4)
# travelhinc <- travelhinc[, 1]
# travelinvt <- matrix(ModeChoice$invt, byrow = T, ncol = 4)
# colnames(travelinvt) <- c("invtair", "invttrain", "invtbus", "invtcar")
# travelgc <- matrix(ModeChoice$gc, byrow = T, ncol = 4)
# colnames(travelgc) <- c("gcair", "gctrain", "gcbus", "gccar")
# travelinvt <- sweep(travelinvt[, -1], 1, travelinvt[, 1])
# travelgc <- sweep(travelgc[, -1], 1, travelgc[, 1])
# Invt <- travelinvt[, 1]
# Gc <- travelgc[, 1]
# traveldat <- cbind(travelhinc, travelinvt, Invt, travelgc, Gc)
# traveldat <- as.data.frame(traveldat)
#
# head(traveldat)
# (fit <- vglm(travelmode ~ Invt + Gc + travelhinc, multinomial(parallel = FALSE ~ travelhinc, refLevel = 1),
#              xij = list(Invt ~ invttrain + invtbus + invtcar, Gc ~ gctrain + gcbus + gccar),
#              form2 = ~ Invt + invttrain + invtbus + invtcar + Gc + gctrain + gcbus + gccar + travelhinc,
#              data = traveldat, trace = TRUE
# ))
#
#
# # ECONOMETRIC TRAVEL CHOICE -----------------------------------------------
#
# # MY DATA
# {
#   library(mlogit)
#   data(ModeChoice, package = "Ecdat")
#   head(ModeChoice)
#   travel.long <- mlogit.data(ModeChoice, choice = "mode", shape = "long", alt.levels = c("air", "train", "bus", "car"))
#   head(travel.long)
#   choice <- sub(".*\\.", "", rownames(travel.long))
#   indv <- sub("\\..*", "", rownames(travel.long))
#   travel.long88 <- as.data.frame(cbind(indv, choice, travel.long))
# }
# head(travel.long88, 5)
# str(travel.long88)
# travel.long88$choice
# # library(plyr)
# # travel.long88$choice <- revalue(travel.long88$choice, c("air"="a", "train"="t", "bus"="b", "car"="c"))
# # dist3 <- new(ReferenceF)
# # (exp_8_3 <- dist3$GLMref_ec(
# #   response = "choice", actual_response = "mode",
# #   individuals = "indv",
# #   explanatory_complete = c("intercept", "hinc"),
# #   depend_y = c("gc", "invt"),
# #   distribution = "logistic", categories_order = c("t", "b", "c", "a"), dataframe = travel.long88,
# #   design = "tutz"
# # ))
#
# dist3 <- new(ReferenceF)
# (exp_8_3 <- dist3$GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc"),
#   depend_y = c("gc", "invt"),
#   distribution = "logistic", categories_order = c("train", "bus", "car", "air"), dataframe = travel.long88,
#   design = "tutz"
# ))
# exp_8_3$`Log-likelihood`
# exp_8_3$Coefficients
#
#
# # Robustness of Student link function in multinomial choice models
# # The log-likelihood obtained with the MNL is −185.91 as obtained by Louviere et al. (2000) page 157.
# (table3 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "logistic", categories_order = c("air", "train", "bus", "car"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 2.0
# ))
#
# # The log-likelihoods obtained with the (reference, F ν ∗ , Z) j 0
# # models were −185.65, −183.79, −142, −183.49 respectively with
# # the four reference alternatives j 0 =air, j 0 =bus, j 0 =car, j 0 =train and
# # correspondind degree of freedom ν ∗ = 3,
# # ν ∗ = 30, ν ∗ = 0.2, ν ∗ = 1.35.
#
# # j_0=air, v* = 3, ll = −185.65
# # DOES NOT WORK
# (table4 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("train", "car", "bus", "air"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 3.0
# ))
#
# (table4 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("train", "bus", "car", "air"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 3
# ))
#
# (table4 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("bus", "train", "car", "air"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 3.0
# ))
#
# (table4 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("bus", "car", "train", "air"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 3.0
# ))
#
# (table4 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("car", "train", "bus", "air"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 3.0
# ))
#
# (table4 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("car", "bus", "train", "air"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 3.0
# ))
#
# # j_0=bus, v* = 30, ll = −183.79
# (table4 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("air", "train", "car", "bus"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 30
# ))
#
# # j_0=car, v* = 0.2, ll = −142
# ## DOES NOT WORK
# (train_1.35 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("air", "train", "bus", "car"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 0.2
# ))
# # j_0=car, v* = 1.0, ll = −142
# (train_1.35 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("air", "train", "bus", "car"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 1
# ))
#
# # j_0=train, v* = 1.35, ll = −183.49
# (train_1.35 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("air", "bus", "car", "train"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 1.35
# ))
#
# (table4 <- dist3$GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept"),
#   depend_y = c("ttme"),
#   distribution = "logistic", categories_order = c("air", "train", "bus", "car"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 2.0
# ))
#
#
#
#
#
# # ANOTHER EXAMPLE FOR ECONOMETRIC MODEL -----------------------------------
#
#
# data("Heating", package = "Ecdat")
#
# # Heating is a "horizontal" data.frame with three choice-specific
# # variables (ic: investment cost, oc: operating cost) and some
# # individual-specific variables (income, region, rooms)
#
# H <- mlogit.data(Heating, shape = "wide", choice = "depvar", varying = c(3:12))
# head(H)
#
# (mi2 <- mlogit(depvar ~ oc + ic | income, H, reflevel = "hp"))
# logLik(mi2)
#
# choice <- sub(".*\\.", "", rownames(H))
# indv <- sub("\\..*", "", rownames(H))
#
# dat4 <- as.data.frame(cbind(indv, choice, H))
# head(dat4, 8)
# str(dat4)
# dist3 <- new(ReferenceF)
# (A98 <- dist3$GLMref_ec(
#   response = "choice", actual_response = "depvar",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "income"),
#   depend_y = c("oc", "ic"),
#   distribution = "logistic",
#   categories_order = c("ec", "er", "gc", "gr", "hp"),
#   dataframe = dat,
#   design = "tutz"
# ))
#
#
# # JSS ---------------------------------------------------------------------
#
# {
#   library(mlogit)
#   data(ModeChoice, package = "Ecdat")
#   head(ModeChoice)
#   travel.long <- mlogit.data(ModeChoice, choice = "mode", shape = "long", alt.levels = c("air", "train", "bus", "car"))
#   head(travel.long)
#   choice <- sub(".*\\.", "", rownames(travel.long))
#   indv <- sub("\\..*", "", rownames(travel.long))
#   travel.long88 <- as.data.frame(cbind(indv, choice, travel.long))
# }
# head(travel.long88, 5)
# str(travel.long88)
#
# travel_dat1 <- travel.long88[travel.long88$mode == T, ]
# head(travel_dat1)
#
#
# # e : REFERENCE, LOGISTIC, COMPLETE
# (l_1 <- GLMref(
#   response = "choice",
#   explanatory_complete = c("intercept", "psize"),
#   explanatory_proportional = c("hinc"),
#   distribution = "logistic",
#   categories_order = c("air", "train", "bus", "car"),
#   dataframe = travel_dat1
# ))
#
# (l_1 <- GLMcum(
#   response = "choice",
#   explanatory_complete = c("intercept", "psize"),
#   explanatory_proportional = c("hinc"),
#   distribution = "normal",
#   categories_order = c("air", "train", "bus", "car"),
#   dataframe = travel_dat1
# ))
#
# (l_1 <- GLMcum(
#   response = "choice",
#   explanatory_complete = c("intercept", "psize"),
#   explanatory_proportional = c("hinc"),
#   distribution = "cauchit",
#   categories_order = c("air", "train", "bus", "car"),
#   dataframe = travel_dat1
# ))
#
# (l_1 <- GLMcum(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc"),
#   distribution = "gompertz",
#   categories_order = c("air", "train", "bus", "car"),
#   dataframe = travel_dat1
# ))
#
#
# # e : REFERENCE, LOGISTIC, PROPORTIONAL
# (l_2 <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "logistic",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# # e : REFERENCE, CAUCHIT, COMPLETE
# (l_3 <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "cauchit",
#   categories_order = c("air", "train", "bus", "car"),
#   dataframe = travel_dat1
# ))
#
# # Then we change the reference category and estimate again the three reference models:
#
# # e : REFERENCE, LOGISTIC, COMPLETE
# (l_1prime <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "logistic",
#   categories_order = c("air", "train", "car", "bus"),
#   dataframe = travel_dat1
# ))
#
#
# # e : REFERENCE, LOGISTIC, PROPORTIONAL
# (l_2prime <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "logistic",
#   categories_order = c("air", "train", "car", "bus"),
#   dataframe = travel_dat1
# ))
#
# # e : REFERENCE, CAUCHIT, COMPLETE
# (l_3prime <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "cauchit",
#   categories_order = c("air", "train", "car", "bus"),
#   dataframe = travel_dat1
# ))
#
#
# # The log-likelihoods l_1 and l_1prime are equal (since the canonical model is invariant under
# # all permutations) whereas the log-likelihoods l_2 and l_2prime are different (respectively
# # l_3 and l_3prime).
#
#
# # 4. Adjacent models for ordinal response
#
# # Equivalence between (adjacent, logistic, complete) and (reference, logistic, complete) models.
#
# # e : REFERENCE, LOGISTIC, COMPLETE
# dist1 <- new(ReferenceF)
# (estimation <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("hinc"),
#   distribution = "logistic",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# # e : ADJACENT, LOGISTIC, COMPLETE
# dist2 <- new(AdjacentR)
# (estimation_prime <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "logistic",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
#
#
# # Try equality with defined matrix A
#
# # Invariance under permutations
#
# # ADJACENT, CAUCHY, COMPLETE
#
# (estimation_1 <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "cauchit",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# # ADJACENT, GOMPERTZ, PROPORTIONAL
# (estimation_2 <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "gompertz",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# # ADJACENT, CAUCHY, COMPLETE (Reverse O)
#
# (estimation_1r <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "cauchit",
#   categories_order = c("train", "car", "bus", "air"),
#   dataframe = travel_dat1
# ))
#
# # ADJACENT, GOMPERTZ, PROPORTIONAL (Reverse O)
# (estimation_2r <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "gompertz",
#   categories_order = c("train", "car", "bus", "air"),
#   dataframe = travel_dat1
# ))
#
# # ADJACENT, GUMBEL, PROPORTIONAL (Reverse O)
# (estimation_3r <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "gumbel",
#   categories_order = c("train", "car", "bus", "air"),
#   dataframe = travel_dat1
# ))
#
# # Remark that the log-likelihoods l_1 and l_1r are equal since the Cauchy distribution is sym-
# # metric whereas the log-likelihoods l_2 and l_2r are different since the Gompertz distribution
# # is not symmetric. Moreover if the Gumbel distribution is used we the reverse order then
# # the log-likelihoods l_2 and l_3r are equal since the Gumbel distribution is the symmetric
# # of the Gompertz distribution.
#
# # 5. Cumulative models for ordinal response
#
# # The equivalence between the (cumulative, Gompertz, proportional)
# # and (sequential, Gompertz, proportional)
# # models has been demonstrated by Läärä and Matthews (1985).
#
# # Cumulative, Gompertz, Proportional
# (estimation <- GLMcum(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "gompertz",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# dat_lok_1 <- data.frame("LogLikIter" = estimation$LogLikIter, "iteration" = 1:nrow(dat_lok_1))
#
# library(ggplot2)
# ggplot(dat_lok_1[-1, ], aes(x = iteration, y = LogLikIter)) +
#   geom_point()
#
#
# # Sequential, Gompertz, Proportional
# dist2 <- new(SequentialR)
# (estimation_prime <- dist2$GLMseq(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "gompertz",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# # KNEE DATASET ------------------------------------------------------------
# # DEMOSTRAR EQUIVALENCIA CON OTRO DATASET
#
# # CUMULATIVE SI PARECE ESTAR BIEN DEFINIDO,
# # PERO QUIZAS LA FORMA EN QUE SE INICIALIZA NO ESTA BIEN
# library(catdata)
# data(knee)
# knee$R4 <- as.factor(knee$R4)
# dist1 <- new(CumulativeR)
# (estimation <- dist1$GLMcum(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "logit",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee
# ))
#
#
#
# # NO ES ESTIMABLE CON EL OTRO MODELO
# (estimation <- dist1$GLMcum(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee
# ))
#
# (estimation_prime <- dist2$GLMseq(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee
# ))
#
# (estimation <- dist1$GLMcum(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "logistic",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee
# ))
#
# (estimation_prime <- dist2$GLMseq(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("S·", "Sex", "Age"),
#   distribution = "gumbel",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee
# ))
#
# dat <- read.csv2("~/Desktop/Test package/dat_cum.csv")
# summary(dat)
# str(dat)
# head(dat)
#
# # FUNCIONA BIEN
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("NA"),
#   distribution = "logistic",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("NA"),
#   distribution = "gompertz",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # MOD3
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("public"),
#   distribution = "gompertz",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # MOD4
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("public"),
#   distribution = "normal",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # MOD5
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("pared", "gpa"),
#   explanatory_proportional = c("intercept", "public"),
#   distribution = "normal",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# # MOD6
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("public", "gpa"),
#   distribution = "normal",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat, beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # MOD7
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("public"),
#   distribution = "cauchit",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat, beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# # MOD 7
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("public"),
#   distribution = "cauchit",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("NA"),
#   distribution = "cauchit",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat, beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
#
#
# # ORDINAL LIBRARY ---------------------------------------------------------
# # DATA
# {
#   library("ordinal")
#   head(wine)
#   wine$temp <- as.numeric((wine$temp))
#   wine$contact <- as.numeric((wine$contact))
#   str(wine)
# }
#
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("NA"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # NO
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept", "temp"),
#   explanatory_proportional = c("NA"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # NO
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept", "temp", "contact"),
#   explanatory_proportional = c("NA"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp"),
#   distribution = "cauchit",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp"),
#   distribution = "logistic",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# # e : ADJACENT, LOGISTIC, COMPLETE
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp", "contact"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("TRUE", "sd"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp", "contact"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine
# ))
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp", "contact"),
#   distribution = "cauchit",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine
# ))
#
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp", "contact"),
#   distribution = "cauchit",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine
# ))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # CREATE VIGNETTE ---------------------------------------------------------
# # library(rmarkdown); library(devtools)
# #
# # usethis::use_vignette("my-vignette")
#
# ratio_ref <- new(ReferenceF)
# ratio_adj <- new(AdjacentR)
# ratio_cum <- new(CumulativeR)
# ratio_seq <- new(SequentialR)
#
# # Libraries ---------------------------------------------------------------
# library(devtools)
# library(tidyverse)
# library(fs)
# library(Rcpp)
# library(RcppArmadillo)
# library(RcppEigen)
# library(dplyr)
# library(tidyr)
# library(nnet)
# library(varhandle)
# library(catdata)
# library(magic)
# library(VGAM)
# library(dobson)
# library(gtools) # For permutations
#
# # Initial configuration ---------------------------------------------------
#
# # load_all()
# # check()
# # use_mit_license("Lorena LEON")
# # check()
# # document()
#
# # usethis::use_rcpp()
# # NO funciona timestwo y he creado el paquete siguiendo los pasos del libro
# # timesTwo(5)
# # Cuando hago esto por segunda vez, funciona.
# # usethis::use_rcpp()
#
# # use_package("forcats")
# # use_package("RcppArmadillo",type = "LinkingTo")
# # use_package("RcppEigen",type = "LinkingTo")
# # Problema de rcpparmadillo not include solo se cambia el orden como sigue:
# # LinkingTo:
# # Rcpp,
# # RcppArmadillo,
# # RcppEigen
# # use_package("RcppEigen",type = "Imports")
# # use_package("RcppArmadillo",type = "Imports")
# # error "The file 'Rcpp.h' should not be included. Please correct to include only 'RcppArmadillo.h'."
# # Borre del distribution.h el include rcpp
# getLoadedDLLs()
# # Using modules
# # We need to tell R that we want to use a C++11 compiler
# Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
# pkgbuild::compile_dll()
# compileAttributes()
# load_all()
#
# # Binomial case -----------------------------------------------------------
#
# # # DATA
# {
#   beetle
#   i <- 1
#   beetle_ext <- data.frame(x = as.matrix(rep(beetle[i, 1], beetle[i, 2])), y = c(rep(1, beetle[i, 3]), rep(0, beetle[i, 2] - beetle[i, 3])))
#   for (i in 2:nrow(beetle)) {
#     beetle_ext <- rbind(beetle_ext, data.frame(x = as.matrix(rep(beetle[i, 1], beetle[i, 2])), y = c(rep(1, beetle[i, 3]), rep(0, beetle[i, 2] - beetle[i, 3]))))
#   }
#   beetle_ext <- as.data.frame(beetle_ext)
#   names(beetle_ext) <- c("x", "y")
#   head(beetle_ext)
#
#   # Matrix and vectors
#   X <- as.matrix(data.frame(x0 = as.vector(rep(1, nrow(beetle_ext))), x1 = as.vector(unlist(as.vector(beetle_ext$x)))))
#   Y <- beetle_ext$y
#   Y <- matrix(Y)
#   beta <- as.matrix(rep(0, 2))
#   mu <- as.matrix(X) %*% beta
#
#   df_beetle <- data.frame(t(data.frame(matrix(unlist(beetle_ext), nrow = length(beetle_ext), byrow = T))))
#   colnames(df_beetle) <- c("x_beetle", "y_beetle")
# }
#
# # Pruebas con librerìa usual
#
# {
#   y <- c(6, 13, 18, 28, 52, 53, 61, 60)
#   n <- c(59, 60, 62, 56, 63, 59, 62, 60)
#   x <- c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
#   n_y <- n - y
#   beetle.mat <- cbind(y, n_y)
#   glm(beetle.mat ~ x, family = binomial(link = "logit"))
#   glm(beetle.mat ~ x, family = binomial(link = "normal"))
#   glm(beetle.mat ~ x, family = binomial(link = "cauchit"))
#   glm(beetle.mat ~ x, family = binomial(link = "cloglog"))
# }
#
# # Logit
# FisherScoring
# dist3 <- new(FisherScoring)
# dist3$GLMm(
#   X_M = X,
#   Y_M = Y,
#   link = "logistic"
# )
#
#
# df_beetle$y_beetle <- as.factor(df_beetle$y_beetle)
#
# library(plyr)
# df_beetle$y_beetle <- revalue(df_beetle$y_beetle, c("1" = "one", "0" = "zero"))
#
# str((df_beetle))
# summary(df_beetle)
#
# dim(df_beetle)
#
# dist1 <- new(ReferenceF)
# dist1$GLMref(
#   response = "y_beetle",
#   explanatory_complete = c("intercept", "x_beetle"),
#   explanatory_proportional = NA,
#   distribution = "logistic",
#   categories_order = c("zero", "one"),
#   dataframe = df_beetle
# )
#
# # For same results than previous one
# dist1 <- new(ReferenceF)
# dist1$GLMref(
#   response = "y_beetle",
#   explanatory_complete = NA,
#   explanatory_proportional = c("intercept", "x_beetle"),
#   distribution = "logistic",
#   categories_order = c(0, 1),
#   dataframe = df_beetle
# )
#
# # normal
# dist3 <- new(FisherScoring)
# dist3$GLMm(
#   X_M = as.matrix(X),
#   Y_M = as.matrix(Y),
#   link = "normal"
# )
#
#
#
# # Multicategorical response -----------------------------------------------
#
# # ADICTION DATA -----------------------------------------------------------
#
# # DATA
# {
#   data(addiction)
#   summary(addiction)
#   head(addiction)
#   # vignette("multinomial-addiction1")
#   addiction$ill <- as.factor(addiction$ill)
#   data1 <- addiction[, c("ill", "gender", "university", "age")]
#   data2 <- na.omit(data1)
# }
# summary(data2)
# colnames(data2)
# str(data2) # RESPONSE FACTOR. COV AS INT
# library(plyr)
# # data2$ill <- revalue(data2$ill, c("0"="d", "1"="o", "2"="z"))
#
# dist1 <- new(ReferenceF)
# (mod1 <- dist1$GLMref(
#   response = "ill",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("age"),
#   distribution = "logistic",
#   categories_order = c("0", "1", "2"),
#   dataframe = dat
# ))
#
# (mod13b <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender", "university"), explanatory_proportional = c("intercept", "age"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# ))
#
# dist1 <- new(AdjacentR)
# (mod1 <- dist1$GLMadj(
#   response = "ill",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("NA"),
#   distribution = "gompetz",
#   categories_order = c("0", "1", "2"),
#   dataframe = dat
# ))
#
# head(data2)
#
# ratio_cum <- new(CumulativeR)
# (mod1 <- GLMcum(
#   response = "ill",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "age"),
#   distribution = "logistic",
#   categories_order = c("0", "1", "2"),
#   dataframe = dat
# ))
#
# ratio_seq <- new(SequentialR)
# (mod1 <- ratio_seq$GLMseq(
#   response = "ill",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "age"),
#   distribution = "gompertz",
#   categories_order = c("0", "1", "2"),
#   dataframe = dat
# ))
#
# (mod2 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c(1, 2, 0), dataframe = dat
# )) # -724.8664
#
# (mod3 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) # -724.8664
# (mod4 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# )) # -737.5639
# (mod5 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#   distribution = "logistic", categories_order = c("0", "2", "1"), dataframe = dat
# )) # -726.4392
# (mod6 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) # -746.9176
#
# (mod7 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# )) #  -730.6547
# (mod8 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "2", "1"), dataframe = dat
# )) #  -730.6547
# (mod8 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("2", "0", "1"), dataframe = dat
# )) #
# (mod9 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) #
#
# summary(data2)
#
# (mod10 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# ))
#
# (mod11 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept", "gender", "university"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# ))
#
# (mod12 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept", "gender", "university"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# )) # -698.9851
# # NO LO DEJA EL OTRO
# (mod13 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept", "gender", "university"), explanatory_proportional = c("age"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) # -698.9851
# # FUNCIONA
# (mod13a <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept", "age", "gender", "university"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) # -679.9081
# (mod13a_1 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept", "age", "gender", "university"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "2", "1"), dataframe = dat
# )) # -679.9081
# (mod13b <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender", "university"), explanatory_proportional = c("intercept", "age"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# )) # -679.9081
#
# # With categorical variables: Not working with just 2 categories
# data2$gender <- as.factor(data2$gender)
# data2$university <- as.factor(data2$university)
# str(data2)
#
# (mod7 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# ))
# (mod8 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "2", "1"), dataframe = dat
# ))
# (mod9 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "0"), dataframe = dat
# ))
# (mod11 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept", "gender", "university"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# ))
# (mod12 <- dist1$GLMref(
#   response = "ill", explanatory_complete = c("intercept", "gender", "university"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("0", "1", "2"), dataframe = dat
# ))
#
#
# # POLITICAL VIEW DATA -----------------------------------------------------
#
# # DATA 2
# Polviews2 <- read.table("http://www.stat.ufl.edu/~aa/cat/data/Polviews2.dat", header = TRUE)
#
# str(Polviews2)
# M2 <- as.data.frame(sapply(Polviews2[, c("ideology", "party", "gender")], unclass))
# M2 <- as.data.frame(M2)
# sum(complete.cases(M2))
# M2$ideology <- as.factor(M2$ideology)
# str(M2)
# summary(M2)
#
# dist1 <- new(ReferenceF)
# (mod1_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -980.4022
# (mod2_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "5", "4"), dataframe = M2
# )) # -980.4022
# (mod3_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2
# )) # -980.4022
# (mod4_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -1043.362
# (mod5_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -977.8485
# (mod6_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "gender", "party"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -774.6079
# (mod7_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("party"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -977.8485
# (mod8_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("intercept", "party"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -996.6733
#
# # NO FUNCIONA ALLA
# (mod9_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("party"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -977.8485
#
# (mod10_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2
# )) # -977.8485
#
# # ALGO RARO TIENE PARTY
# (mod11_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2
# )) # -775.7563
#
#
# (mod6_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#   distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2
# ))
#
# # SI FUNCIONA
# (mod14_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party", "gender"), explanatory_proportional = c("NA"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -774.6079
# (mod15_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2
# )) # -775.5648
# (mod16_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("5", "4", "3", "2", "1"), dataframe = M2
# )) # -775.5648
#
# # TIENE OTRO ORDEN
# (mod16_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2
# )) # -775.5648
# (mod16_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("5", "2", "3", "4", "1"), dataframe = M2
# )) # -775.5648
#
# # ALGO RARO CON PARTY
# (mod16_2 <- dist1$GLMref(
#   response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#   distribution = "logistic", categories_order = c("5", "2", "3", "4", "1"), dataframe = M2
# )) # -775.5648
#
# # GENDER WITHOUT INTERCEPT
#
# # All posible permutations
# a1 <- dist1$GLMref(
#   response = "ideology",
#   explanatory_complete = c("intercept", "party"),
#   explanatory_proportional = "gender",
#   distribution = "logistic",
#   categories_order = c(0, 1, 3, 2, 4),
#   dataframe = M2
# )
#
# ## Reference
# ### Invariance under permutations
# # Multinomial logit model is invariant under all permutations of the response categories
#
# dist1 <- new(ReferenceF)
# all_permutations <- permutations(v = c(0, 1, 2, 3, 4), repeats.allowed = F, n = 5, r = 5)
# Log_lik_Vec <- NA
# for (element in 1:nrow(all_permutations)) {
#   a1 <- dist1$GLMref(
#     response = "ideology",
#     explanatory_complete = c("intercept", "party", "gender"),
#     explanatory_proportional = NA,
#     distribution = "logistic",
#     categories_order = all_permutations[element, ],
#     dataframe = M2
#   )
#   Log_lik_Vec[element] <- a1$`Log-likelihood`
# }
# Log_lik_Vec
#
#
# # MULTINOMIAL PARTY EX 8.3 TUTZ -------------------------------------------
# # DATA
# {
#   partypref <- matrix(data = c(
#     114, 10, 53, 224, 134, 9, 42, 226, 114, 8, 23, 174, 339, 30, 13,
#     414, 42, 5, 44, 161, 88, 10, 60, 171, 90, 8, 31, 168, 413, 23, 14, 375
#   ), nrow = 8, byrow = TRUE)
#   partydat <- data.frame(
#     party = c(
#       rep("CDU", sum(partypref[, 1])), rep("SPD", sum(partypref[, 4])),
#       rep("The Liberals", sum(partypref[, 2])), rep("The Greens", sum(partypref[, 3]))
#     ),
#     sex = c(
#       rep(0, sum(partypref[1:4, 1])), rep(1, sum(partypref[5:8, 1])),
#       rep(0, sum(partypref[1:4, 4])), rep(1, sum(partypref[5:8, 4])),
#       rep(0, sum(partypref[1:4, 2])), rep(1, sum(partypref[5:8, 2])),
#       rep(0, sum(partypref[1:4, 3])), rep(1, sum(partypref[5:8, 3]))
#     ),
#     age = c(
#       rep(c(1:4, 1:4), partypref[, 1]), rep(c(1:4, 1:4), partypref[, 4]),
#       rep(c(1:4, 1:4), partypref[, 2]), rep(c(1:4, 1:4), partypref[, 3])
#     )
#   )
#   partydat$age <- as.factor(partydat$age)
# }
# head(partydat)
# str(partydat)
# summary(partydat)
#
# (m_party_1 <- dist1$GLMref(
#   response = "party", explanatory_complete = c("intercept", "sex", "age"), explanatory_proportional = NA,
#   distribution = "logistic", categories_order = c("SPD", "The Greens", "The Liberals", "CDU"), dataframe = partydat
# )) # -3521.484
#
# # MULTINOMIAL TRAVEL EX 8.3 TUTZ ------------------------------------------
# # vignette("multinomial-travel")
#
# library(mlogit)
# data(ModeChoice, package = "Ecdat")
# head(ModeChoice)
# travel.long <- mlogit.data(ModeChoice, choice = "mode", shape = "long", alt.levels = c("air", "train", "bus", "car"))
#
# head(travel.long)
# travel.kat.id <- mlogit(mode ~ invt + gc | hinc, data = travel.long)
# logLik(travel.kat.id)
#
# # Now with VGAM
# travelmode <- matrix(ModeChoice$mode, byrow = T, ncol = 4)
# colnames(travelmode) <- c("air", "train", "bus", "car")
# travelhinc <- matrix(ModeChoice$hinc, byrow = T, ncol = 4)
# travelhinc <- travelhinc[, 1]
# travelinvt <- matrix(ModeChoice$invt, byrow = T, ncol = 4)
# colnames(travelinvt) <- c("invtair", "invttrain", "invtbus", "invtcar")
# travelgc <- matrix(ModeChoice$gc, byrow = T, ncol = 4)
# colnames(travelgc) <- c("gcair", "gctrain", "gcbus", "gccar")
# travelinvt <- sweep(travelinvt[, -1], 1, travelinvt[, 1])
# travelgc <- sweep(travelgc[, -1], 1, travelgc[, 1])
# Invt <- travelinvt[, 1]
# Gc <- travelgc[, 1]
# traveldat <- cbind(travelhinc, travelinvt, Invt, travelgc, Gc)
# traveldat <- as.data.frame(traveldat)
#
# head(traveldat)
# (fit <- vglm(travelmode ~ Invt + Gc + travelhinc, multinomial(parallel = FALSE ~ travelhinc, refLevel = 1),
#              xij = list(Invt ~ invttrain + invtbus + invtcar, Gc ~ gctrain + gcbus + gccar),
#              form2 = ~ Invt + invttrain + invtbus + invtcar + Gc + gctrain + gcbus + gccar + travelhinc,
#              data = traveldat, trace = TRUE
# ))
#
#
# # ECONOMETRIC TRAVEL CHOICE -----------------------------------------------
#
# # MY DATA
# {
#   library(mlogit)
#   data(ModeChoice, package = "Ecdat")
#   head(ModeChoice)
#   travel.long <- mlogit.data(ModeChoice, choice = "mode", shape = "long", alt.levels = c("air", "train", "bus", "car"))
#   head(travel.long)
#   choice <- sub(".*\\.", "", rownames(travel.long))
#   indv <- sub("\\..*", "", rownames(travel.long))
#   travel.long88 <- as.data.frame(cbind(indv, choice, travel.long))
# }
# head(travel.long88, 5)
# str(travel.long88)
# travel.long88$choice
# # library(plyr)
# # travel.long88$choice <- revalue(travel.long88$choice, c("air"="a", "train"="t", "bus"="b", "car"="c"))
# # dist3 <- new(ReferenceF)
# # (exp_8_3 <- dist3$GLMref_ec(
# #   response = "choice", actual_response = "mode",
# #   individuals = "indv",
# #   explanatory_complete = c("intercept", "hinc"),
# #   depend_y = c("gc", "invt"),
# #   distribution = "logistic", categories_order = c("t", "b", "c", "a"), dataframe = travel.long88,
# #   design = "tutz"
# # ))
#
# dist3 <- new(ReferenceF)
# (exp_8_3 <- dist3$GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc"),
#   depend_y = c("gc", "invt"),
#   distribution = "logistic", categories_order = c("train", "bus", "car", "air"), dataframe = travel.long88,
#   design = "tutz"
# ))
# exp_8_3$`Log-likelihood`
# exp_8_3$Coefficients
#
#
# # Robustness of Student link function in multinomial choice models
# # The log-likelihood obtained with the MNL is −185.91 as obtained by Louviere et al. (2000) page 157.
# (table3 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "logistic", categories_order = c("air", "train", "bus", "car"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 2.0
# ))
#
# # The log-likelihoods obtained with the (reference, F ν ∗ , Z) j 0
# # models were −185.65, −183.79, −142, −183.49 respectively with
# # the four reference alternatives j 0 =air, j 0 =bus, j 0 =car, j 0 =train and
# # correspondind degree of freedom ν ∗ = 3,
# # ν ∗ = 30, ν ∗ = 0.2, ν ∗ = 1.35.
#
# # j_0=air, v* = 3, ll = −185.65
# # DOES NOT WORK
# (table4 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("train", "bus", "car", "air"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 3.0
# ))
#
# # j_0=bus, v* = 30, ll =  −183.79
# (table4 <- GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("air", "train", "car", "bus"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 30
# ))
#
# # j_0=car, v* = 0.2, ll = −142
# ## DOES NOT WORK
# (train_1.35 <- dist3$GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("air", "train", "bus", "car"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 0.7
# ))
# # j_0=car, v* = 1.0, ll = −142
# (train_1.35 <- dist3$GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("air", "train", "bus", "car"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 1
# ))
#
# # j_0=train, v* = 1.35, ll = −183.49
# (train_1.35 <- dist3$GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   depend_y = c("gc", "ttme"),
#   distribution = "student", categories_order = c("air", "bus", "car", "train"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 1.35
# ))
#
# (table4 <- dist3$GLMref_ec(
#   response = "choice", actual_response = "mode",
#   individuals = "indv",
#   explanatory_complete = c("intercept"),
#   depend_y = c("ttme"),
#   distribution = "logistic", categories_order = c("air", "train", "bus", "car"), dataframe = travel.long88,
#   design = "louviere", freedom_degrees = 2.0
# ))
# table4$`Log-likelihood`
# table4$Coefficients
#
#
#
#
# # ANOTHER EXAMPLE FOR ECONOMETRIC MODEL -----------------------------------
#
#
# data("Heating", package = "Ecdat")
#
# # Heating is a "horizontal" data.frame with three choice-specific
# # variables (ic: investment cost, oc: operating cost) and some
# # individual-specific variables (income, region, rooms)
#
# H <- mlogit.data(Heating, shape = "wide", choice = "depvar", varying = c(3:12))
# head(H)
#
# (mi2 <- mlogit(depvar ~ oc + ic | income, H, reflevel = "hp"))
# logLik(mi2)
#
# choice <- sub(".*\\.", "", rownames(H))
# indv <- sub("\\..*", "", rownames(H))
#
# dat4 <- as.data.frame(cbind(indv, choice, H))
# head(dat4, 8)
# str(dat4)
# dist3 <- new(ReferenceF)
# (A98 <- dist3$GLMref_ec(
#   response = "choice", actual_response = "depvar",
#   individuals = "indv",
#   explanatory_complete = c("intercept", "income"),
#   depend_y = c("oc", "ic"),
#   distribution = "logistic",
#   categories_order = c("ec", "er", "gc", "gr", "hp"),
#   dataframe = dat,
#   design = "tutz"
# ))
#
#
# # JSS ---------------------------------------------------------------------
#
# {
#   library(mlogit)
#   data(ModeChoice, package = "Ecdat")
#   head(ModeChoice)
#   travel.long <- mlogit.data(ModeChoice, choice = "mode", shape = "long", alt.levels = c("air", "train", "bus", "car"))
#   head(travel.long)
#   choice <- sub(".*\\.", "", rownames(travel.long))
#   indv <- sub("\\..*", "", rownames(travel.long))
#   travel.long88 <- as.data.frame(cbind(indv, choice, travel.long))
# }
# head(travel.long88, 5)
# str(travel.long88)
#
# travel_dat1 <- travel.long88[travel.long88$mode == T, ]
# head(travel_dat1)
#
#
# # e : REFERENCE, LOGISTIC, COMPLETE
# (l_1 <- GLMref(
#   response = "choice",
#   explanatory_complete = c("intercept", "psize"),
#   explanatory_proportional = c("hinc"),
#   distribution = "logistic",
#   categories_order = c("air", "train", "bus", "car"),
#   dataframe = travel_dat1
# ))
#
# (l_1 <- GLMcum(
#   response = "choice",
#   explanatory_complete = c("intercept", "psize"),
#   explanatory_proportional = c("hinc"),
#   distribution = "normal",
#   categories_order = c("air", "train", "bus", "car"),
#   dataframe = travel_dat1
# ))
#
# (l_1 <- GLMcum(
#   response = "choice",
#   explanatory_complete = c("intercept", "psize"),
#   explanatory_proportional = c("hinc"),
#   distribution = "cauchit",
#   categories_order = c("air", "train", "bus", "car"),
#   dataframe = travel_dat1
# ))
#
# (l_1 <- GLMcum(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc"),
#   distribution = "gompertz",
#   categories_order = c("air", "train", "bus", "car"),
#   dataframe = travel_dat1
# ))
#
#
# # e : REFERENCE, LOGISTIC, PROPORTIONAL
# (l_2 <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "logistic",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# # e : REFERENCE, CAUCHIT, COMPLETE
# (l_3 <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "cauchit",
#   categories_order = c("air", "train", "bus", "car"),
#   dataframe = travel_dat1
# ))
#
# # Then we change the reference category and estimate again the three reference models:
#
# # e : REFERENCE, LOGISTIC, COMPLETE
# (l_1prime <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "logistic",
#   categories_order = c("air", "train", "car", "bus"),
#   dataframe = travel_dat1
# ))
#
#
# # e : REFERENCE, LOGISTIC, PROPORTIONAL
# (l_2prime <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "logistic",
#   categories_order = c("air", "train", "car", "bus"),
#   dataframe = travel_dat1
# ))
#
# # e : REFERENCE, CAUCHIT, COMPLETE
# (l_3prime <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "cauchit",
#   categories_order = c("air", "train", "car", "bus"),
#   dataframe = travel_dat1
# ))
#
#
# # The log-likelihoods l_1 and l_1prime are equal (since the canonical model is invariant under
# # all permutations) whereas the log-likelihoods l_2 and l_2prime are different (respectively
# # l_3 and l_3prime).
#
#
# # 4. Adjacent models for ordinal response
#
# # Equivalence between (adjacent, logistic, complete) and (reference, logistic, complete) models.
#
# # e : REFERENCE, LOGISTIC, COMPLETE
# dist1 <- new(ReferenceF)
# (estimation <- dist1$GLMref(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("hinc"),
#   distribution = "logistic",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# # e : ADJACENT, LOGISTIC, COMPLETE
# dist2 <- new(AdjacentR)
# (estimation_prime <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "logistic",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
#
#
# # Try equality with defined matrix A
#
# # Invariance under permutations
#
# # ADJACENT, CAUCHY, COMPLETE
#
# (estimation_1 <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "cauchit",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# # ADJACENT, GOMPERTZ, PROPORTIONAL
# (estimation_2 <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "gompertz",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# # ADJACENT, CAUCHY, COMPLETE (Reverse O)
#
# (estimation_1r <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("intercept", "hinc", "psize"),
#   explanatory_proportional = c("NA"),
#   distribution = "cauchit",
#   categories_order = c("train", "car", "bus", "air"),
#   dataframe = travel_dat1
# ))
#
# # ADJACENT, GOMPERTZ, PROPORTIONAL (Reverse O)
# (estimation_2r <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "gompertz",
#   categories_order = c("train", "car", "bus", "air"),
#   dataframe = travel_dat1
# ))
#
# # ADJACENT, GUMBEL, PROPORTIONAL (Reverse O)
# (estimation_3r <- dist2$GLMadj(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "gumbel",
#   categories_order = c("train", "car", "bus", "air"),
#   dataframe = travel_dat1
# ))
#
# # Remark that the log-likelihoods l_1 and l_1r are equal since the Cauchy distribution is sym-
# # metric whereas the log-likelihoods l_2 and l_2r are different since the Gompertz distribution
# # is not symmetric. Moreover if the Gumbel distribution is used we the reverse order then
# # the log-likelihoods l_2 and l_3r are equal since the Gumbel distribution is the symmetric
# # of the Gompertz distribution.
#
# # 5. Cumulative models for ordinal response
#
# # The equivalence between the (cumulative, Gompertz, proportional)
# # and (sequential, Gompertz, proportional)
# # models has been demonstrated by Läärä and Matthews (1985).
#
# # Cumulative, Gompertz, Proportional
# (estimation <- GLMcum(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "gompertz",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# dat_lok_1 <- data.frame("LogLikIter" = estimation$LogLikIter, "iteration" = 1:nrow(dat_lok_1))
#
# library(ggplot2)
# ggplot(dat_lok_1[-1, ], aes(x = iteration, y = LogLikIter)) +
#   geom_point()
#
#
# # Sequential, Gompertz, Proportional
# dist2 <- new(SequentialR)
# (estimation_prime <- dist2$GLMseq(
#   response = "choice",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept", "hinc", "psize"),
#   distribution = "gompertz",
#   categories_order = c("air", "bus", "car", "train"),
#   dataframe = travel_dat1
# ))
#
# # Severity of Disturbed Dreams (Anderson) ---------------------------------
# {
#   dreams_d <- read.csv("~/Desktop/Test package/data/Severity of Disturbed Dreams.csv")
#   head(dreams_d)
#
#   # Wide to long
#   library(tidyr)
#   dreams_d1 <- gather(dreams_d, Level, Total, Not.severe:Very.severe)
#
#   # Grouped to ungrouped
#   library(vcdExtra)
#   dreams_d1 <- expand.dft(dreams_d1, freq = "Total")
#   head(dreams_d1)
#   summary(dreams_d1)
#
#   summary.fastLm <- function(object, ...) {
#     coef <- object$coefficients
#     se <- object$stderr
#     tval <- coef / se
#
#     object$coefficients <- cbind(
#       "Estimate" = coef,
#       "Std. Error" = se,
#       "z value" = tval,
#       "Pr(>|z|)" = 2 * pnorm(-abs(tval))
#     )
#     colnames(object$coefficients) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
#
#     # cf src/stats/R/lm.R and case with no weights and an intercept
#     # f <- object$fitted.values
#     # r <- object$residuals
#     # mss <- sum((f - mean(f))^2)
#     # mss <- if (object$intercept) sum((f - mean(f))^2) else sum(f^2)
#     # rss <- sum(r^2)
#     #
#     # object$r.squared <- mss/(mss + rss)
#     # df.int <- if (object$intercept) 1L else 0L
#     # n <- length(f)
#     # rdf <- object$df
#     # object$adj.r.squared <- 1 - (1 - object$r.squared) * ((n - df.int)/rdf)
#     class(object) <- "summary.fastLm"
#     object
#   }
#
#   dist1 <- new(ReferenceF)
#   all_permutations <- permutations(v = c(0, 1, 2, 3, 4), repeats.allowed = F, n = 5, r = 5)
#   Log_lik_Vec <- NA
#   all_permutations <- permutations(
#     v = c("Very.severe", "Severe.1", "Severe.2", "Not.severe"),
#     repeats.allowed = F, n = 4, r = 4
#   )
#
# }
# # REFERENCE, LOGISTIC, COMPLETE
#
# summary.pcglm <- function(object, ...) {
#   coef <- object$coefficients
#   se   <- object$stderr
#   tval <- coef/se
#
#   object$coefficients <- cbind("Estimate"     = coef,
#                                "Std. Error" = se,
#                                "z value"    = tval,
#                                "Pr(>|z|)"   = 2*pnorm(-abs(tval)))
#   colnames(object$coefficients) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
#   printCoefmat(object$coefficients, P.values=TRUE, has.Pvalue=TRUE, ...)
#   # cf src/stats/R/lm.R and case with no weights and an intercept
#   # f <- object$fitted.values
#   # r <- object$residuals
#   #mss <- sum((f - mean(f))^2)
#   # mss <- if (object$intercept) sum((f - mean(f))^2) else sum(f^2)
#   # rss <- sum(r^2)
#   #
#   # object$r.squared <- mss/(mss + rss)
#   # df.int <- if (object$intercept) 1L else 0L
#   # n <- length(f)
#   # rdf <- object$df
#   # object$adj.r.squared <- 1 - (1 - object$r.squared) * ((n - df.int)/rdf)
#   class(object) <- "summary.pcglm"
#   object
# }
# # MODELS
# {
#   (l_1 <- GLMref(
#     response = "Level",
#     explanatory_complete = c("intercept", "Age"),
#     explanatory_proportional = c("NA"),
#     distribution = "logistic",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1
#   ))
#
#   summary.pcglm(l_1)$coefficients[1]
#
#
#   (l_1 <- GLMref(
#     response = "Level",
#     explanatory_complete = c("intercept", "Age"),
#     explanatory_proportional = c("NA"),
#     distribution = "student",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1,
#     30
#   ))
#
#   summary.fastLm(l_1)$coefficients
#   l_1$deviance
#   l_1$`Log-likelihood`
#
#   # REFERENCE, LOGISTIC, PROPORTIONAL
#
#   l_2 <- GLMref(
#     response = "Level",
#     explanatory_complete = c("NA"),
#     explanatory_proportional = c("intercept", "Age"),
#     distribution = "logistic",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(l_2)$coefficients
#   l_2$deviance
#   l_2$`Log-likelihood`
#
#   # REFERENCE, CAUCHIT, COMPLETE
#
#   l_3 <- GLMref(
#     response = "Level",
#     explanatory_complete = c("intercept", "Age"),
#     explanatory_proportional = c("NA"),
#     distribution = "cauchit",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(l_3)$coefficients
#   l_3$deviance
#   l_3$`Log-likelihood`
#
#   # Then we change the reference category (Severe.2) and estimate again the three reference models:
#
#   l_1prime <- GLMref(
#     response = "Level",
#     explanatory_complete = c("intercept", "Age"),
#     explanatory_proportional = c("NA"),
#     distribution = "logistic",
#     categories_order = c("Not.severe", "Severe.1", "Very.severe", "Severe.2"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(l_1prime)$coefficients
#   l_1prime$deviance
#   l_1prime$`Log-likelihood`
#
#   # REFERENCE, LOGISTIC, PROPORTIONAL
#   l_2prime <- GLMref(
#     response = "Level",
#     explanatory_complete = c("NA"),
#     explanatory_proportional = c("intercept", "Age"),
#     distribution = "logistic",
#     categories_order = c("Not.severe", "Severe.1", "Very.severe", "Severe.2"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(l_2prime)$coefficients
#   l_2prime$deviance
#   l_2prime$`Log-likelihood`
#
#   # REFERENCE, CAUCHIT, COMPLETE
#   l_3prime <- GLMref(
#     response = "Level",
#     explanatory_complete = c("intercept", "Age"),
#     explanatory_proportional = c("NA"),
#     distribution = "cauchit",
#     categories_order = c("Not.severe", "Severe.1", "Very.severe", "Severe.2"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(l_3prime)$coefficients
#   l_3prime$deviance
#   l_3prime$`Log-likelihood`
#
#   # ADJACENT, LOGISTIC, COMPLETE
#   lprime <- GLMadj(
#     response = "Level",
#     explanatory_complete = c("intercept", "Age"),
#     explanatory_proportional = c("NA"),
#     distribution = "logistic",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(lprime)$coefficients
#   lprime$deviance
#   lprime$`Log-likelihood`
#
#   # ADJACENT, CAUCHY, COMPLETE
#   estimation_1 <- GLMadj(
#     response = "Level",
#     explanatory_complete = c("intercept", "Age"),
#     explanatory_proportional = c("NA"),
#     distribution = "cauchit",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(estimation_1)$coefficients
#   estimation_1$deviance
#   estimation_1$`Log-likelihood`
#
#   # ADJACENT, GOMPERTZ, PROPORTIONAL
#   estimation_2 <- GLMadj(
#     response = "Level",
#     explanatory_complete = c("NA"),
#     explanatory_proportional = c("intercept", "Age"),
#     distribution = "gompertz",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(estimation_2)$coefficients
#   estimation_2$deviance
#   estimation_2$`Log-likelihood`
#
#   # ADJACENT, CAUCHY, COMPLETE (Reverse order)
#   estimation_1r <- GLMadj(
#     response = "Level",
#     explanatory_complete = c("intercept", "Age"),
#     explanatory_proportional = c("NA"),
#     distribution = "cauchit",
#     categories_order = c("Very.severe", "Severe.2", "Severe.1", "Not.severe"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(estimation_1r)$coefficients
#   estimation_1r$deviance
#   estimation_1r$`Log-likelihood`
#
#   # ADJACENT, GOMPERTZ, PROPORTIONAL (Reverse order)
#   estimation_2r <- GLMadj(
#     response = "Level",
#     explanatory_complete = c("NA"),
#     explanatory_proportional = c("intercept", "Age"),
#     distribution = "gompertz",
#     categories_order = c("Very.severe", "Severe.2", "Severe.1", "Not.severe"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(estimation_2r)$coefficients
#   estimation_2r$deviance
#   estimation_2r$`Log-likelihood`
#
#   # ADJACENT, GUMBEL, PROPORTIONAL (Reverse order)
#   estimation_3r <- GLMadj(
#     response = "Level",
#     explanatory_complete = c("NA"),
#     explanatory_proportional = c("intercept", "Age"),
#     distribution = "gumbel",
#     categories_order = c("Very.severe", "Severe.2", "Severe.1", "Not.severe"),
#     dataframe = dreams_d1
#   )
#   summary.fastLm(estimation_3r)$coefficients
#   estimation_3r$deviance
#   estimation_3r$`Log-likelihood`
#
#   # SEQUENTIAL, GOMPERTZ, PROPORTIONAL
#   (estimation_prime <- ratio_seq$GLMseq(
#     response = "Level",
#     explanatory_complete = c("NA"),
#     explanatory_proportional = c("intercept", "Age"),
#     distribution = "gompertz",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1
#   ))
#
#   (estimation_prime <- GLMseq(
#     response = "Level",
#     explanatory_complete = c("intercept"),
#     explanatory_proportional = c("Age"),
#     distribution = "gompertz",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1
#   ))
#
#   summary.fastLm(estimation_prime)$coefficients
#   estimation_prime$deviance
#   estimation_prime$`Log-likelihood`
#
#   # CUMULATIVE, GOMPERTZ, PROPORTIONAL
#
#   (estimation <- GLMcum(
#     response = "Level",
#     explanatory_complete = c("NA"),
#     explanatory_proportional = c("intercept", "Age"),
#     distribution = "normal",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1,
#     beta_t = c("FALSE"), beta_init = c(-2.1150969, 0.2739375)
#   ))
#
#   (estimation <- GLMcum(
#     response = "Level",
#     explanatory_complete = c("intercept"),
#     explanatory_proportional = c("Age"),
#     distribution = "gompertz",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1,
#     beta_t = c("FALSE"), beta_init = c(-2.1150969, 0.2739375)
#   ))
#
#   # Invariance under permutations
#
#   (l <- GLMcum(
#     response = "Level",
#     explanatory_complete = c("intercept"),
#     explanatory_proportional = c("Age"),
#     distribution = "gompertz",
#     categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#     dataframe = dreams_d1,
#     beta_t = c("FALSE"), beta_init = c(-2.1150969, 0.2739375)
#   ))
#
#   (l <- GLMcum(
#     response = "Level",
#     explanatory_complete = c("intercept"),
#     explanatory_proportional = c("Age"),
#     distribution = "cauchit",
#     categories_order = c("Very.severe", "Severe.2", "Severe.1", "Not.severe"),
#     dataframe = dreams_d1,
#     beta_t = c("FALSE"), beta_init = c(-2.1150969, 0.2739375)
#   ))
#
#   (l <- GLMcum(
#     response = "Level",
#     explanatory_complete = c("intercept"),
#     explanatory_proportional = c("Age"),
#     distribution = "logistic",
#     categories_order = c("Very.severe", "Severe.1", "Severe.2", "Not.severe"),
#     dataframe = dreams_d1,
#     beta_t = c("FALSE"), beta_init = c(-2.1150969, 0.2739375)
#   ))
# }
# # PLOTS
#
# # PLOT ADJACENT
# plot <- list()
#
# for(dist in c("logistic", "normal", "cauchit", "gompertz") ){
#   Log_lik_Vec <- NA
#   for (element in 1:nrow(all_permutations)) {
#     l <- GLMadj(
#       response = "Level",
#       explanatory_complete = c("intercept", "Age"),
#       explanatory_proportional = c("NA"),
#       distribution = dist,
#       categories_order = all_permutations[element, ],
#       dataframe = dreams_d1
#     )
#     Log_lik_Vec[element] <- l$`Log-likelihood`
#   }
#   Log_lik_Vec
#   all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
#   names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
#   to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
#   to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
#   to_plot$Permutation <- as.factor(to_plot$Permutation)
#   to_plot$Distribution <- dist
#   groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
#   to_plot <- left_join(to_plot, groups)
#   to_plot <- to_plot[to_plot$LogLik >= -1000,]
#   title <- str_c("Adjacent, ", dist, ", complete")
#   plot[[dist]] <- to_plot %>%
#     arrange(-LogLik) %>%
#     mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
#     ggplot(aes(x = Permutation, y = LogLik)) +
#     geom_point() +
#     geom_line(aes(group = gn)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     ggtitle(title)+
#     xlab("")+
#     ylab("")
# }
#
#
# plot_adj <- grid.arrange(plot[["logistic"]], plot[["normal"]],
#                          plot[["cauchit"]], plot[["gompertz"]],
#                          ncol=2, nrow=2)
#
# # PLOT ADJACENT Proportional
#
# plot <- list()
#
# for(dist in c("logistic", "normal", "cauchit", "gompertz") ){
#   Log_lik_Vec <- NA
#   for (element in 1:nrow(all_permutations)) {
#     l <- GLMadj(
#       response = "Level",
#       explanatory_complete = c("NA"),
#       explanatory_proportional = c("intercept", "Age"),
#       distribution = dist,
#       categories_order = all_permutations[element, ],
#       dataframe = dreams_d1
#     )
#     Log_lik_Vec[element] <- l$`Log-likelihood`
#   }
#   Log_lik_Vec
#   all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
#   names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
#   to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
#   to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
#   to_plot$Permutation <- as.factor(to_plot$Permutation)
#   to_plot$Distribution <- dist
#   groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
#   to_plot <- left_join(to_plot, groups)
#   to_plot <- to_plot[to_plot$LogLik >= -1000,]
#   title <- str_c("Adjacent, ", dist, ", proportional")
#   plot[[dist]] <- to_plot %>%
#     arrange(-LogLik) %>%
#     mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
#     ggplot(aes(x = Permutation, y = LogLik)) +
#     geom_point() +
#     geom_line(aes(group = gn)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     ggtitle(title)+
#     xlab("")+
#     ylab("")
# }
#
#
# plot_adj_pro <- grid.arrange(plot[["logistic"]], plot[["normal"]],
#                              plot[["cauchit"]], plot[["gompertz"]],
#                              ncol=2, nrow=2)
#
#
# # PLOT CUMULATIVE
#
# plot1 <- list()
#
# for(dist in c("logistic", "normal", "cauchit", "gompertz") ){
#   Log_lik_Vec <- NA
#   for (element in 1:nrow(all_permutations)) {
#     skip_to_next <- FALSE
#     tryCatch({
#
#       l <- GLMcum(
#         response = "Level",
#         explanatory_complete = c("intercept", "Age"),
#         explanatory_proportional = c("NA"),
#         distribution = dist,
#         categories_order = all_permutations[element, ],
#         dataframe = dreams_d1)
#       Log_lik_Vec[element] <- l$`Log-likelihood`
#       all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
#       names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
#       to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
#       to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
#       to_plot$Permutation <- as.factor(to_plot$Permutation)
#       to_plot$Distribution <- dist
#       groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
#       to_plot <- left_join(to_plot, groups)
#       to_plot <- to_plot[to_plot$LogLik >= -1000,]
#       title <- str_c("Cumulative, ", dist, ", complete")
#       plot1[[dist]] <- to_plot %>%
#         arrange(-LogLik) %>%
#         mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
#         ggplot(aes(x = Permutation, y = LogLik)) +
#         geom_point() +
#         geom_line(aes(group = gn)) +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#         ggtitle(title)+
#         xlab("")+
#         ylab("") },error = function(e) {skip_to_next <<- TRUE})
#     if(skip_to_next) { next }
#   }
# }
#
# plot_cum <- grid.arrange(plot1[["logistic"]], plot1[["normal"]],
#                          plot1[["cauchit"]], plot1[["gompertz"]],
#                          ncol=2, nrow=2)
#
# # PLOT ADJACENT Proportional
#
# plot <- list()
#
# for(dist in c("logistic", "normal", "cauchit", "gompertz") ){
#   Log_lik_Vec <- NA
#   for (element in 1:nrow(all_permutations)) {
#     l <- GLMadj(
#       response = "Level",
#       explanatory_complete = c("NA"),
#       explanatory_proportional = c("intercept", "Age"),
#       distribution = dist,
#       categories_order = all_permutations[element, ],
#       dataframe = dreams_d1
#     )
#     Log_lik_Vec[element] <- l$`Log-likelihood`
#   }
#   Log_lik_Vec
#   all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
#   names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
#   to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
#   to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
#   to_plot$Permutation <- as.factor(to_plot$Permutation)
#   to_plot$Distribution <- dist
#   groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
#   to_plot <- left_join(to_plot, groups)
#   to_plot <- to_plot[to_plot$LogLik >= -1000,]
#   title <- str_c("Adjacent, ", dist, ", proportional")
#   plot[[dist]] <- to_plot %>%
#     arrange(-LogLik) %>%
#     mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
#     ggplot(aes(x = Permutation, y = LogLik)) +
#     geom_point() +
#     geom_line(aes(group = gn)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     ggtitle(title)+
#     xlab("")+
#     ylab("")
# }
#
#
# plot_adj_pro <- grid.arrange(plot[["logistic"]], plot[["normal"]],
#                              plot[["cauchit"]], plot[["gompertz"]],
#                              ncol=2, nrow=2)
#
#
# # REFERENCE PLOT
#
# plot_r <- list()
#
# for(dist in c("logistic", "normal", "cauchit", "gompertz") ){
#   Log_lik_Vec <- NA
#   for (element in 1:nrow(all_permutations)) {
#     skip_to_next <- FALSE
#     tryCatch({
#
#       l <- GLMcum(
#         response = "Level",
#         explanatory_complete = c("intercept", "Age"),
#         explanatory_proportional = c("NA"),
#         distribution = dist,
#         categories_order = all_permutations[element, ],
#         dataframe = dreams_d1)
#       Log_lik_Vec[element] <- l$`Log-likelihood`
#       all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
#       names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
#       to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
#       to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
#       to_plot$Permutation <- as.factor(to_plot$Permutation)
#       to_plot$Distribution <- dist
#       groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
#       to_plot <- left_join(to_plot, groups)
#       to_plot <- to_plot[to_plot$LogLik >= -1000,]
#       title <- str_c("Cumulative, ", dist, ", complete")
#       plot_r[[dist]] <- to_plot %>%
#         arrange(-LogLik) %>%
#         mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
#         ggplot(aes(x = Permutation, y = LogLik)) +
#         geom_point() +
#         geom_line(aes(group = gn)) +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#         ggtitle(title)+
#         xlab("")+
#         ylab("") },error = function(e) {skip_to_next <<- TRUE})
#     if(skip_to_next) { next }
#   }
# }
#
# plot_ref <- grid.arrange(plot_r[["logistic"]], plot_r[["normal"]],
#                          plot_r[["cauchit"]], plot_r[["gompertz"]],
#                          ncol=2, nrow=2)
#
#
#
# # KNEE DATASET ------------------------------------------------------------
# # DEMOSTRAR EQUIVALENCIA CON OTRO DATASET
#
# # CUMULATIVE SI PARECE ESTAR BIEN DEFINIDO,
# # PERO QUIZAS LA FORMA EN QUE SE INICIALIZA NO ESTA BIEN
#
# {
#   library(catdata)
#   data(knee)
#   knee$R4 <- as.factor(knee$R4)
# }
#
# estimation <- GLMcum(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "logistic",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee,
#   beta_t = c("FALSE"), beta_init = c(-0.08, 0.95, -0.84, -0.46, -0.18)
# )
#
#
# # NO ES ESTIMABLE CON EL OTRO MODELO
# (estimation <- dist1$GLMcum(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee
# ))
#
# (estimation_prime <- dist2$GLMseq(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee
# ))
#
# (estimation <- dist1$GLMcum(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "logistic",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee
# ))
#
# (estimation_prime <- dist2$GLMseq(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("S·", "Sex", "Age"),
#   distribution = "gumbel",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee
# ))
#
# dat <- read.csv2("~/Desktop/Test package/dat_cum.csv")
# summary(dat)
# str(dat)
# head(dat)
#
# # FUNCIONA BIEN
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("NA"),
#   distribution = "logistic",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("NA"),
#   distribution = "gompertz",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # MOD3
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("public"),
#   distribution = "gompertz",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # MOD4
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("public"),
#   distribution = "normal",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # MOD5
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("pared", "gpa"),
#   explanatory_proportional = c("intercept", "public"),
#   distribution = "normal",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# # MOD6
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("public", "gpa"),
#   distribution = "normal",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat, beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # MOD7
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("public"),
#   distribution = "cauchit",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat, beta_t = c("F"), beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# # MOD 7
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("public"),
#   distribution = "cauchit",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# (estimation_prime <- GLMcum(
#   response = "apply",
#   explanatory_complete = c("intercept", "pared"),
#   explanatory_proportional = c("NA"),
#   distribution = "cauchit",
#   categories_order = c("unlikely", "somewhat likely", "very likely"),
#   dataframe = dat, beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
#
#
# # ORDINAL LIBRARY ---------------------------------------------------------
# # DATA
# {
#   library("ordinal")
#   head(wine)
#   wine$temp <- as.numeric((wine$temp))
#   wine$contact <- as.numeric((wine$contact))
#   str(wine)
# }
#
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("NA"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # NO
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept", "temp"),
#   explanatory_proportional = c("NA"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# # NO
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept", "temp", "contact"),
#   explanatory_proportional = c("NA"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp"),
#   distribution = "cauchit",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp"),
#   distribution = "logistic",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("F"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
#
# # e : ADJACENT, LOGISTIC, COMPLETE
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp", "contact"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine,
#   beta_t = c("TRUE", "sd"),
#   beta_init = c(0.72, 2.76, 4.19, 5.06, -1.6, -0.85)
# ))
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp", "contact"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine
# ))
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp", "contact"),
#   distribution = "cauchit",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine
# ))
#
#
# (estimation_prime <- GLMcum(
#   response = "rating",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("temp", "contact"),
#   distribution = "cauchit",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = wine
# ))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # CUMULATIVE RATIO --------------------------------------------------------
# {
#   library(catdata)
#   data(knee)
#   knee$R4 <- as.factor(knee$R4)
# }
#
# # FUNCIONA
# (estimation <- GLMcum(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "logistic",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee,
#   beta_t = c("FALSE"), beta_init = c(-0.08, 0.95, -0.84, -0.46, -0.18)
# ))
#
# # FUNCIONA
# (estimation <- GLMcum(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "normal",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee,
#   beta_t = c("FALSE"), beta_init = c(-0.08, 0.95, -0.84, -0.46, -0.18)
# ))
#
# # FUNCIONA
# (estimation <- GLMcum(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "cauchit",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee,
#   beta_t = c("FALSE"), beta_init = c(-0.08, 0.95, -0.84, -0.46, -0.18)
# ))
#
# # FUNCIONA
# (estimation <- GLMcum(
#   response = "R4",
#   explanatory_complete = c("intercept"),
#   explanatory_proportional = c("Th", "Sex", "Age"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee,
#   beta_t = c("FALSE"), beta_init = c(-0.08, 0.95, -0.84, -0.46, -0.18)
# ))
#
# (estimation <- GLMcum(
#   response = "R4",
#   explanatory_complete = c("NA"),
#   explanatory_proportional = c("intercept"),
#   distribution = "gompertz",
#   categories_order = c("1", "2", "3", "4", "5"),
#   dataframe = knee,
#   beta_t = c("FALSE"), beta_init = c(-0.08, 0.95, -0.84, -0.46, -0.18)
# ))
#
#
#
#
#
#
#
# # CREATE VIGNETTE ---------------------------------------------------------
# # library(rmarkdown); library(devtools)
# #
# # usethis::use_vignette("my-vignette")
#
#
# # pkgdir("~/Desktop/pack_200320")
