#
# # Libraries ---------------------------------------------------------------
# library(devtools);library(tidyverse);library(fs);library(Rcpp);library(RcppArmadillo);library(RcppEigen)
# library(dplyr); library(tidyr); library(nnet); library(varhandle);library(catdata); library(magic)
# library(VGAM) ; library(dobson)
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
# #error "The file 'Rcpp.h' should not be included. Please correct to include only 'RcppArmadillo.h'."
# # Borre del distribution.h el include rcpp
# getLoadedDLLs()
# # Using modules
# # We need to tell R that we want to use a C++11 compiler
# Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# pkgbuild::compile_dll()
# compileAttributes()
# load_all()
#
# # Binomial case -----------------------------------------------------------
#
# # # DATA
# {
#   beetle
#   i=1
#   beetle_ext <- data.frame(x = as.matrix(rep(beetle[i,1],beetle[i,2])), y= c(rep(1,beetle[i,3]), rep(0,beetle[i,2]-beetle[i,3]) ))
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
#
#   df_beetle <- data.frame(t(data.frame(matrix(unlist(beetle_ext), nrow=length(beetle_ext), byrow=T))))
#   colnames(df_beetle) = c("x_beetle", "y_beetle")
#
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
#   glm(beetle.mat~x, family = binomial(link = "normal"))
#   glm(beetle.mat~x, family = binomial(link = "cauchit"))
#   glm(beetle.mat~x, family = binomial(link = "cloglog"))
# }
#
# # Logit
# FisherScoring
# dist3 <- new(FisherScoring)
# dist3$GLMm(X_M = X,
#            Y_M = Y,
#            link = "logistic")
#
#
# df_beetle$y_beetle <- as.factor(df_beetle$y_beetle)
#
# library(plyr)
# df_beetle$y_beetle <- revalue(df_beetle$y_beetle, c("1"="one", "0"="zero"))
#
# str((df_beetle))
# summary(df_beetle)
#
# dim(df_beetle)
#
# dist1 <- new(ReferenceF)
# dist1$GLMref(response = "y_beetle",
#              explanatory_complete = c("intercept","x_beetle"),
#              explanatory_proportional = NA,
#              distribution = "logistic",
#              categories_order = c("zero", "one"),
#              dataframe = df_beetle )
#
# # For same results than previous one
# dist1 <- new(ReferenceF)
# dist1$GLMref(response = "y_beetle",
#              explanatory_complete = NA,
#              explanatory_proportional = c("intercept","x_beetle"),
#              distribution = "logistic",
#              categories_order = c(0,1),
#              dataframe = df_beetle )
#
# # normal
# dist3 <- new(FisherScoring)
# dist3$GLMm(X_M = as.matrix(X),
#            Y_M = as.matrix(Y),
#            link = "normal")
#
# # Cauchit
# dist3 <- new(FisherScoring)
# dist3$GLMm(X_M = X,
#            Y_M = Y,
#            link = "cauchit")
#
# # Student 2 freedom degrees
# FisherScoring
# dist3 <- new(FisherScoring)
# dist3$GLMm(X_M = X,
#            Y_M = Y,
#            link = "student")
#
# # anadir otro parametro para los grados de libertad
#
# # Gumbel
# dist3$GLMm(X_M = X,
#            Y_M = Y,
#            link = "gumbel")
#
# # Gompertz
# dist3$GLMm(X_M = X,
#            Y_M = Y,
#            link = "gompertz")
#
#
# # Multicategorical response -----------------------------------------------
# # DATA
# {
#   data(addiction)
#   summary(addiction)
#   head(addiction)
#   # vignette("multinomial-addiction1")
#   attach(addiction)
#   ill <- as.factor(ill)
#   addiction$ill<-as.factor(addiction$ill)
#   data1 <- addiction[,c("ill","gender","university", "age")]
#   data2 <- na.omit(data1)
#
# }
# summary(data2)
# colnames(data2)
# str(data2) # RESPONSE FACTOR. COV AS INT
#
# dist1 <- new(ReferenceF)
# (mod1 <- dist1$GLMref(response = "ill", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#                       distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 )) # -724.8664
# (mod2 <- dist1$GLMref(response = "ill", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#                       distribution = "logistic", categories_order = c("0","2","1"), dataframe = data2 )) # -724.8664
# (mod3 <- dist1$GLMref(response = "ill", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#                       distribution = "logistic", categories_order = c("1","2","0"), dataframe = data2 )) # -724.8664
# (mod4 <- dist1$GLMref(response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#                       distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 )) # -737.5639
# (mod5 <- dist1$GLMref(response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#                       distribution = "logistic", categories_order = c("0","2","1"), dataframe = data2 )) # -726.4392
# (mod6 <- dist1$GLMref(response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#                       distribution = "logistic", categories_order = c("1","2","0"), dataframe = data2 )) # -746.9176
#
# (mod7 <- dist1$GLMref(response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#                       distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 )) #  -730.6547
# (mod8 <- dist1$GLMref(response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#                       distribution = "logistic", categories_order = c("0","2","1"), dataframe = data2 )) #  -730.6547
# (mod8 <- dist1$GLMref(response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#                       distribution = "logistic", categories_order = c("2","0","1"), dataframe = data2 )) #
# (mod9 <- dist1$GLMref(response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#                       distribution = "logistic", categories_order = c("1","2","0"), dataframe = data2 )) #
#
# (mod10 <- dist1$GLMref(response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("gender"),
#                        distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 )) # -742.7698
# (mod11 <- dist1$GLMref(response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept", "gender", "university"),
#                        distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 )) # -737.4464
# (mod12 <- dist1$GLMref(response = "ill", explanatory_complete = c("intercept", "gender", "university"), explanatory_proportional = c("NA"),
#                        distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 )) # -698.9851
# # NO LO DEJA EL OTRO
# (mod13 <- dist1$GLMref(response = "ill", explanatory_complete = c("intercept", "gender", "university"), explanatory_proportional = c("age"),
#                        distribution = "logistic", categories_order = c("1","2","0"), dataframe = data2 )) # -698.9851
# # FUNCIONA
# (mod13a <- dist1$GLMref(response = "ill", explanatory_complete = c("intercept","age", "gender", "university"), explanatory_proportional = c("NA"),
#                         distribution = "logistic", categories_order = c("1","2","0"), dataframe = data2 )) # -679.9081
# (mod13a_1 <- dist1$GLMref(response = "ill", explanatory_complete = c("intercept","age", "gender", "university"), explanatory_proportional = c("NA"),
#                           distribution = "logistic", categories_order = c("0","2","1"), dataframe = data2 )) # -679.9081
# (mod13b <- dist1$GLMref(response = "ill", explanatory_complete = c("gender", "university"), explanatory_proportional = c("intercept","age"),
#                         distribution = "logistic", categories_order = c("1","2","0"), dataframe = data2 )) # -679.9081
#
#
# # With categorical variables: Not working with just 2 categories
# data2$gender <- as.factor(data2$gender)
# data2$university <- as.factor(data2$university)
# str(data2)
#
# (mod7 <- dist1$GLMref(response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#                       distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 )) # -732.0838
# (mod8 <- dist1$GLMref(response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#                       distribution = "logistic", categories_order = c("0","2","1"), dataframe = data2 )) # -732.0838
# (mod9 <- dist1$GLMref(response = "ill", explanatory_complete = c("gender"), explanatory_proportional = c("NA"),
#                       distribution = "logistic", categories_order = c("1","2","0"), dataframe = data2 )) # -732.0838
# (mod10 <- dist1$GLMref(response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("gender"),
#                        distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 )) # -742.7698
# (mod11 <- dist1$GLMref(response = "ill", explanatory_complete = c("NA"), explanatory_proportional = c("intercept", "gender", "university"),
#                        distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 )) # -737.4464
# (mod12 <- dist1$GLMref(response = "ill", explanatory_complete = c("intercept", "gender", "university"), explanatory_proportional = c("NA"),
#                        distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 )) # -698.9851
#
#
# (mod12 <- dist1$GLMref(response = "ill", explanatory_complete = c("intercept", "gender", "university"), explanatory_proportional = c("NA"),
#                        distribution = "logistic", categories_order = c("0","1","2"), dataframe = data2 ))
#
# # DATA 2
# Polviews2 <- read.table("http://www.stat.ufl.edu/~aa/cat/data/Polviews2.dat", header=TRUE)
#
# str(Polviews2)
# M2<-as.data.frame(sapply(Polviews2[,c("ideology","party", "gender" )] ,unclass))
# M2 <- as.data.frame(M2)
# sum(complete.cases(M2))
# M2$ideology <- as.factor(M2$ideology)
# str(M2); summary(M2)
#
# dist1 <- new(ReferenceF)
# (mod1_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#                         distribution = "logistic", categories_order = c("1","2", "3", "4", "5"), dataframe = M2 )) # -980.4022
# (mod2_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#                         distribution = "logistic", categories_order = c("1","2", "3", "5","4"), dataframe = M2 )) # -980.4022
# (mod3_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("NA"),
#                         distribution = "logistic", categories_order = c("2", "3", "4", "5","1"), dataframe = M2 )) # -980.4022
# (mod4_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#                         distribution = "logistic", categories_order = c("1","2", "3", "4", "5"), dataframe = M2 )) # -1043.362
# (mod5_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept","party"), explanatory_proportional = c("NA"),
#                         distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2 )) # -977.8485
# (mod6_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept","gender", "party"), explanatory_proportional = c("NA"),
#                         distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2 )) # -774.6079
# (mod7_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("party"),
#                         distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2 )) # -977.8485
# (mod8_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("intercept", "party"),
#                         distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2 )) # -996.6733
#
# # NO FUNCIONA ALLA
# (mod9_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept"), explanatory_proportional = c("party"),
#                         distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2 )) # -977.8485
#
# (mod10_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept","gender"), explanatory_proportional = c("NA"),
#                          distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2 )) # -977.8485
#
# # ALGO RARO TIENE PARTY
# (mod11_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept","party"), explanatory_proportional = c("NA"),
#                          distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2 )) # -775.7563
#
#
#
# (mod6_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("NA"), explanatory_proportional = c("intercept"),
#                         distribution = "logistic", categories_order = c("2", "3", "4", "5", "1"), dataframe = M2 ))
#
# # SI FUNCIONA
# (mod14_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept", "party" ,"gender"), explanatory_proportional = c("NA"),
#                          distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2 )) # -774.6079
# (mod15_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#                          distribution = "logistic", categories_order = c("1", "2", "3", "4", "5"), dataframe = M2 )) # -775.5648
# (mod16_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#                          distribution = "logistic", categories_order =  c("5", "4", "3", "2","1"), dataframe = M2 )) # -775.5648
#
# # TIENE OTRO ORDEN
# (mod16_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#                          distribution = "logistic", categories_order =  c("2", "3", "4", "5","1"), dataframe = M2 )) # -775.5648
# (mod16_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#                          distribution = "logistic", categories_order =  c("5", "2", "3", "4","1"), dataframe = M2 )) # -775.5648
#
# # ALGO RARO CON PARTY
# (mod16_2 <- dist1$GLMref(response = "ideology", explanatory_complete = c("intercept", "party"), explanatory_proportional = c("gender"),
#                          distribution = "logistic", categories_order =  c("5", "2", "3", "4", "1"), dataframe = M2 )) # -775.5648
#
# # GENDER WITHOUT INTERCEPT
#
# # All posible permutations
# a1 = dist1$GLMref(response = "ideology",
#                   explanatory_complete = c("intercept","party"),
#                   explanatory_proportional = "gender",
#                   distribution = "logistic",
#                   categories_order = c(0,1,3,2,4),
#                   dataframe = M2 )
#
# ## Reference
# ### Invariance under permutations
# # Multinomial logit model is invariant under all permutations of the response categories
#
# dist1 <- new(ReferenceF)
# all_permutations = permutations(v=c(0,1,2,3,4),repeats.allowed=F, n = 5, r = 5)
# Log_lik_Vec = NA
# for (element in 1:nrow(all_permutations)){
#   a1 = dist1$GLMref(response = "ideology",
#                     explanatory_complete = c("intercept","party","gender"),
#                     explanatory_proportional = NA,
#                     distribution = "logistic",
#                     categories_order = all_permutations[element,],
#                     dataframe = M2 )
#   Log_lik_Vec[element] = a1$`Log-likelihood`
# }
# Log_lik_Vec
#
#
# # MULTINOMIAL PARTY EX 8.3 TUTZ -------------------------------------------
# # DATA
# {partypref <- matrix(data=c(114, 10, 53, 224,134,9,42,226,114,8,23,174,339,30,13,
#                             414,42,5,44,161,88,10,60,171,90,8,31,168, 413,23,14,375), nrow=8, byrow=TRUE)
# partydat<-data.frame(
#   party=c(rep("CDU",sum(partypref[,1])),rep("SPD",sum(partypref[,4])),
#           rep("The Liberals",sum(partypref[,2])),rep("The Greens",sum(partypref[,3]))),
#   sex=c(rep(0,sum(partypref[1:4,1])),rep(1,sum(partypref[5:8,1])),
#         rep(0,sum(partypref[1:4,4])),rep(1,sum(partypref[5:8,4])),
#         rep(0,sum(partypref[1:4,2])),rep(1,sum(partypref[5:8,2])),
#         rep(0,sum(partypref[1:4,3])),rep(1,sum(partypref[5:8,3]))),
#   age=c(rep(c(1:4,1:4), partypref[,1]),rep(c(1:4,1:4), partypref[,4]),
#         rep(c(1:4,1:4), partypref[,2]),rep(c(1:4,1:4), partypref[,3])))
# partydat$age <- as.factor(partydat$age)}
# head(partydat)
# str(partydat)
# summary(partydat)
#
# (m_party_1 <- dist1$GLMref(response = "party", explanatory_complete = c("intercept", "sex", "age"), explanatory_proportional = NA,
#                            distribution = "logistic", categories_order =  c("SPD", "The Greens", "The Liberals", "CDU"), dataframe = partydat)) # -3521.484
#
# # MULTINOMIAL TRAVEL EX 8.3 TUTZ ------------------------------------------
# # vignette("multinomial-travel")
#
# library(mlogit)
# data(ModeChoice, package="Ecdat")
# head(ModeChoice)
# travel.long <- mlogit.data(ModeChoice, choice="mode", shape="long", alt.levels=c("air","train","bus","car"))
#
# head(travel.long)
# travel.kat.id <- mlogit(mode ~ invt + gc|hinc, data=travel.long)
#
# # Now with VGAM
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
# fit <- vglm(travelmode ~ Invt + Gc + travelhinc, multinomial(parallel = FALSE ~ travelhinc, refLevel = 1),
#             xij = list(Invt ~ invttrain + invtbus + invtcar,Gc ~ gctrain + gcbus + gccar),
#             form2 = ~ Invt + invttrain + invtbus + invtcar + Gc + gctrain + gcbus + gccar + travelhinc,
#             data = traveldat, trace = TRUE)
#
# # DATOS PARA MI EC MODELO -------------------------------------------------
#
# head(travel.long)
# choice <- sub('.*\\.', '', rownames(travel.long))
# indv <- sub("\\..*","", rownames(travel.long))
#
# travel.long87 <- cbind(indv, choice, travel.long)
# travel.long88 <- as.data.frame(travel.long87)
#
# head(travel.long88,8)
# tail(travel.long88)
#
# dist3 <- new(ReferenceF)
# A98 = dist3$GLMref_ec(response = "choice", actual_response = "mode",
#                       individuals = "indv",
#                       explanatory_complete = c("intercept","hinc"),
#                       depend_y = c("invt", "gc"),
#                       distribution = "logistic", categories_order =  c("train", "bus", "car", "air"), dataframe = travel.long88)
# dim(A98$Y_init)
# dim(A98$X_EXT)
# A98$Coefficients
#
# A98 = dist3$GLMref_ec(response = "choice", actual_response = "mode",
#                       individuals = "indv",
#                       explanatory_complete = c("intercept","hinc"),
#                       depend_y = c("invt", "gc"),
#                       distribution = "normal", categories_order =  c("bus", "car", "air", "train"), dataframe = travel.long88)
# dim(A98$Y_init)
# dim(A98$X_EXT)
# A98$Coefficients
#
#
#
# # ANOTHER EXAMPLE FOR ECONOMETRIC MODEL -----------------------------------
#
#
# data("Heating",package="Ecdat")
#
# # Heating is a "horizontal" data.frame with three choice-specific
# # variables (ic: investment cost, oc: operating cost) and some
# # individual-specific variables (income, region, rooms)
#
# H <- mlogit.data(Heating, shape="wide", choice="depvar", varying=c(3:12))
# head(H)
#
# m <- mlogit(depvar~ic+oc|0, H)
#
#
# head(H)
# choice <- sub('.*\\.', '', rownames(H))
# indv <- sub("\\..*","", rownames(H))
#
# dat4 <- cbind(indv, choice, H)
# dat4 <- as.data.frame(dat4)
#
# head(dat4,8)
# str(dat4)
# dist3 <- new(ReferenceF)
#
# mi2 <- mlogit(depvar~oc+ic|income, H, reflevel="hp")
# A98 = dist3$GLMref_ec(response = "choice", actual_response = "depvar",
#                       individuals = "indv",
#                       explanatory_complete = c("intercept","income"),
#                       depend_y = c("oc","ic"),
#                       distribution = "logistic", categories_order =  c("ec", "er", "gc", "gr", "hp"), dataframe = dat4)
# A98$Coefficients
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
