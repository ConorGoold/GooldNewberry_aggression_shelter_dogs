#*********************************************************************************************
# R code for:
# Goold & Newberry (2017) Aggressiveness as a latent personality trait of domestic dogs:
#                         testing local independence and measurement invariance
# prepint: http://biorxiv.org/content/early/2017/03/18/117440
#*********************************************************************************************

# load relevant packages
# NB: install packages if required -- use normal install.packages() unless otherwise stated below!

require(Amelia)
require(psych)
#install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
require(lavaan)
#devtools::install_github("simsem/semTools/semTools")
require(semTools)
require(reshape2)
require(data.table)
# install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
require(loo)
require(ggplot2)
require(rethinking)

#========================================================
# set your working directory!
#========================================================

# setwd("your_path")

#=================================================================================
# load the reliability/validity session data
#=================================================================================

rel_data <- read.csv("Reliability_validity_data.csv")

# proportion for each answer:
#           interacting_people2 = the aggression towards people example
#           kennel2 = the aggression towards dogs example
apply(rel_data[,-1], 2, function(z) table(z) / nrow(rel_data) )

#=================================================================================
# load the behavioual data
#=================================================================================

raw_data <- read.csv("raw_data.csv")

# COLNAMES:  1. day (day since arrival to shelter)
#            2. handling (handling context)
#            3. kennel_assess (in kennel context)
#            4. out_kennel (out of kennel context)
#            5. res_food (interactions around food context)
#            6. res_toys (interactions around toys context)
#            7. soc_known (interacting with known people context)
#            8. soc_unknown (interaction with unfamiliar people context)
#            9. dogs_male (interacting with male dogs context)
#            10. dogs_female (interacting with female dogs context)
#            11. id variable
#            12. sex (1 = males; 2 = females)
#            13. mean_age (average age whilst at the shelter)
#            14. mean_weight (average weight whilst at the shelter)
#            15. site (which shelter site: London, Brands Hatch, Old Windsor)
#            16. neutered (1 = before arrival; 2 = at the shelter; 3 = not)
#            17. total_days (total number of days at the shelter)
#            18. source_type (1 = gift dogs; 2 = returned dogs; 3 = stray dogs)

#========================================================
# organise the data for imputation
#========================================================
# The context variables are coded with numbers representing the behavioural codes for the observations on on each day.
# the following extracts the observations of aggression, denoted by a 1, and makes all others 0

# function to extract aggression code and aggregate the data
getDich <- function(data, x, codeNum) {
  zeroOneNA <- data.frame(values = ifelse(data[,x] == codeNum, 1, ifelse(data[,x] == "NA", NA, 0)),
                          id = data$id)
  aggZeroOneNA <- aggregate( values ~ id, data=zeroOneNA,
                             FUN=function(z) ifelse(all(is.na(z)) == FALSE , sum(z,na.rm=TRUE), NA), na.action = na.pass)
  ifelse(aggZeroOneNA[,2] > 0, 1, ifelse(is.na(aggZeroOneNA[,2]), NA, 0))
}

data_to_impute <- data.frame( aggregate(cbind(sex, mean_age, mean_weight, site, neutered,
                                          total_days, source_type
                                          ) ~ id, data = raw_data, FUN = function(z) head(z,1),
                                    na.action = na.pass
                                    ),

                          handling = getDich(raw_data, "handling", 11),
                          soc_known = getDich(raw_data, "soc_known", 13),
                          soc_unknown = getDich(raw_data, "soc_unknown", 13),

                          kennel_dogs =  getDich(raw_data, "kennel_assess", 9), # split kennel
                          kennel_people = getDich(raw_data, "kennel_assess", 11), # assessment category into 2

                          out_kennel_dogs = getDich(raw_data, "out_kennel", 6), # split out of kennel
                          out_kennel_people = getDich(raw_data, "out_kennel", 10), # category into 2

                          res_food = getDich(raw_data, "res_food", 14),
                          res_toys = getDich(raw_data, "res_toys", 12),
                          dogs_male = getDich(raw_data, "dogs_male", 16),
                          dogs_female = getDich(raw_data, "dogs_female", 16)
                          )

#========================================================
# impute the missing data
#========================================================

# a 3-column matrix of logical bounds for the imputed variables (column number, lower bound, upper bound)
bounds_matrix <- matrix(c(4, min(raw_data$mean_weight, na.rm = T), max(raw_data$mean_weight, na.rm = T) ,
                         6, 1, 3,
                         9, 0, 1,
                         10, 0, 1,
                         11, 0, 1,
                         12, 0, 1,
                         13, 0, 1,
                         14, 0, 1,
                         15, 0, 1,
                         16, 0, 1,
                         17, 0, 1,
                         18, 0, 1,
                         19, 0, 1
                      ),
                      nrow = 13, ncol = 3, byrow = TRUE
                  )

set.seed(12345)
impute_data <- amelia(x=data_to_impute, m=5, cs=c("id"),
                      ords=names(data_to_impute)[c(9:19)], noms = c("neutered"),
                      bounds = bounds_matrix)

# get each imputed data set and re-set column names
impute1 <- as.data.frame(impute_data$imputations[[1]])
colnames(impute1) <- colnames(data_to_impute)
impute2 <- as.data.frame(impute_data$imputations[[2]])
colnames(impute2) <- colnames(data_to_impute)
impute3 <- as.data.frame(impute_data$imputations[[3]])
colnames(impute3) <- colnames(data_to_impute)
impute4 <- as.data.frame(impute_data$imputations[[4]])
colnames(impute4) <- colnames(data_to_impute)
impute5 <- as.data.frame(impute_data$imputations[[5]])
colnames(impute5) <- colnames(data_to_impute)

# complete data frame of all imputations
imputes_all <- cbind(rbind(impute1, impute2, impute3, impute4, impute5),
                     impute_id = rep(1:5,each=4743))

#========================================================================
# Structural equation model / confirmatory factor analysis
#========================================================================

# create a list of imputed data frames
impute_list = split(imputes_all, imputes_all$impute_id)

## set up model in lavaan syntax

two_fac_model <- '

  LV1 =~ handling + kennel_people + out_kennel_people + soc_known + soc_unknown + res_food + res_toys
  LV2 =~ dogs_female + dogs_male + kennel_dogs + out_kennel_dogs
'


varNames <- c("handling", "kennel_people", "out_kennel_people", 'soc_known', "soc_unknown",
              "res_food", "res_toys", "dogs_female", "dogs_male", "kennel_dogs" ,"out_kennel_dogs")

# run CFAs using the runMI() function to combine data sets

two_fac_mod_orthogonal <- runMI(model=two_fac_model, data = impute_list, ordered = varNames, fun="cfa",
                              estimator="WLSMV" , parameterization = "theta", std.lv = TRUE,
                              orthogonal=TRUE)

# check fit indices
anova(two_fac_mod_orthogonal, test = "D2", indices = "all")

two_fac_mod_correlated <- runMI(model=two_fac_model, data = impute_list, ordered = varNames, fun="cfa",
                            estimator="WLSMV" , parameterization = "theta", std.lv = TRUE)

anova(two_fac_mod_correlated, test = "D2", indices = "all")
summary(two_fac_mod_correlated)

# test for local independence between hypothesised variables

LI_model <- '
  # latent variables
  LV1 =~ handling + kennel_people + out_kennel_people + soc_known + soc_unknown + res_food + res_toys
  LV2 =~ dogs_female + dogs_male + kennel_dogs + out_kennel_dogs

  # residual correlations
  handling ~~ kennel_people + soc_known + soc_unknown + res_toys
  out_kennel_people ~~ soc_known + soc_unknown
  soc_known ~~ soc_unknown
  dogs_female ~~ dogs_male
  out_kennel_dogs ~~ dogs_female + dogs_male
'

LI_model_fit <- runMI(model=LI_model, data = impute_list, ordered = varNames, fun="cfa",
                   estimator="WLSMV" , parameterization = "theta", std.lv = TRUE)

anova(LI_model_fit, test = "D2", indices = "all")
summary(LI_model_fit)

#========================================================================
# Compute the Stan hierarchical logistic regressions
#========================================================================

# Organise the data

LV1_data <- imputes_all[, c("impute_id","id","mean_age","mean_weight","site","sex",
                         "neutered","total_days","source_type",
                         "handling","soc_known","soc_unknown",
                         "res_food","res_toys","kennel_people","out_kennel_people")]


LV1_melt <- melt( data = LV1_data ,
                      id.vars = c("impute_id", "id", "mean_age", "sex",
                                  "mean_weight","total_days","site","neutered","source_type") ,
                      measured.vars = c("handling", "kennel_people","out_kennel_people",
                                        "res_food" , "res_toys", "soc_known", "soc_unknown"),
                      value.name = "response")

LV1_melt$age_groups = with(LV1_melt,
                              ifelse(mean_age <= 0.83 , 1,
                                     ifelse(mean_age > 0.83 & mean_age <= 3, 2,
                                            ifelse(mean_age > 3 & mean_age <= 6, 3 , 4 )
                                     )))

LV1_d <- LV1_melt[with(LV1_melt, order(id)), ]
rownames(LV1_d) <- 1:nrow(LV1_d)

LV2_data <- imputes_all[, c("impute_id","id","mean_age","mean_weight","site","sex",
                            "neutered","total_days","source_type",
                         "dogs_female","dogs_male","kennel_dogs","out_kennel_dogs")]

LV2_melt <- melt( data = LV2_data ,
                      id.vars = c("impute_id","id","mean_age","mean_weight","site","sex",
                                  "neutered","total_days","source_type") ,
                      measured.vars = c("dogs_female","dogs_male","kennel_dogs","out_kennel_dogs"),
                      value.name = "response")

LV2_melt$age_groups = with(LV2_melt,
                              ifelse(mean_age <= 0.83 , 1,
                                     ifelse(mean_age > 0.83 & mean_age <= 3, 2,
                                            ifelse(mean_age > 3 & mean_age <= 6, 3 , 4 )
                                     )))

LV2_d <- LV2_melt[with(LV2_melt, order(id)), ]
rownames(LV2_d) <- 1:nrow(LV2_d)

#=================================================================================
# runStanMI: a function for computing the Stan models on the imputed data sets
#            for each of the different models, combining the posterior distributions,
#             and saving the output to a working directory
#=================================================================================


runStanMI <- function( MIdata, modelType, modelFile,
                       nChains, nIter, nWarmup, nCores,
                       displayParameters = c("a","Beta", "sigmaID"),
                       outputSaveName, # file name to hold summary output file
                       mcmcSaveName, # file name to hold MCMC matrix
                       combSaveName  # file name to hold combined MCMC matrix
                       )
  {
  mcmcList <- rep(list(list()), max(MIdata$impute))
  for(i in 1:max(MIdata$impute)){
    getData <- MIdata[MIdata$impute_id==i, ]
    getData$total_daysZ = (getData$total_days - (mean(getData$total_days)))/(2 * sd(getData$total_days))
    getData$weightZ <- (getData$mean_weight - (mean(getData$mean_weight)))/(2 * sd(getData$mean_weight))
    getData[,c("variable","age_groups",
               "sex","neutered","source_type","site")] <- apply(getData[,c("variable","age_groups",
                                                                           "sex","neutered","source_type","site")],
                                                                2,as.factor)

    if (modelType == "fullMI") {
      designMatrix <- model.matrix( ~ variable * age_groups + variable * sex +
                                      total_daysZ + weightZ + neutered + source_type + site ,
                                    data = getData,
                                    contrasts.arg = list(variable = "contr.sum",age_groups="contr.sum",sex="contr.sum",
                                                         neutered = "contr.sum",site = "contr.sum", source_type = "contr.sum") )
    }

    if (modelType == "ageMI") {
      designMatrix <- model.matrix( ~ variable * age_groups + sex +
                                      total_daysZ + weightZ + neutered + source_type + site ,
                                    data = getData,
                                    contrasts.arg = list(variable = "contr.sum",age_groups="contr.sum",sex="contr.sum",
                                                         neutered = "contr.sum",site = "contr.sum", source_type = "contr.sum") )
    }

    if (modelType == "sexMI") {
      designMatrix <- model.matrix( ~ age_groups + variable * sex +
                                      total_daysZ + weightZ + neutered + source_type + site ,
                                    data = getData,
                                    contrasts.arg = list(variable = "contr.sum",age_groups="contr.sum",sex="contr.sum",
                                                         neutered = "contr.sum",site = "contr.sum", source_type = "contr.sum") )
    }

    if (modelType == "noMI") {
      designMatrix <- model.matrix( ~ variable + age_groups + sex +
                                      total_daysZ + weightZ + neutered + source_type + site ,
                                    data = getData,
                                    contrasts.arg = list(variable = "contr.sum",age_groups="contr.sum",sex="contr.sum",
                                                         neutered = "contr.sum",site = "contr.sum", source_type = "contr.sum") )
    }

    dataList <- list( y = getData$response, N = nrow(getData), K = ncol(designMatrix[,-1]),
                      id = getData$id, Nid = length(unique(getData$id)),
                      X = designMatrix[,-1])

    startTime = proc.time()
    fitStan <- stan(file = modelFile, data = dataList,
                    chains = nChains, warmup = nWarmup, iter = nIter, cores = nCores)
    proc.time() - startTime
    print(paste("*******Complete Stan model of imputation ", i, "...saving summary***********", sep = ""))
    capture.output(print(fitStan, pars=displayParameters, digits=3, probs = c(0.025,0.975)),
                   file = paste(outputSaveName, i, ".csv", sep=""))
    print(paste("*******Saving MCMC matrix of imputation ", i, sep = ""))
    fwrite(x = as.data.frame(fitStan), file = paste(mcmcSaveName, i, ".csv", sep=""))
    print(paste("*******Combining the list of MCMC imputations***********", sep = ""))
    mcmcList[[i]] <- as.matrix(fitStan)
    print(paste("*******Complete Stan model of imputation ", i, "***********", sep = ""))
  }
  print(paste("*******Saving combined MCMC matrix****************", sep = ""))
  fwrite(x = as.data.frame(Reduce("+", mcmcList) / length(mcmcList)), paste(combSaveName, ".csv", sep=""), row.names = FALSE)
}

# which models to compute (all of them)
which_models <- c("fullMI", # testing age and sex measurement invariance
                  "ageMI", # testing age measurement invariance
                  "sexMI", # testing sex measurement invariance
                  "noMI" # assuming complete measurement invariance
                  )

#==================================================
# run the runStanMI() function on the imputed data sets for the specified models.
# with all the models, there are 40 models to run on each of the
# aggressiveness towards people and dogs, respectively, variables

# aggressiveness towards people
for( m in which_models) {
  setwd("")  # set working directory
  runStanMI(MIdata = LV1_d, modelFile = "Stan_hier_logistic.stan", modelType = m,
            nChains = 4, nIter = 2000, nWarmup = 1000, nCores = 4,
            displayParameters = c("a","Beta","sigmaID"),
            outputSaveName = paste("StanSummary_LV1_", m, "_", sep=""),
            mcmcSaveName = paste("StanMCMC_LV1_", m, "_", sep=""),
            combSaveName = paste("combinedMCMC_LV1_", m, sep=""))
  print(paste("***********Complete ", m, "...on to the next********************", sep=""))
}

# aggressiveness towards dogs
for( m in which_models) {
  setwd("")  # set working directory
  runStanMI(MIdata = LV2_d, modelFile = "Stan_hier_logistic.stan", modelType = m,
            nChains = 4, nIter = 2000, nWarmup = 1000, nCores = 4,
            displayParameters = c("a","Beta","sigmaID"),
            outputSaveName = paste("StanSummary_LV2_", m, "_", sep=""),
            mcmcSaveName = paste("StanMCMC_LV2_", m, "_", sep=""),
            combSaveName = paste("combinedMCMC_LV2_", m, sep=""))
  print(paste("***********Complete ", m, "...on to the next********************", sep=""))
}

#=================================================================================
# Post-processing of Stan results
#=================================================================================

#=================================================================================
# LV1 data
#=================================================================================

# load MCMC matrices, get log-likelihoods and get WAIC estimates
fullMI <- as.data.frame(fread("CombinedMCMC_LV1_fullMI.csv"))
ageMI <- as.data.frame(fread("CombinedMCMC_LV1_ageMI.csv"))
sexMI <- as.data.frame(fread("CombinedMCMC_LV1_sexMI.csv"))
noMI <- as.data.frame(fread("CombinedMCMC_LV1_noMI.csv"))

fullMI_waic <- waic(as.matrix(fullMI[, grep("log_lik",colnames(fullMI))]))
ageMI_waic <- waic(as.matrix(ageMI[, grep("log_lik",colnames(ageMI))]))
sexMI_waic <- waic(as.matrix(sexMI[, grep("log_lik",colnames(sexMI))]))
noMI_waic <- waic(as.matrix(noMI[, grep("log_lik",colnames(noMI))]))

WAICtable_LV1 <- data.frame( model = c("Model 1", "Model 2", "Model 3", "Model 4"),
                             WAIC = c(fullMI_waic$waic,ageMI_waic$waic,sexMI_waic$waic,noMI_waic$waic),
                             SE = c(fullMI_waic$se_waic,ageMI_waic$se_waic,sexMI_waic$se_waic,noMI_waic$se_waic)
                             )

# Now, plot the fixed sex and age relationships as in Figure 1
mcmcMat <- as.data.frame(fread("CombinedMCMC_LV1_fullMI.csv"))
fixedEffects <- mcmcMat[ , c(1, grep("Beta", colnames(mcmcMat)))]
colnames(fixedEffects) <- c( "intercept", "HND", "KP", "OKP", "RF", "RT", "FPL",
                             "age1", "age2", "age3", "male", "totalDays", "weight",
                             "notNeutered", "neuteredOnSite", "gift", "return", "BDCH", "BBH",
                             "HND_age1", "KP_age1", "OKP_age1", "RF_age1", "RT_age1", "FPL_age1",
                             "HND_age2", "KP_age2", "OKP_age2", "RF_age2", "RT_age2", "FPL_age2",
                             "HND_age3", "KP_age3", "OKP_age3", "RF_age3", "RT_age3", "FPL_age3",
                             "HND_male", "KP_male", "OKP_male", "RF_male", "RT_male", "FPL_male")

# convert dummy coded variables (sum coding) so that all estimates for fixed effects are available
fixedEffects$NPL <- with(fixedEffects, (HND + KP + OKP + RF + RT + FPL)*-1)
fixedEffects$age4 <- with(fixedEffects, (age1+age2+age3)*-1)
fixedEffects$female <- with(fixedEffects, male*-1)
fixedEffects$neutered <- with(fixedEffects, (notNeutered+neuteredOnSite)*-1)
fixedEffects$stray <- with(fixedEffects, (gift+return)*-1)
fixedEffects$BOW <- with(fixedEffects, (BDCH+BBH)*-1)
fixedEffects$HND_age4 <- with(fixedEffects, (HND_age1+HND_age2+HND_age3)*-1)
fixedEffects$KP_age4 <- with(fixedEffects, (KP_age1+KP_age2+KP_age3)*-1)
fixedEffects$OKP_age4 <- with(fixedEffects, (OKP_age1+OKP_age2+OKP_age3)*-1)
fixedEffects$RF_age4 <- with(fixedEffects, (RF_age1+RF_age2+RF_age3)*-1)
fixedEffects$RT_age4 <- with(fixedEffects, (RT_age1+RT_age2+RT_age3)*-1)
fixedEffects$FPL_age4 <- with(fixedEffects, (FPL_age1+FPL_age2+FPL_age3)*-1)
fixedEffects$HND_female <- with(fixedEffects, HND_male*-1)
fixedEffects$KP_female <- with(fixedEffects, KP_male*-1)
fixedEffects$OKP_female <- with(fixedEffects, OKP_male*-1)
fixedEffects$RF_female <- with(fixedEffects, RF_male*-1)
fixedEffects$RT_female <- with(fixedEffects, RT_male*-1)
fixedEffects$FPL_female <- with(fixedEffects, FPL_male*-1)
fixedEffects$NPL_age1 <- with(fixedEffects, (HND_age1 + KP_age1 + OKP_age1 + RF_age1 + RT_age1 + FPL_age1)*-1)
fixedEffects$NPL_age2 <- with(fixedEffects, (HND_age2 + KP_age2 + OKP_age2 + RF_age2 + RT_age2 + FPL_age2)*-1)
fixedEffects$NPL_age3 <- with(fixedEffects, (HND_age3 + KP_age3 + OKP_age3 + RF_age3 + RT_age3 + FPL_age3)*-1)
fixedEffects$NPL_age4 <- with( fixedEffects, (NPL_age1+NPL_age2+NPL_age3)*-1)
fixedEffects$NPL_male <- with(fixedEffects, (HND_male + KP_male + OKP_male + RF_male + RT_male + FPL_male)*-1)
fixedEffects$NPL_female <- with(fixedEffects, NPL_male*-1)

# get the predicted log odds of aggression
predictedLogOdds <- with(fixedEffects,
                         data.frame( "HND" = intercept + HND ,"KP" = intercept + KP, "OKP" = intercept + OKP, "RF" = intercept + RF,
                                     "RT" = intercept + RT, "FPL" = intercept + FPL ,"NPL" = intercept + NPL, "age1" = intercept + age1,
                                     "age2" = intercept + age2,"age3" = intercept + age3, "age4" = intercept + age4,
                                     "male" = intercept + male, "female" = intercept + female, "notNeutered" =intercept+notNeutered,
                                     "neuteredOnSite"=intercept+neuteredOnSite, "neutered"=intercept+ neutered,
                                     "gift" = intercept+gift, "return"=intercept+return, "stray"=intercept+stray,
                                     "BDCH"=intercept+BDCH,"BBH"=intercept+BBH,"BOW"=intercept+BOW,
                                     "HND_age1"=intercept+HND+age1+HND_age1, "KP_age1"=intercept+KP+age1+KP_age1,
                                     "OKP_age1"=intercept+OKP+age1+OKP_age1,"RF_age1"=intercept+RF+age1+RF_age1,
                                     "RT_age1"=intercept+RT+age1+RT_age1,"FPL_age1"=intercept+FPL+age1+FPL_age1, "NPL_age1"=intercept+age1+NPL,
                                     "HND_age2"=intercept+HND+age2+HND_age2, "KP_age2"=intercept+KP+age2+KP_age2,
                                     "OKP_age2"=intercept+OKP+age2+OKP_age2,"RF_age2"=intercept+RF+age2+RF_age2,
                                     "RT_age2"=intercept+RT+age2+RT_age2,"FPL_age2"=intercept+FPL+age2+FPL_age2, "NPL_age2"=intercept+NPL+age2+NPL_age2,
                                     "HND_age3"=intercept+HND+age3+HND_age3, "KP_age3"=intercept+KP+age3+KP_age3,
                                     "OKP_age3"=intercept+OKP+age3+OKP_age3,"RF_age3"=intercept+RF+age3+RF_age3,
                                     "RT_age3"=intercept+RT+age3+RT_age3,"FPL_age3"=intercept+FPL+age3+FPL_age3, "NPL_age3"=intercept+NPL+age3+NPL_age3,
                                     "HND_age4"=intercept+HND+age4+HND_age4, "KP_age4"=intercept+KP+age4+KP_age4,
                                     "OKP_age4"=intercept+OKP+age4+OKP_age4,"RF_age4"=intercept+RF+age4+RF_age4,
                                     "RT_age4"=intercept+RT+age4+RT_age4,"FPL_age4"=intercept+FPL+age4+FPL_age4, "NPL_age4"=intercept+NPL+age4+NPL_age4,
                                     "HND_male"=intercept+HND+male+HND_male, "KP_male"=intercept+KP+male+KP_male,
                                     "OKP_male" =intercept+OKP+male+OKP_male, "RF_male"=intercept+RF+male+RF_male,
                                     "RT_male"=intercept+RT+male+RT_male, "FPL_male"=intercept+FPL+male+FPL_male, "NPL_male"=intercept+NPL+male+NPL_male,
                                     "HND_female"=intercept+HND+female+HND_female, "KP_female"=intercept+KP+female+KP_female,
                                     "OKP_female" =intercept+OKP+female+OKP_female, "RF_female"=intercept+RF+female+RF_female,
                                     "RT_female"=intercept+RT+female+RT_female, "FPL_female"=intercept+FPL+female+FPL_female, "NPL_female"=intercept+NPL+female+NPL_female
                         ))

# marginalise over the random effects and convert to probabilities for plotting
c2 <- ((16*sqrt(3))/(15*pi))^2
predictedLogOdds_approx <- apply(predictedLogOdds, 2 , function(x) x/sqrt(1 + c2 * mcmcMat[,"sigmaID"]^2))
predictedProbs <- apply(predictedLogOdds_approx, 2, function(x) logistic(x) )

d <- LV1_d

# get raw probabilities of aggression by age groups and sex
rawProbsAge = data.frame("HND_age1" = sum(d[d$variable=="handling"&d$age_groups==1, "response"])/
                           length(d[d$variable=="handling"&d$age_groups==1, "response"]),
                         "KP_age1" = sum(d[d$variable=="kennel_people"&d$age_groups==1, "response"])/
                           length(d[d$variable=="kennel_people"&d$age_groups==1, "response"]),
                         "OKP_age1" = sum(d[d$variable=="out_kennel_people"&d$age_groups==1, "response"])/
                           length(d[d$variable=="out_kennel_people"&d$age_groups==1, "response"]),
                         "RF_age1"= sum(d[d$variable=="res_food"&d$age_groups==1, "response"])/
                           length(d[d$variable=="res_food"&d$age_groups==1, "response"]),
                         "RT_age1"= sum(d[d$variable=="res_toys"&d$age_groups==1, "response"])/
                           length(d[d$variable=="res_toys"&d$age_groups==1, "response"]),
                         "FPL_age1"= sum(d[d$variable=="soc_known"&d$age_groups==1, "response"])/
                           length(d[d$variable=="soc_known"&d$age_groups==1, "response"]),
                         "NPL_age1"= sum(d[d$variable=="soc_unknown"&d$age_groups==1, "response"])/
                           length(d[d$variable=="soc_unknown"&d$age_groups==1, "response"]),
                         "HND_age2" = sum(d[d$variable=="handling"&d$age_groups==2, "response"])/
                           length(d[d$variable=="handling"&d$age_groups==2, "response"]),
                         "KP_age2" = sum(d[d$variable=="kennel_people"&d$age_groups==2, "response"])/
                           length(d[d$variable=="kennel_people"&d$age_groups==2, "response"]),
                         "OKP_age2" = sum(d[d$variable=="out_kennel_people"&d$age_groups==2, "response"])/
                           length(d[d$variable=="out_kennel_people"&d$age_groups==2, "response"]),
                         "RF_age2"= sum(d[d$variable=="res_food"&d$age_groups==2, "response"])/
                           length(d[d$variable=="res_food"&d$age_groups==2, "response"]),
                         "RT_age2"= sum(d[d$variable=="res_toys"&d$age_groups==2, "response"])/
                           length(d[d$variable=="res_toys"&d$age_groups==2, "response"]),
                         "FPL_age2"= sum(d[d$variable=="soc_known"&d$age_groups==2, "response"])/
                           length(d[d$variable=="soc_known"&d$age_groups==2, "response"]),
                         "NPL_age2"= sum(d[d$variable=="soc_unknown"&d$age_groups==2, "response"])/
                           length(d[d$variable=="soc_unknown"&d$age_groups==2, "response"]),
                         "HND_age3" = sum(d[d$variable=="handling"&d$age_groups==3, "response"])/
                           length(d[d$variable=="handling"&d$age_groups==3, "response"]),
                         "KP_age3" = sum(d[d$variable=="kennel_people"&d$age_groups==3, "response"])/
                           length(d[d$variable=="kennel_people"&d$age_groups==3, "response"]),
                         "OKP_age3" = sum(d[d$variable=="out_kennel_people"&d$age_groups==3, "response"])/
                           length(d[d$variable=="out_kennel_people"&d$age_groups==3, "response"]),
                         "RF_age3"= sum(d[d$variable=="res_food"&d$age_groups==3, "response"])/
                           length(d[d$variable=="res_food"&d$age_groups==3, "response"]),
                         "RT_age3"= sum(d[d$variable=="res_toys"&d$age_groups==3, "response"])/
                           length(d[d$variable=="res_toys"&d$age_groups==3, "response"]),
                         "FPL_age3"= sum(d[d$variable=="soc_known"&d$age_groups==3, "response"])/
                           length(d[d$variable=="soc_known"&d$age_groups==3, "response"]),
                         "NPL_age3"= sum(d[d$variable=="soc_unknown"&d$age_groups==3, "response"])/
                           length(d[d$variable=="soc_unknown"&d$age_groups==3, "response"]),
                         "HND_age4" = sum(d[d$variable=="handling"&d$age_groups==4, "response"])/
                           length(d[d$variable=="handling"&d$age_groups==4, "response"]),
                         "KP_age4" = sum(d[d$variable=="kennel_people"&d$age_groups==4, "response"])/
                           length(d[d$variable=="kennel_people"&d$age_groups==4, "response"]),
                         "OKP_age4" = sum(d[d$variable=="out_kennel_people"&d$age_groups==4, "response"])/
                           length(d[d$variable=="out_kennel_people"&d$age_groups==4, "response"]),
                         "RF_age4"= sum(d[d$variable=="res_food"&d$age_groups==4, "response"])/
                           length(d[d$variable=="res_food"&d$age_groups==4, "response"]),
                         "RT_age4"= sum(d[d$variable=="res_toys"&d$age_groups==4, "response"])/
                           length(d[d$variable=="res_toys"&d$age_groups==4, "response"]),
                         "FPL_age4"= sum(d[d$variable=="soc_known"&d$age_groups==4, "response"])/
                           length(d[d$variable=="soc_known"&d$age_groups==4, "response"]),
                         "NPL_age4"= sum(d[d$variable=="soc_unknown"&d$age_groups==4, "response"])/
                           length(d[d$variable=="soc_unknown"&d$age_groups==4, "response"])
)

rawProbsSex = data.frame("HND_male" = sum(d[d$variable=="handling"&d$sex==1, "response"])/
                           length(d[d$variable=="handling"&d$sex==1, "response"]),
                         "KP_male" = sum(d[d$variable=="kennel_people"&d$sex==1, "response"])/
                           length(d[d$variable=="kennel_people"&d$sex==1, "response"]),
                         "OKP_male" = sum(d[d$variable=="out_kennel_people"&d$sex==1, "response"])/
                           length(d[d$variable=="out_kennel_people"&d$sex==1, "response"]),
                         "RF_male"= sum(d[d$variable=="res_food"&d$sex==1, "response"])/
                           length(d[d$variable=="res_food"&d$sex==1, "response"]),
                         "RT_male"= sum(d[d$variable=="res_toys"&d$sex==1, "response"])/
                           length(d[d$variable=="res_toys"&d$sex==1, "response"]),
                         "FPL_male"= sum(d[d$variable=="soc_known"&d$sex==1, "response"])/
                           length(d[d$variable=="soc_known"&d$sex==1, "response"]),
                         "NPL_male"= sum(d[d$variable=="soc_unknown"&d$sex==1, "response"])/
                           length(d[d$variable=="soc_unknown"&d$sex==1, "response"]),
                         "HND_female" = sum(d[d$variable=="handling"&d$sex==2, "response"])/
                           length(d[d$variable=="handling"&d$sex==2, "response"]),
                         "KP_female" = sum(d[d$variable=="kennel_people"&d$sex==2, "response"])/
                           length(d[d$variable=="kennel_people"&d$sex==2, "response"]),
                         "OKP_female" = sum(d[d$variable=="out_kennel_people"&d$sex==2, "response"])/
                           length(d[d$variable=="out_kennel_people"&d$sex==2, "response"]),
                         "RF_female"= sum(d[d$variable=="res_food"&d$sex==2, "response"])/
                           length(d[d$variable=="res_food"&d$sex==2, "response"]),
                         "RT_female"= sum(d[d$variable=="res_toys"&d$sex==2, "response"])/
                           length(d[d$variable=="res_toys"&d$sex==2, "response"]),
                         "FPL_female"= sum(d[d$variable=="soc_known"&d$sex==2, "response"])/
                           length(d[d$variable=="soc_known"&d$sex==2, "response"]),
                         "NPL_female"= sum(d[d$variable=="soc_unknown"&d$sex==2, "response"])/
                           length(d[d$variable=="soc_unknown"&d$sex==2, "response"])
)

rawProbsAge <- apply(rawProbsAge, 2, function(x) as.numeric(x))
rawProbsSex <- apply(rawProbsSex, 2, function(x) as.numeric(x))

plotDF_sex <- data.frame( parameter = rep(c("HND","KP","OKP","EAT","TOY","FPL","UPL"), 2),
                          sex =  rep(c("Male", "Female"), each = 7, length=14) ,
                          mu = apply(predictedProbs[,51:64], 2 , mean ),
                          HDIlow = apply(predictedProbs[,51:64], 2, function(x) HPDI(x, 0.95) )[1,],
                          HDIhigh = apply(predictedProbs[,51:64], 2, function(x) HPDI(x, 0.95) )[2,],
                          rawProbs = rawProbsSex,
                          meanMu = rep(apply(predictedProbs[,12:13], 2, mean), each=7),
                          meanHDIlow = rep(apply(predictedProbs[,12:13], 2, function(z) HPDI(z, 0.95))[,1], each=7),
                          meanHDIhigh = rep(apply(predictedProbs[,12:13], 2, function(z) HPDI(z, 0.95))[,2], each=7)
)

plotDF_sex$sex <- factor(plotDF_sex$sex, levels = c("Female","Male"))
plotDF_sex$parameter <- factor(plotDF_sex$parameter, levels=c("HND","KP","OKP","FPL","UPL","EAT","TOY"))

p1 <- ggplot( plotDF_sex, aes(parameter, mu, label=parameter)) +
      geom_text(nudge_y = -0.02, size=5) +
      geom_point(size=2, position=position_dodge(0.5)) +
      geom_errorbar( aes(ymin = HDIlow, ymax = HDIhigh), lwd = 0.6, width=0, position=position_dodge(0.5)) +
      geom_point( aes(parameter, rawProbs), col = rangi2, size = 2,shape=17, position=position_dodge(0.5)) +
      theme_bw(base_family = "Times") +
      facet_wrap(~sex, nrow = 1) +
      ylab("Probability\n") +
      xlab("") +
      ggtitle("A") +
      theme( plot.title = element_text(size=20, hjust = 0),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.text = element_text(size=25, colour = "black"),
             axis.text.y = element_text(size=20, colour = "black"),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.title.y = element_text(size=20, colour = "black"))

#ggsave( "Fig_1a.tiff", last_plot(), width=9, height=5)


plotDF_age <- data.frame( parameter = rep(c("HND","KP","OKP","EAT","TOY","FPL","UPL"), 4),
                          ageGroup =  rep(c("4 to 10 months", "10 months to 3 years", "3 to 6 years", "Over 6 years"),
                                          each = 7, length=28) ,
                          mu = apply(predictedProbs[,23:50], 2 , mean ),
                          HDIlow = apply(predictedProbs[,23:50], 2, function(x) HPDI(x, 0.95) )[1,],
                          HDIhigh = apply(predictedProbs[,23:50], 2, function(x) HPDI(x, 0.95) )[2,],
                          rawProbs = rawProbsAge,
                          meanMu = rep(apply(predictedProbs[,8:11], 2, mean), each=7),
                          meanHDIlow = rep(apply(predictedProbs[,8:11], 2, function(z) HPDI(z, 0.95))[,1], each=7),
                          meanHDIhigh = rep(apply(predictedProbs[,8:11], 2, function(z) HPDI(z, 0.95))[,2], each=7)
)

plotDF_age$ageGroup <- factor(plotDF_age$ageGroup, levels = unique(plotDF_age$ageGroup))
plotDF_age$parameter <- factor(plotDF_age$parameter, levels=c("HND","KP","OKP","FPL","UPL","EAT","TOY"))

p2 <- ggplot( plotDF_age, aes(parameter, mu, label=parameter)) +
      geom_text(nudge_y = -0.02, size=5) +
      geom_point(size=2, position=position_dodge(0.5)) +
      geom_errorbar( aes(ymin = HDIlow, ymax = HDIhigh), lwd = 0.6, width=0, position=position_dodge(0.5)) +
      geom_point( aes(parameter, rawProbs), col = rangi2, size = 2,shape=17, position=position_dodge(0.5)) +
      theme_bw(base_family = "Times") +
      facet_wrap(~ageGroup, nrow = 1) +
      ylab("Probability\n") +
      xlab("") +
      ggtitle("B") +
      theme( plot.title = element_text(size=20, hjust = 0),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.text = element_text(size=25, colour = "black"),
             axis.text.y = element_text(size=20, colour = "black"),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.title.y = element_text(size=20, colour = "black"))

#ggsave( "Fig_1b.tiff", last_plot(), width=14, height=6)

g <- arrangeGrob( p1, p2, nrow=2 )
ggsave("Fig1.tiff", g, width = 15, height = 10, units = "in", dpi = 500)

#===========================================================================
# Now compute the odds ratios for all different comparisons and interactions (in Supporting Material tables)
#===========================================================================

# convert the predictedLogOdds to odds
odds <- as.data.frame(apply(predictedLogOdds, 2 , exp))

comparisons <- with(odds,
                    data.frame("Weight" = exp(fixedEffects$weight), "totalDays" = exp(fixedEffects$totalDays),
                               "BDCH.v.BBH" = BDCH / BBH, "BDCH.v.BOW" = BDCH/BOW, "BOW.v.BBH"=BOW/BBH,
                               "NN.v.NOS" = notNeutered/neuteredOnSite, "NN.v.N" =notNeutered/neutered, "N.v.NOS"=neutered / neuteredOnSite,
                               "Gift.v.Return"=gift/return, "Gift.v.Stray"=gift/stray, "Stray.v.return"=stray/return, "females.v.male" = female / male,
                               # contexts
                               "HND.v.KP" = HND/KP,"HND.v.OKP" = HND/OKP,"HND.v.RF"=HND/RF, "HND.v.RT"=HND/RT,"HND.V.FPL"=HND/FPL,"HND.V.NPL"=HND/NPL,
                               "KP.V.OKP" = KP/OKP,"KP.V.RF"=HND/RF,"KP.V.RT"=KP/RT,"KP.V.FPL"=KP/FPL, "KP.V.NPL"=KP/NPL,
                               "OKP.V.RF"=OKP/RF, "OKP.V.RT"=OKP/RT, "OKP.V.FPL"=OKP/FPL, "OKP.V.NPL"=OKP/NPL,
                               "RF.V.RT"=RF/RT, "RF.V.FPL"=RF/FPL, "RF.V.NPL"=RF/NPL,
                               "RT.V.FPL"=RT/FPL, "RT.V.NPL"=RT/NPL, "FPL.V.NPL"=FPL/NPL,
                               # ages
                               "age1.v.age2"=age1/age2,"age1.v.age3"=age1/age3,"age1.v.age4"=age1/age4,"age2.v.age3"=age2/age3,
                               "age2.v.age4"=age2/age4,"age3.v.age4"=age3/age4,
                               # two way interactions: contexts and sex
                               "HND.V.KP_female.v.male" = (HND_female/KP_female) / (HND_male/KP_male),
                               "HND.V.OKP_female.v.male" = (HND_female/OKP_female) / (HND_male/OKP_male),
                               "HND.V.RF_female.v.male" = (HND_female/RF_female) / (HND_male/RF_male),
                               "HND.V.RT_female.v.male" = (HND_female/RT_female) / (HND_male/RT_male),
                               "HND.V.FPL_female.v.male" = (HND_female/FPL_female) / (HND_male/FPL_male),
                               "HND.V.NPL_female.v.male" = (HND_female/NPL_female) / (HND_male/NPL_male),
                               "KP.V.OKP_female.v.male" = (KP_female/OKP_female) / (KP_male/OKP_male),
                               "KP.V.RF_female.v.male" = (KP_female/RF_female) / (KP_male/RF_male),
                               "KP.V.RT_female.v.male" = (KP_female/RT_female) / (KP_male/RT_male),
                               "KP.V.FPL_female.v.male" = (KP_female/FPL_female) / (KP_male/FPL_male),
                               "KP.V.NPL_female.v.male" = (KP_female/NPL_female) / (KP_male/NPL_male),
                               "OKP.V.RF_female.v.male" = (OKP_female/RF_female) / (OKP_male/RF_male),
                               "OKP.V.RT_female.v.male" = (OKP_female/RT_female) / (OKP_male/RT_male),
                               "OKP.V.FPL_female.v.male" = (OKP_female/FPL_female) / (OKP_male/FPL_male),
                               "OKP.V.NPL_female.v.male" = (OKP_female/NPL_female) / (OKP_male/NPL_male),
                               "RF.V.RT_female.v.male" = (RF_female/RT_female) / (RF_male/RT_male),
                               "RF.V.FPL_female.v.male" = (RF_female/FPL_female) / (RF_male/FPL_male),
                               "RF.V.NPL_female.v.male" = (RF_female/NPL_female) / (RF_male/NPL_male),
                               "RT.V.FPL_female.v.male" = (RT_female/FPL_female) / (RT_male/FPL_male),
                               "RT.V.NPL_female.v.male" = (RT_female/NPL_female) / (RT_male/NPL_male),
                               "FPL.V.NPL_female.v.male" = (FPL_female/NPL_female) / (FPL_male/NPL_male),

                               # two way interactions: contexts and age
                               "HND.V.KP_age1.v.age2" = (HND_age1/KP_age1) / (HND_age2/KP_age2),
                               "HND.V.OKP_age1.v.age2" = (HND_age1/OKP_age1) / (HND_age2/OKP_age2),
                               "HND.V.RF_age1.v.age2" = (HND_age1/RF_age1) / (HND_age2/RF_age2),
                               "HND.V.RT_age1.v.age2" = (HND_age1/RT_age1) / (HND_age2/RT_age2),
                               "HND.V.FPL_age1.v.age2" = (HND_age1/FPL_age1) / (HND_age2/FPL_age2),
                               "HND.V.NPL_age1.v.age2" = (HND_age1/NPL_age1) / (HND_age2/NPL_age2),
                               "KP.V.OKP_age1.v.age2" = (KP_age1/OKP_age1) / (KP_age2/OKP_age2),
                               "KP.V.RF_age1.v.age2" = (KP_age1/RF_age1) / (KP_age2/RF_age2),
                               "KP.V.RT_age1.v.age2" = (KP_age1/RT_age1) / (KP_age2/RT_age2),
                               "KP.V.FPL_age1.v.age2" = (KP_age1/FPL_age1) / (KP_age2/FPL_age2),
                               "KP.V.NPL_age1.v.age2" = (KP_age1/NPL_age1) / (KP_age2/NPL_age2),
                               "OKP.V.RF_age1.v.age2" = (OKP_age1/RF_age1) / (OKP_age2/RF_age2),
                               "OKP.V.RT_age1.v.age2" = (OKP_age1/RT_age1) / (OKP_age2/RT_age2),
                               "OKP.V.FPL_age1.v.age2" = (OKP_age1/FPL_age1) / (OKP_age2/FPL_age2),
                               "OKP.V.NPL_age1.v.age2" = (OKP_age1/NPL_age1) / (OKP_age2/NPL_age2),
                               "RF.V.RT_age1.v.age2" = (RF_age1/RT_age1) / (RF_age2/RT_age2),
                               "RF.V.FPL_age1.v.age2" = (RF_age1/FPL_age1) / (RF_age2/FPL_age2),
                               "RF.V.NPL_age1.v.age2" = (RF_age1/NPL_age1) / (RF_age2/NPL_age2),
                               "RT.V.FPL_age1.v.age2" = (RT_age1/FPL_age1) / (RT_age2/FPL_age2),
                               "RT.V.NPL_age1.v.age2" = (RT_age1/NPL_age1) / (RT_age2/NPL_age2),
                               "FPL.V.NPL_age1.v.age2" = (FPL_age1/NPL_age1) / (FPL_age2/NPL_age2),

                               "HND.V.KP_age1.v.age3" = (HND_age1/KP_age1) / (HND_age3/KP_age3),
                               "HND.V.OKP_age1.v.age3" = (HND_age1/OKP_age1) / (HND_age3/OKP_age3),
                               "HND.V.RF_age1.v.age3" = (HND_age1/RF_age1) / (HND_age3/RF_age3),
                               "HND.V.RT_age1.v.age3" = (HND_age1/RT_age1) / (HND_age3/RT_age3),
                               "HND.V.FPL_age1.v.age3" = (HND_age1/FPL_age1) / (HND_age3/FPL_age3),
                               "HND.V.NPL_age1.v.age3" = (HND_age1/NPL_age1) / (HND_age3/NPL_age3),
                               "KP.V.OKP_age1.v.age3" = (KP_age1/OKP_age1) / (KP_age3/OKP_age3),
                               "KP.V.RF_age1.v.age3" = (KP_age1/RF_age1) / (KP_age3/RF_age3),
                               "KP.V.RT_age1.v.age3" = (KP_age1/RT_age1) / (KP_age3/RT_age3),
                               "KP.V.FPL_age1.v.age3" = (KP_age1/FPL_age1) / (KP_age3/FPL_age3),
                               "KP.V.NPL_age1.v.age3" = (KP_age1/NPL_age1) / (KP_age3/NPL_age3),
                               "OKP.V.RF_age1.v.age3" = (OKP_age1/RF_age1) / (OKP_age3/RF_age3),
                               "OKP.V.RT_age1.v.age3" = (OKP_age1/RT_age1) / (OKP_age3/RT_age3),
                               "OKP.V.FPL_age1.v.age3" = (OKP_age1/FPL_age1) / (OKP_age3/FPL_age3),
                               "OKP.V.NPL_age1.v.age3" = (OKP_age1/NPL_age1) / (OKP_age3/NPL_age3),
                               "RF.V.RT_age1.v.age3" = (RF_age1/RT_age1) / (RF_age3/RT_age3),
                               "RF.V.FPL_age1.v.age3" = (RF_age1/FPL_age1) / (RF_age3/FPL_age3),
                               "RF.V.NPL_age1.v.age3" = (RF_age1/NPL_age1) / (RF_age3/NPL_age3),
                               "RT.V.FPL_age1.v.age3" = (RT_age1/FPL_age1) / (RT_age3/FPL_age3),
                               "RT.V.NPL_age1.v.age3" = (RT_age1/NPL_age1) / (RT_age3/NPL_age3),
                               "FPL.V.NPL_age1.v.age3" = (FPL_age1/NPL_age1) / (FPL_age3/NPL_age3),

                               "HND.V.KP_age1.v.age4" = (HND_age1/KP_age1) / (HND_age4/KP_age4),
                               "HND.V.OKP_age1.v.age4" = (HND_age1/OKP_age1) / (HND_age4/OKP_age4),
                               "HND.V.RF_age1.v.age4" = (HND_age1/RF_age1) / (HND_age4/RF_age4),
                               "HND.V.RT_age1.v.age4" = (HND_age1/RT_age1) / (HND_age4/RT_age4),
                               "HND.V.FPL_age1.v.age4" = (HND_age1/FPL_age1) / (HND_age4/FPL_age4),
                               "HND.V.NPL_age1.v.age4" = (HND_age1/NPL_age1) / (HND_age4/NPL_age4),
                               "KP.V.OKP_age1.v.age4" = (KP_age1/OKP_age1) / (KP_age4/OKP_age4),
                               "KP.V.RF_age1.v.age4" = (KP_age1/RF_age1) / (KP_age4/RF_age4),
                               "KP.V.RT_age1.v.age4" = (KP_age1/RT_age1) / (KP_age4/RT_age4),
                               "KP.V.FPL_age1.v.age4" = (KP_age1/FPL_age1) / (KP_age4/FPL_age4),
                               "KP.V.NPL_age1.v.age4" = (KP_age1/NPL_age1) / (KP_age4/NPL_age4),
                               "OKP.V.RF_age1.v.age4" = (OKP_age1/RF_age1) / (OKP_age4/RF_age4),
                               "OKP.V.RT_age1.v.age4" = (OKP_age1/RT_age1) / (OKP_age4/RT_age4),
                               "OKP.V.FPL_age1.v.age4" = (OKP_age1/FPL_age1) / (OKP_age4/FPL_age4),
                               "OKP.V.NPL_age1.v.age4" = (OKP_age1/NPL_age1) / (OKP_age4/NPL_age4),
                               "RF.V.RT_age1.v.age4" = (RF_age1/RT_age1) / (RF_age4/RT_age4),
                               "RF.V.FPL_age1.v.age4" = (RF_age1/FPL_age1) / (RF_age4/FPL_age4),
                               "RF.V.NPL_age1.v.age4" = (RF_age1/NPL_age1) / (RF_age4/NPL_age4),
                               "RT.V.FPL_age1.v.age4" = (RT_age1/FPL_age1) / (RT_age4/FPL_age4),
                               "RT.V.NPL_age1.v.age4" = (RT_age1/NPL_age1) / (RT_age4/NPL_age4),
                               "FPL.V.NPL_age1.v.age4" = (FPL_age1/NPL_age1) / (FPL_age4/NPL_age4),

                               "HND.V.KP_age2.v.age3" = (HND_age2/KP_age2) / (HND_age3/KP_age3),
                               "HND.V.OKP_age2.v.age3" = (HND_age2/OKP_age2) / (HND_age3/OKP_age3),
                               "HND.V.RF_age2.v.age3" = (HND_age2/RF_age2) / (HND_age3/RF_age3),
                               "HND.V.RT_age2.v.age3" = (HND_age2/RT_age2) / (HND_age3/RT_age3),
                               "HND.V.FPL_age2.v.age3" = (HND_age2/FPL_age2) / (HND_age3/FPL_age3),
                               "HND.V.NPL_age2.v.age3" = (HND_age2/NPL_age2) / (HND_age3/NPL_age3),
                               "KP.V.OKP_age2.v.age3" = (KP_age2/OKP_age2) / (KP_age3/OKP_age3),
                               "KP.V.RF_age2.v.age3" = (KP_age2/RF_age2) / (KP_age3/RF_age3),
                               "KP.V.RT_age2.v.age3" = (KP_age2/RT_age2) / (KP_age3/RT_age3),
                               "KP.V.FPL_age2.v.age3" = (KP_age2/FPL_age2) / (KP_age3/FPL_age3),
                               "KP.V.NPL_age2.v.age3" = (KP_age2/NPL_age2) / (KP_age3/NPL_age3),
                               "OKP.V.RF_age2.v.age3" = (OKP_age2/RF_age2) / (OKP_age3/RF_age3),
                               "OKP.V.RT_age2.v.age3" = (OKP_age2/RT_age2) / (OKP_age3/RT_age3),
                               "OKP.V.FPL_age2.v.age3" = (OKP_age2/FPL_age2) / (OKP_age3/FPL_age3),
                               "OKP.V.NPL_age2.v.age3" = (OKP_age2/NPL_age2) / (OKP_age3/NPL_age3),
                               "RF.V.RT_age2.v.age3" = (RF_age2/RT_age2) / (RF_age3/RT_age3),
                               "RF.V.FPL_age2.v.age3" = (RF_age2/FPL_age2) / (RF_age3/FPL_age3),
                               "RF.V.NPL_age2.v.age3" = (RF_age2/NPL_age2) / (RF_age3/NPL_age3),
                               "RT.V.FPL_age2.v.age3" = (RT_age2/FPL_age2) / (RT_age3/FPL_age3),
                               "RT.V.NPL_age2.v.age3" = (RT_age2/NPL_age2) / (RT_age3/NPL_age3),
                               "FPL.V.NPL_age2.v.age3" = (FPL_age2/NPL_age2) / (FPL_age3/NPL_age3),

                               "HND.V.KP_age2.v.age4" = (HND_age2/KP_age2) / (HND_age4/KP_age4),
                               "HND.V.OKP_age2.v.age4" = (HND_age2/OKP_age2) / (HND_age4/OKP_age4),
                               "HND.V.RF_age2.v.age4" = (HND_age2/RF_age2) / (HND_age4/RF_age4),
                               "HND.V.RT_age2.v.age4" = (HND_age2/RT_age2) / (HND_age4/RT_age4),
                               "HND.V.FPL_age2.v.age4" = (HND_age2/FPL_age2) / (HND_age4/FPL_age4),
                               "HND.V.NPL_age2.v.age4" = (HND_age2/NPL_age2) / (HND_age4/NPL_age4),
                               "KP.V.OKP_age2.v.age4" = (KP_age2/OKP_age2) / (KP_age4/OKP_age4),
                               "KP.V.RF_age2.v.age4" = (KP_age2/RF_age2) / (KP_age4/RF_age4),
                               "KP.V.RT_age2.v.age4" = (KP_age2/RT_age2) / (KP_age4/RT_age4),
                               "KP.V.FPL_age2.v.age4" = (KP_age2/FPL_age2) / (KP_age4/FPL_age4),
                               "KP.V.NPL_age2.v.age4" = (KP_age2/NPL_age2) / (KP_age4/NPL_age4),
                               "OKP.V.RF_age2.v.age4" = (OKP_age2/RF_age2) / (OKP_age4/RF_age4),
                               "OKP.V.RT_age2.v.age4" = (OKP_age2/RT_age2) / (OKP_age4/RT_age4),
                               "OKP.V.FPL_age2.v.age4" = (OKP_age2/FPL_age2) / (OKP_age4/FPL_age4),
                               "OKP.V.NPL_age2.v.age4" = (OKP_age2/NPL_age2) / (OKP_age4/NPL_age4),
                               "RF.V.RT_age2.v.age4" = (RF_age2/RT_age2) / (RF_age4/RT_age4),
                               "RF.V.FPL_age2.v.age4" = (RF_age2/FPL_age2) / (RF_age4/FPL_age4),
                               "RF.V.NPL_age2.v.age4" = (RF_age2/NPL_age2) / (RF_age4/NPL_age4),
                               "RT.V.FPL_age2.v.age4" = (RT_age2/FPL_age2) / (RT_age4/FPL_age4),
                               "RT.V.NPL_age2.v.age4" = (RT_age2/NPL_age2) / (RT_age4/NPL_age4),
                               "FPL.V.NPL_age2.v.age4" = (FPL_age2/NPL_age2) / (FPL_age4/NPL_age4),

                               "HND.V.KP_age3.v.age4" = (HND_age3/KP_age3) / (HND_age4/KP_age4),
                               "HND.V.OKP_age3.v.age4" = (HND_age3/OKP_age3) / (HND_age4/OKP_age4),
                               "HND.V.RF_age3.v.age4" = (HND_age3/RF_age3) / (HND_age4/RF_age4),
                               "HND.V.RT_age3.v.age4" = (HND_age3/RT_age3) / (HND_age4/RT_age4),
                               "HND.V.FPL_age3.v.age4" = (HND_age3/FPL_age3) / (HND_age4/FPL_age4),
                               "HND.V.NPL_age3.v.age4" = (HND_age3/NPL_age3) / (HND_age4/NPL_age4),
                               "KP.V.OKP_age3.v.age4" = (KP_age3/OKP_age3) / (KP_age4/OKP_age4),
                               "KP.V.RF_age3.v.age4" = (KP_age3/RF_age3) / (KP_age4/RF_age4),
                               "KP.V.RT_age3.v.age4" = (KP_age3/RT_age3) / (KP_age4/RT_age4),
                               "KP.V.FPL_age3.v.age4" = (KP_age3/FPL_age3) / (KP_age4/FPL_age4),
                               "KP.V.NPL_age3.v.age4" = (KP_age3/NPL_age3) / (KP_age4/NPL_age4),
                               "OKP.V.RF_age3.v.age4" = (OKP_age3/RF_age3) / (OKP_age4/RF_age4),
                               "OKP.V.RT_age3.v.age4" = (OKP_age3/RT_age3) / (OKP_age4/RT_age4),
                               "OKP.V.FPL_age3.v.age4" = (OKP_age3/FPL_age3) / (OKP_age4/FPL_age4),
                               "OKP.V.NPL_age3.v.age4" = (OKP_age3/NPL_age3) / (OKP_age4/NPL_age4),
                               "RF.V.RT_age3.v.age4" = (RF_age3/RT_age3) / (RF_age4/RT_age4),
                               "RF.V.FPL_age3.v.age4" = (RF_age3/FPL_age3) / (RF_age4/FPL_age4),
                               "RF.V.NPL_age3.v.age4" = (RF_age3/NPL_age3) / (RF_age4/NPL_age4),
                               "RT.V.FPL_age3.v.age4" = (RT_age3/FPL_age3) / (RT_age4/FPL_age4),
                               "RT.V.NPL_age3.v.age4" = (RT_age3/NPL_age3) / (RT_age4/NPL_age4),
                               "FPL.V.NPL_age3.v.age4" = (FPL_age3/NPL_age3) / (FPL_age4/NPL_age4)
                    ))

summCompars <- data.frame(paramComp = colnames(comparisons),
                          muDiff = apply(comparisons, 2 , mean),
                          HDIlow = apply(comparisons, 2, function(x) HPDI(x,prob=0.95))[1,] ,
                          HDIhigh = apply(comparisons, 2, function(x) HPDI(x, prob=0.95))[2,])
rownames(summCompars) <- 1:ncol(comparisons)
summCompars$cred <- ifelse(summCompars$HDIhigh < (1/1.25) | summCompars$HDIlow > 1.25 , "YES",
                           ifelse(summCompars$HDIlow >= (1/1.25) & summCompars$HDIhigh <= 1.25 , "NULL", "ROPE"))
summCompars$cred[1:12] <- "No inferential decision" # these effects are not interpreted inferentially
summCompars[,2:4] <- apply(summCompars[,2:4], 2, function(z) round(z,3))

write.csv(summCompars, "MultipleComparisons_LV1data.csv", row.names=F)

#=================================================================================
# Do the same for the LV2 data
#=================================================================================

# load MCMC matrices, get log-likelihoods and get WAIC estimates
fullMI <- as.data.frame(fread("CombinedMCMC_LV2_fullMI.csv"))
ageMI <- as.data.frame(fread("CombinedMCMC_LV2_ageMI.csv"))
sexMI <- as.data.frame(fread("CombinedMCMC_LV2_sexMI.csv"))
noMI <- as.data.frame(fread("CombinedMCMC_LV2_noMI.csv"))

fullMI_waic <- waic(as.matrix(fullMI[, grep("log_lik",colnames(fullMI))]))
ageMI_waic <- waic(as.matrix(ageMI[, grep("log_lik",colnames(ageMI))]))
sexMI_waic <- waic(as.matrix(sexMI[, grep("log_lik",colnames(sexMI))]))
noMI_waic <- waic(as.matrix(noMI[, grep("log_lik",colnames(noMI))]))

WAICtable_LV2 <- data.frame( model = c("Model 1", "Model 2", "Model 3", "Model 4"),
                             WAIC = c(fullMI_waic$waic,ageMI_waic$waic,sexMI_waic$waic,noMI_waic$waic),
                             SE = c(fullMI_waic$se_waic,ageMI_waic$se_waic,sexMI_waic$se_waic,noMI_waic$se_waic)
                            )

# Now, plot the fixed sex and age relationships as in Figure 2
mcmcMat <- as.data.frame(fread("CombinedMCMC_LV2_fullMI.csv"))
fixedEffects <- mcmcMat[ , c(1, grep("Beta", colnames(mcmcMat)))]
colnames(fixedEffects) <- c( "intercept", "DF", "DM", "KD",
                             "age1", "age2", "age3", "male", "totalDays", "weight",
                             "notNeutered", "neuteredOnSite", "gift", "return", "BDCH", "BBH",
                             "DF_age1", "DM_age1", "KD_age1",
                             "DF_age2", "DM_age2", "KD_age2",
                             "DF_age3", "DM_age3", "KD_age3",
                             "DF_male", "DM_male", "KD_male")

# convert dummy coded variables (sum coding) so that all estimates for fixed effects are available
fixedEffects$OKD <- with(fixedEffects, (DF + DM + KD)*-1)
fixedEffects$age4 <- with(fixedEffects, (age1+age2+age3)*-1)
fixedEffects$female <- with(fixedEffects, male*-1)
fixedEffects$neutered <- with(fixedEffects, (notNeutered+neuteredOnSite)*-1)
fixedEffects$stray <- with(fixedEffects, (gift+return)*-1)
fixedEffects$BOW <- with(fixedEffects, (BDCH+BBH)*-1)
fixedEffects$OKD_age1 <- with(fixedEffects, (DM_age1 + DF_age1 + KD_age1)*-1)
fixedEffects$OKD_age2 <- with(fixedEffects, (DM_age2 + DF_age2 + KD_age2)*-1)
fixedEffects$OKD_age3 <- with(fixedEffects, (DM_age3 + DF_age3 + KD_age3)*-1)
fixedEffects$OKD_male <- with(fixedEffects, (DF_male + DM_male + KD_male)*-1)
fixedEffects$DF_age4 <- with(fixedEffects, (DF_age1+DF_age2+DF_age3)*-1)
fixedEffects$DM_age4 <- with(fixedEffects, (DM_age1+DM_age2+DM_age3)*-1)
fixedEffects$KD_age4 <- with(fixedEffects, (KD_age1+KD_age2+KD_age3)*-1)
fixedEffects$OKD_age4 <- with(fixedEffects, (OKD_age1+OKD_age2+OKD_age3)*-1)
fixedEffects$DF_female <- with(fixedEffects, DF_male*-1)
fixedEffects$DM_female <- with(fixedEffects, DM_male*-1)
fixedEffects$KD_female <- with(fixedEffects, KD_male*-1)
fixedEffects$OKD_female <- with(fixedEffects, OKD_male*-1)

# get the predicted log odds of aggression
predictedLogOdds <- with(fixedEffects,
                         data.frame( "DF" = intercept + DF ,"DM" = intercept + DM, "KD" = intercept + KD, "OKD" = intercept + OKD,
                                     "age1" = intercept + age1,
                                     "age2" = intercept + age2,"age3" = intercept + age3, "age4" = intercept + age4,
                                     "male" = intercept + male, "female" = intercept + female, "notNeutered" =intercept+notNeutered,
                                     "neuteredOnSite"=intercept+neuteredOnSite, "neutered"=intercept+ neutered,
                                     "gift" = intercept+gift, "return"=intercept+return, "stray"=intercept+stray,
                                     "BDCH"=intercept+BDCH,"BBH"=intercept+BBH,"BOW"=intercept+BOW,
                                     "DF_age1"=intercept+DF+age1+DF_age1, "DM_age1"=intercept+DM+age1+DM_age1,
                                     "KD_age1"=intercept+KD+age1+KD_age1,"OKD_age1"=intercept+OKD+age1+OKD_age1,
                                     "DF_age2"=intercept+DF+age2+DF_age2, "DM_age2"=intercept+DM+age2+DM_age2,
                                     "KD_age2"=intercept+KD+age2+KD_age2,"OKD_age2"=intercept+OKD+age2+OKD_age2,
                                     "DF_age3"=intercept+DF+age3+DF_age3, "DM_age3"=intercept+DM+age3+DM_age3,
                                     "KD_age3"=intercept+KD+age3+KD_age3,"OKD_age3"=intercept+OKD+age3+OKD_age3,
                                     "DF_age4"=intercept+DF+age4+DF_age4, "DM_age4"=intercept+DM+age4+DM_age4,
                                     "KD_age4"=intercept+KD+age4+KD_age4,"OKD_age4"=intercept+OKD+age4+OKD_age4,
                                     "DF_male"=intercept+DF+male+DF_male, "DM_male"=intercept+DM+male+DM_male,
                                     "KD_male" =intercept+KD+male+KD_male, "OKD_male"=intercept+OKD+male+OKD_male,
                                     "DF_female"=intercept+DF+female+DF_female, "DM_female"=intercept+DM+female+DM_female,
                                     "KD_female" =intercept+KD+female+KD_female, "OKD_female"=intercept+OKD+female+OKD_female
                         ))


# marginalise over the random effects and convert to probabilities for plotting
c2 <- ((16*sqrt(3))/(15*pi))^2
predictedLogOdds_approx <- apply(predictedLogOdds, 2 , function(x) x/sqrt(1 + c2 * mcmcMat[,"sigmaID"]^2))
predictedProbs <- apply(predictedLogOdds_approx, 2, function(x) logistic(x) )

d <- LV2_d

# get raw probabilities of aggression by age groups and sex
rawProbsAge = data.frame("DF_age1" = sum(d[d$variable=="dogs_female"&d$age_groups==1, "response"])/
                           length(d[d$variable=="dogs_female"&d$age_groups==1, "response"]),
                         "DM_age1" = sum(d[d$variable=="dogs_male"&d$age_groups==1, "response"])/
                           length(d[d$variable=="dogs_male"&d$age_groups==1, "response"]),
                         "KD_age1" = sum(d[d$variable=="kennel_dogs"&d$age_groups==1, "response"])/
                           length(d[d$variable=="kennel_dogs"&d$age_groups==1, "response"]),
                         "OKD_age1"= sum(d[d$variable=="out_kennel_dogs"&d$age_groups==1, "response"])/
                           length(d[d$variable=="out_kennel_dogs"&d$age_groups==1, "response"]),
                         "DF_age2" = sum(d[d$variable=="dogs_female"&d$age_groups==2, "response"])/
                           length(d[d$variable=="dogs_female"&d$age_groups==2, "response"]),
                         "DM_age2" = sum(d[d$variable=="dogs_male"&d$age_groups==2, "response"])/
                           length(d[d$variable=="dogs_male"&d$age_groups==2, "response"]),
                         "KD_age2" = sum(d[d$variable=="kennel_dogs"&d$age_groups==2, "response"])/
                           length(d[d$variable=="kennel_dogs"&d$age_groups==2, "response"]),
                         "OKD_age2"= sum(d[d$variable=="out_kennel_dogs"&d$age_groups==2, "response"])/
                           length(d[d$variable=="out_kennel_dogs"&d$age_groups==2, "response"]),
                         "DF_age3" = sum(d[d$variable=="dogs_female"&d$age_groups==3, "response"])/
                           length(d[d$variable=="dogs_female"&d$age_groups==3, "response"]),
                         "DM_age3" = sum(d[d$variable=="dogs_male"&d$age_groups==3, "response"])/
                           length(d[d$variable=="dogs_male"&d$age_groups==3, "response"]),
                         "KD_age3" = sum(d[d$variable=="kennel_dogs"&d$age_groups==3, "response"])/
                           length(d[d$variable=="kennel_dogs"&d$age_groups==3, "response"]),
                         "OKD_age3"= sum(d[d$variable=="out_kennel_dogs"&d$age_groups==3, "response"])/
                           length(d[d$variable=="out_kennel_dogs"&d$age_groups==3, "response"]),
                         "DF_age4" = sum(d[d$variable=="dogs_female"&d$age_groups==4, "response"])/
                           length(d[d$variable=="dogs_female"&d$age_groups==4, "response"]),
                         "DM_age4" = sum(d[d$variable=="dogs_male"&d$age_groups==4, "response"])/
                           length(d[d$variable=="dogs_male"&d$age_groups==4, "response"]),
                         "KD_age4" = sum(d[d$variable=="kennel_dogs"&d$age_groups==4, "response"])/
                           length(d[d$variable=="kennel_dogs"&d$age_groups==4, "response"]),
                         "OKD_age4"= sum(d[d$variable=="out_kennel_dogs"&d$age_groups==4, "response"])/
                           length(d[d$variable=="out_kennel_dogs"&d$age_groups==4, "response"])
)

rawProbsSex = data.frame("DF_male" = sum(d[d$variable=="dogs_female"&d$sex==1, "response"])/
                           length(d[d$variable=="dogs_female"&d$sex==1, "response"]),
                         "DM_male" = sum(d[d$variable=="dogs_male"&d$sex==1, "response"])/
                           length(d[d$variable=="dogs_male"&d$sex==1, "response"]),
                         "KD_male" = sum(d[d$variable=="kennel_dogs"&d$sex==1, "response"])/
                           length(d[d$variable=="kennel_dogs"&d$sex==1, "response"]),
                         "OKD_male"= sum(d[d$variable=="out_kennel_dogs"&d$sex==1, "response"])/
                           length(d[d$variable=="out_kennel_dogs"&d$sex==1, "response"]),
                         "DF_female" = sum(d[d$variable=="dogs_female"&d$sex==2, "response"])/
                           length(d[d$variable=="dogs_female"&d$sex==2, "response"]),
                         "DM_female" = sum(d[d$variable=="dogs_male"&d$sex==2, "response"])/
                           length(d[d$variable=="dogs_male"&d$sex==2, "response"]),
                         "KD_female" = sum(d[d$variable=="kennel_dogs"&d$sex==2, "response"])/
                           length(d[d$variable=="kennel_dogs"&d$sex==2, "response"]),
                         "OKD_female"= sum(d[d$variable=="out_kennel_dogs"&d$sex==2, "response"])/
                           length(d[d$variable=="out_kennel_dogs"&d$sex==2, "response"])
)

rawProbsAge <- apply(rawProbsAge, 2, function(x) as.numeric(x))
rawProbsSex <- apply(rawProbsSex, 2, function(x) as.numeric(x))

plotDF_sex <- data.frame( parameter = rep(c("DF","DM","KD","OKD"), 2),
                          sex =  rep(c("Male", "Female"), each = 4) ,
                          mu = apply(predictedProbs[,36:43], 2 , mean ),
                          HDIlow = apply(predictedProbs[,36:43], 2, function(x) HPDI(x, 0.95) )[1,],
                          HDIhigh = apply(predictedProbs[,36:43], 2, function(x) HPDI(x, 0.95) )[2,],
                          rawProbs = rawProbsSex,
                          meanMu = rep(apply(predictedProbs[,9:10], 2, mean), each=4),
                          meanHDIlow = rep(apply(predictedProbs[,9:10], 2, function(z) HPDI(z, 0.95))[,1], each=4),
                          meanHDIhigh = rep(apply(predictedProbs[,9:10], 2, function(z) HPDI(z, 0.95))[,2], each=4)
)

plotDF_sex$sex <- factor(plotDF_sex$sex, levels = c("Female","Male"))
plotDF_sex$parameter <- factor(plotDF_sex$parameter,levels=c("KD","OKD","DF","DM"))

p1 <- ggplot( plotDF_sex, aes(parameter, mu, label=parameter)) +
      geom_text(nudge_y = -0.02, size=5) +
      geom_point(size=2, position=position_dodge(0.5)) +
      geom_errorbar( aes(ymin = HDIlow, ymax = HDIhigh), lwd = 0.6, width=0, position=position_dodge(0.5)) +
      geom_point( aes(parameter, rawProbs), col = rangi2, size = 2,shape=17, position=position_dodge(0.5)) +
      theme_bw(base_family = "Times") +
      facet_wrap(~sex, nrow = 1) +
      ylab("Probability\n") +
      xlab("") +
      ggtitle("A") +
      theme( plot.title = element_text(size=20, hjust = 0),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.text = element_text(size=25, colour = "black"),
             axis.text.y = element_text(size=20, colour = "black"),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.title.y = element_text(size=20, colour = "black"))

#ggsave( "Fig_2a.tiff", last_plot(), width=9, height=5)

plotDF_age <- data.frame( parameter = rep(c("DF","DM","KD","OKD"), 4),
                          ageGroup =  rep(c("4 to 10 months", "10 months to 3 years", "3 to 6 years", "Over 6 years"),
                                          each = 4, length=16) ,
                          mu = apply(predictedProbs[,20:35], 2 , mean ),
                          HDIlow = apply(predictedProbs[,20:35], 2, function(x) HPDI(x, 0.95) )[1,],
                          HDIhigh = apply(predictedProbs[,20:35], 2, function(x) HPDI(x, 0.95) )[2,],
                          rawProbs = rawProbsAge,
                          meanMu = rep(apply(predictedProbs[,5:8], 2, mean), each=4),
                          meanHDIlow = rep(apply(predictedProbs[,5:8], 2, function(z) HPDI(z, 0.95))[,1], each=4),
                          meanHDIhigh = rep(apply(predictedProbs[,5:8], 2, function(z) HPDI(z, 0.95))[,2], each=4)
)

plotDF_age$ageGroup <- factor(plotDF_age$ageGroup, levels = unique(plotDF_age$ageGroup))
plotDF_age$parameter <- factor(plotDF_age$parameter, levels=c("KD","OKD","DF","DM"))

p2 <- ggplot( plotDF_age, aes(parameter, mu, label=parameter)) +
      geom_text(nudge_y = 0.035, size=5) +
      geom_point(size=2, position=position_dodge(0.5)) +
      geom_errorbar( aes(ymin = HDIlow, ymax = HDIhigh), lwd = 0.6, width=0, position=position_dodge(0.5)) +
      geom_point( aes(parameter, rawProbs), col = rangi2, size = 2,shape=17, position=position_dodge(0.5)) +
      theme_bw(base_family = "Times") +
      facet_wrap(~ageGroup, nrow = 1) +
      ylab("Probability\n") +
      xlab("") +
      ggtitle("B") +
      theme( plot.title = element_text(size=20, hjust = 0),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.text = element_text(size=25, colour = "black"),
             axis.text.y = element_text(size=20, colour = "black"),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.title.y = element_text(size=20, colour = "black"))

#ggsave( "Fig_2b.tiff", last_plot(), width=14, height=6)

g <- arrangeGrob( p1, p2, nrow=2 )
ggsave("Fig2.tiff", g, width = 15, height = 10, units = "in")

# get odds ratio results
odds <- as.data.frame(apply(predictedLogOdds, 2 , exp))

comparisons <- with(odds,
                    data.frame("Weight" = exp(fixedEffects$weight), "totalDays" = exp(fixedEffects$totalDays),
                               "BDCH.v.BBH" = BDCH / BBH, "BDCH.v.BOW" = BDCH/BOW, "BOW.v.BBH"=BOW/BBH,
                               "NN.v.NOS" = notNeutered/neuteredOnSite, "NN.v.N" =notNeutered/neutered, "N.v.NOS"=neutered / neuteredOnSite,
                               "Gift.v.Return"=gift/return, "Gift.v.Stray"=gift/stray, "Stray.v.return"=stray/return, "females.v.male" = female / male,
                               # contexts
                               "DF.v.DM" = DF/DM,"DF.v.KD" = DF/KD,"DF.v.OKD"=DF/OKD, "DM.V.KD" = DM/KD,"DM.V.OKD"=DF/OKD,"KD.V.OKD"=KD/OKD,
                               # ages
                               "age1.v.age2"=age1/age2,"age1.v.age3"=age1/age3,"age1.v.age4"=age1/age4,"age2.v.age3"=age2/age3,
                               "age2.v.age4"=age2/age4,"age3.v.age4"=age3/age4,
                               # two way interactions: contexts and sex
                               "DF.V.DM_female.v.male" = (DF_female/DM_female) / (DF_male/DM_male),
                               "DF.V.KD_female.v.male" = (DF_female/KD_female) / (DF_male/KD_male),
                               "DF.V.OKD_female.v.male" = (DF_female/OKD_female) / (DF_male/OKD_male),
                               "DM.V.KD_female.v.male" = (DM_female/KD_female) / (DM_male/KD_male),
                               "DM.V.OKD_female.v.male" = (DM_female/OKD_female) / (DM_male/OKD_male),
                               "KD.V.OKD_female.v.male" = (KD_female/OKD_female) / (KD_male/OKD_male),

                               # two way interactions: contexts and age
                               "DF.V.DM_age1.v.age2" = (DF_age1/DM_age1) / (DF_age2/DM_age2),
                               "DF.V.KD_age1.v.age2" = (DF_age1/KD_age1) / (DF_age2/KD_age2),
                               "DF.V.OKD_age1.v.age2" = (DF_age1/OKD_age1) / (DF_age2/OKD_age2),
                               "DM.V.KD_age1.v.age2" = (DM_age1/KD_age1) / (DM_age2/KD_age2),
                               "DM.V.OKD_age1.v.age2" = (DM_age1/OKD_age1) / (DM_age2/OKD_age2),
                               "KD.V.OKD_age1.v.age2" = (KD_age1/OKD_age1) / (KD_age2/OKD_age2),

                               "DF.V.DM_age1.v.age3" = (DF_age1/DM_age1) / (DF_age3/DM_age3),
                               "DF.V.KD_age1.v.age3" = (DF_age1/KD_age1) / (DF_age3/KD_age3),
                               "DF.V.OKD_age1.v.age3" = (DF_age1/OKD_age1) / (DF_age3/OKD_age3),
                               "DM.V.KD_age1.v.age3" = (DM_age1/KD_age1) / (DM_age3/KD_age3),
                               "DM.V.OKD_age1.v.age3" = (DM_age1/OKD_age1) / (DM_age3/OKD_age3),
                               "KD.V.OKD_age1.v.age3" = (KD_age1/OKD_age1) / (KD_age3/OKD_age3),

                               "DF.V.DM_age1.v.age4" = (DF_age1/DM_age1) / (DF_age4/DM_age4),
                               "DF.V.KD_age1.v.age4" = (DF_age1/KD_age1) / (DF_age4/KD_age4),
                               "DF.V.OKD_age1.v.age4" = (DF_age1/OKD_age1) / (DF_age4/OKD_age4),
                               "DM.V.KD_age1.v.age4" = (DM_age1/KD_age1) / (DM_age4/KD_age4),
                               "DM.V.OKD_age1.v.age4" = (DM_age1/OKD_age1) / (DM_age4/OKD_age4),
                               "KD.V.OKD_age1.v.age4" = (KD_age1/OKD_age1) / (KD_age4/OKD_age4),

                               "DF.V.DM_age2.v.age3" = (DF_age2/DM_age2) / (DF_age3/DM_age3),
                               "DF.V.KD_age2.v.age3" = (DF_age2/KD_age2) / (DF_age3/KD_age3),
                               "DF.V.OKD_age2.v.age3" = (DF_age2/OKD_age2) / (DF_age3/OKD_age3),
                               "DM.V.KD_age2.v.age3" = (DM_age2/KD_age2) / (DM_age3/KD_age3),
                               "DM.V.OKD_age2.v.age3" = (DM_age2/OKD_age2) / (DM_age3/OKD_age3),
                               "KD.V.OKD_age2.v.age3" = (KD_age2/OKD_age2) / (KD_age3/OKD_age3),

                               "DF.V.DM_age2.v.age4" = (DF_age2/DM_age2) / (DF_age4/DM_age4),
                               "DF.V.KD_age2.v.age4" = (DF_age2/KD_age2) / (DF_age4/KD_age4),
                               "DF.V.OKD_age2.v.age4" = (DF_age2/OKD_age2) / (DF_age4/OKD_age4),
                               "DM.V.KD_age2.v.age4" = (DM_age2/KD_age2) / (DM_age4/KD_age4),
                               "DM.V.OKD_age2.v.age4" = (DM_age2/OKD_age2) / (DM_age4/OKD_age4),
                               "KD.V.OKD_age2.v.age4" = (KD_age2/OKD_age2) / (KD_age4/OKD_age4),

                               "DF.V.DM_age3.v.age4" = (DF_age3/DM_age3) / (DF_age4/DM_age4),
                               "DF.V.KD_age3.v.age4" = (DF_age3/KD_age3) / (DF_age4/KD_age4),
                               "DF.V.OKD_age3.v.age4" = (DF_age3/OKD_age3) / (DF_age4/OKD_age4),
                               "DM.V.KD_age3.v.age4" = (DM_age3/KD_age3) / (DM_age4/KD_age4),
                               "DM.V.OKD_age3.v.age4" = (DM_age3/OKD_age3) / (DM_age4/OKD_age4),
                               "KD.V.OKD_age3.v.age4" = (KD_age3/OKD_age3) / (KD_age4/OKD_age4)
                    ))

summCompars <- data.frame(paramComp = colnames(comparisons),
                          muDiff = apply(comparisons, 2 , mean),
                          HDIlow = apply(comparisons, 2, function(x) HPDI(x,prob=0.95))[1,] ,
                          HDIhigh = apply(comparisons, 2, function(x) HPDI(x, prob=0.95))[2,])
rownames(summCompars) <- 1:ncol(comparisons)
summCompars$cred <- ifelse(summCompars$HDIhigh <= (1/1.25) | summCompars$HDIlow >= 1.25 , "YES",
                           ifelse(summCompars$HDIlow >= (1/1.25) & summCompars$HDIhigh <= 1.25 , "NULL", "ROPE"))
summCompars$cred[1:12] <- "No inferential decision"
summCompars[,2:4] <- apply(summCompars[,2:4], 2, function(z) round(z,3))

write.csv(summCompars, "MultipleComparisons_LV2data.csv", row.names=F)
