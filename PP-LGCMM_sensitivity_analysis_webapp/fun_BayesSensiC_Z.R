
library(mvtnorm)
library(ggplot2)
library(tidyverse)
options(tibble.print_max = Inf, tibble.width = Inf)
library(MplusAutomation)
# 

# 
try <- function(){
  datN200T5=read.csv('~/Library/CloudStorage/Box-Box/Labs/SensitivityAnalysis/BSAgcm/Rfunction/LGCMdata.csv',header = T)
  dat = datN200T5  # data frame containing variables in the LGCM model 
  Ycols = c(1:5)  # columns for repeated measures of outcome from the first to the last time points 
  SY.loadings = c(-4,-3,-2,-1,0)  # loadings of the outcome slope
  Mcols = c(6:10)  # columns for repeated measures of mediator from the first to the last time points
  SM.loadings = c(0,1,2,3,4) # loadings of the mediator slope,
  within_person_residual_cov = "eM1-eM5 (sigma_eM);
 eY1-eY5 (sigma_eY);"
  Xcol = 11 # column for the independent variable 
  Z1cols = c(12) # columns for the level-2 covariates in both the level-2 models of mediator and outcome  
  Z2cols = NULL # columns for the level-2 covariates in the level-2 model of mediator only
  Z3cols = NULL # columns for the level-2 covariates in the level-2 model of outcome only
  alph=0.05 # alpha level per effect
  BITERATIONS=20000
  Bayes_originalmodel= "
IM by M1-M5@1;
SM by M1@0 M2@1 M3@2 M4@3 M5@4;
IY by Y1-Y5@1;
SY by Y1@-4 Y2@-3 Y3@-2 Y4@-1 Y5@0;
[M1-M5@0]; [Y1-Y5@0];
M1-M5 (sigma_eM); Y1-Y5 (sigma_eY);

[IM SM] (a0_IM a0_SM);
[IY SY]  (b0_IY b0_SY);
IM on X  (a1_IM);
SM on X  (a1_SM);
IY on X IM SM Z (b1_IY b2_IY b3_IY b4_IY);
SY on X IM SM Z (b1_SY b2_SY b3_SY b4_SY);
IM (revarIM);SM (revarSM);
IM with SM (IMwSM);
IY (revarIY);SY (revarSY);
IY with SY (IYwSY);
IY SY with IM@0 SM@0;
  "
  priors_originalmodel="
revarIM ~ iw(1,3); revarSM ~iw(1,3);IMwSM ~ iw(0,3);
revarIY ~ iw(1,3); revarSY ~ iw(1,3);IYwSY ~ iw(0,3);"
  priors_muC="
aC_IM ~ N(0.39, 0.25); aC_SM ~ N(0.39, 0.25);
bC_IY ~ N(-0.39, 0.25); bC_SY ~ N(-0.39, 0.25); "
}

fun.BayesC_Z <- function(
  dat = datN200T5, # data frame containing variables in the LGCM model
  Ycols = c(1:5), # columns for repeated measures of outcome from the first to the last time points
  SY.loadings = c(-4,-3,-2,-1,0), # loadings of the outcome slope
  Mcols = c(6:10), # columns for repeated measures of mediator from the first to the last time points
  SM.loadings = c(0,1,2,3,4), # loadings of the mediator slope
  Xcol = 11, # column for the independent variable
  Z1cols = c(12), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
  Z2cols = NULL, # columns for the level-2 covariates in the level-2 model of mediator only
  Z3cols = NULL, # columns for the level-2 covariates in the level-2 model of outcome only
  within_person_residual_cov = "eM1-eM5 (sigma_eM);
 eY1-eY5 (sigma_eY);",
  priors_muC="
aIM_C ~ N(0.39, 0.25); aSM_C ~ N(0.39, 0.25);
bIY_C ~ N(-0.39, 0.25); bSY_C ~ N(-0.39, 0.25);  ",
  BITERATIONS=20000,THIN=5,CHAINS=2,
  
  non_default_priors_original = TRUE,
  Bayes_originalmodel= "
IM by M1-M5@1;
SM by M1@0 M2@1 M3@2 M4@3 M5@4;
IY by Y1-Y5@1;
SY by Y1@-4 Y2@-3 Y3@-2 Y4@-1 Y5@0;
[M1-M5@0]; [Y1-Y5@0]; M1-M5; Y1-Y5;
[IM SM];[IY SY]; IY on X Z; SY on X Z;

IM on X  (aIM_X);
SM on X  (aSM_X);
IY on IM SM (bIY_IM bIY_SM);
SY on IM SM (bSY_IM bSY_SM);

IM (revarIM);SM (revarSM); IM with SM (IMwSM);
IY (revarIY);SY (revarSY); IY with SY (IYwSY);
IY SY with IM@0 SM@0;
  ",
  priors_originalmodel="
revarIM ~ iw(1,3); revarSM ~iw(1,3);IMwSM ~ iw(0,3);
revarIY ~ iw(1,3); revarSY ~ iw(1,3);IYwSY ~ iw(0,3);"
)
{
  
  YTT = length(Ycols); MTT = length(Mcols)
  colnames(dat)[Ycols] = paste0('Y', 1:length(Ycols) )
  colnames(dat)[Mcols] = paste0('M', 1:length(Mcols) )
  colnames(dat)[Xcol] = 'X'
  Z1names=colnames(dat)[Z1cols]
  if( sum(is.na(Z1names))>0 ){Z1names = NULL}
  Z2names=colnames(dat)[Z2cols]
  if( sum(is.na(Z2names))>0 ){Z2names = NULL}
  Z3names=colnames(dat)[Z3cols]
  if( sum(is.na(Z3names))>0 ){Z3names = NULL}
  # model objects for Model M0 and Model M1 using the Bayes estimator in Mplus ----
  
  ## if non-default prior for original model params
  if( non_default_priors_original ) {
    naive.mplusModel.Bayes = mplusObject(
      TITLE = "lgcm",
      ANALYSIS = paste0(
        "MODEL = NOCOVARIANCES;
  ESTIMATOR = BAYES; 
  POINT = MEAN;", "\n",
        "THIN = ", THIN,";\n",
        "CHAINS = " , CHAINS,";\n",
        "PROCESS = 2;
  ALGORITHM = GIBBS; 
  BITERATIONS = ", BITERATIONS,";\n",
        collapse = "" ),
      
      MODEL = Bayes_originalmodel,
      
      MODELPRIORS = paste0( priors_originalmodel ),
      
      MODELCONSTRAINT = paste0(
        # specific
        "NEW (X_IM_IY, X_SM_IY, X_IM_SY, X_SM_SY  );\n", 
        "X_IM_IY = aIM_X*bIY_IM; \n",
        "X_SM_IY = aSM_X*bIY_SM; \n",
        "X_IM_SY = aIM_X*bSY_IM; \n",
        "X_SM_SY = aSM_X*bSY_SM; \n"

      ),
      
      OUTPUT = "TECH1 TECH3 TECH4 TECH8;\n",
      
      usevariables = colnames(dat), autov = F,
      rdata = dat
    )
    
    sensiC.mplusModel.Bayes = mplusObject(
      TITLE = "lgcm",
      ANALYSIS = paste0(
        "MODEL = NOCOVARIANCES;
  ESTIMATOR = BAYES; 
  POINT = MEAN;", "\n",
        "THIN = ", THIN,";\n",
        "CHAINS = " , CHAINS,";\n",
        "PROCESS = 2;
  ALGORITHM = GIBBS; 
  BITERATIONS = ", BITERATIONS,";\n",
        collapse = "" ),
      
      MODEL = paste0(Bayes_originalmodel,
                     "\n",
                     "C by;\n",
                     "[C@0];\n",
                     "C@1;\n",
                     "IM on C (aIM_C);\n",
                     "SM on C (aSM_C);\n",
                     "IY on C (bIY_C);\n",
                     "SY on C (bSY_C);\n"
      ),
      MODELPRIORS = paste0(priors_originalmodel,
                           "\n", priors_muC),
      
      MODELCONSTRAINT = paste0(
        # specific
        "NEW (X_IM_IY, X_SM_IY, X_IM_SY, X_SM_SY  );\n", 
        "X_IM_IY = aIM_X*bIY_IM; \n",
        "X_SM_IY = aSM_X*bIY_SM; \n",
        "X_IM_SY = aIM_X*bSY_IM; \n",
        "X_SM_SY = aSM_X*bSY_SM; \n"
        
      ),
      
      OUTPUT = "TECH1 TECH3 TECH4 TECH8;\n",
      
      usevariables = colnames(dat), autov = F,
      rdata = dat
    )
  }
  
  ## if default Mplus prior 
  if(  !non_default_priors_original ) {
    
    # original model M0 code
    mpluscodeModel_M0= paste0(
      # level 1 model
      "IM by ", paste0("M",1:MTT,"@",1, collapse = "\n"), ";\n",
      "SM by ", paste0("M",1:MTT,"@",SM.loadings, collapse = "\n"), ";\n",
      "IY by ", paste0("Y",1:YTT,"@",1, collapse = "\n"), ";\n",
      "SY by ", paste0("Y",1:YTT,"@",SY.loadings, collapse = "\n"), ";\n",
      
      "[M1-M",MTT,"@0]\n;",
      "[Y1-Y",YTT,"@0]\n;",
      
      "! created within-person centered residual;\n",
      paste0("eM",1:MTT," by M",1:MTT,"@1;\n", collapse = ""),
      paste0("eY",1:YTT," by Y",1:YTT,"@1;\n", collapse = ""),
      "M1-M",MTT,"@0;\n",
      "Y1-Y",YTT,"@0;\n",
      # "! constrained within-person residual;\n",
      # "eM1-eM",MTT," (sigma_eM);\n",
      # "eY1-eY",YTT," (sigma_eY);\n",
      "!user-specified within-person residual covariance structure;\n",
      within_person_residual_cov,
      "\n",
      
      # level-2 model
      "[IM SM] (a0_IM a0_SM);\n",
      "[IY SY] (b0_IY b0_SY);\n",
      "IM on X (aIM_X);\n",
      "SM on X (aSM_X);\n",
      "IY on X IM SM (bIY_X bIY_IM bIY_SM);\n",
      "SY on X IM SM (bSY_X bSY_IM bSY_SM);\n",
      
      "IM (revarIM);\n",
      "SM (revarSM);\n",
      "IM with SM (IMwSM);\n",
      "IY (revarIY);\n",
      "SY (revarSY);\n",
      "IY with SY (IYwSY);\n",
      "IY SY with IM@0 SM@0;\n",
      
      collapse = ""
    )
    ## mplus level-2 model including Z
    mpluscodeModel_M0_Z=mpluscodeModel_M0
    if( length(Z1names) >=1 ){
      mpluscodeModel_M0_Z = paste0(mpluscodeModel_M0, 
                                   paste0('IM SM on ', Z1names, ";\n"),
                                   paste0('IY SY on ', Z1names, ";\n"),
                                   collapse = ""
      )
    }
    if( length(Z2names)>=1 ){
      mpluscodeModel_M0_Z = paste0(mpluscodeModel_M0_Z, 
                                   paste0('IM SM on ', Z2names, ";\n"),
                                   collapse = ""
      )
    }
    if( length(Z3names)>=1 ){
      mpluscodeIYSYonZ3 = paste0('IY SY on ', Z3names, ";\n")
      mpluscodeModel_M0_Z = paste0(mpluscodeModel_M0_Z, 
                                   paste0('IY SY on ', Z3names, ";\n"),
                                   collapse = ""
      )
    }
    
    naive.mplusModel.Bayes = mplusObject(
      TITLE = "lgcm",
      ANALYSIS = paste0(
        "MODEL = NOCOVARIANCES;
  ESTIMATOR = BAYES; 
  POINT = MEAN;", "\n",
  "THIN = ", THIN,";\n",
  "CHAINS = " , CHAINS,";\n",
  "PROCESS = 2;
  ALGORITHM = GIBBS; 
  BITERATIONS = ", BITERATIONS,";\n",
        collapse = "" ),
      
      MODEL = mpluscodeModel_M0_Z,
      
     # MODELPRIORS = paste0( priors_originalmodel ),
      
      MODELCONSTRAINT = paste0(
        # specific
        "NEW (X_IM_IY, X_SM_IY, X_IM_SY, X_SM_SY  );\n", 
        "X_IM_IY = aIM_X*bIY_IM; \n",
        "X_SM_IY = aSM_X*bIY_SM; \n",
        "X_IM_SY = aIM_X*bSY_IM; \n",
        "X_SM_SY = aSM_X*bSY_SM; \n"
        
      ),
      
      OUTPUT = "TECH1 TECH3 TECH4 TECH8;\n",
      
      usevariables = colnames(dat), autov = F,
      rdata = dat
    )
    
    # sensiC.mplusModel.Bayes ----
    sensiC.mplusModel.Bayes = mplusObject(
      TITLE = "lgcm",
      ANALYSIS = paste0(
        "MODEL = NOCOVARIANCES;
  ESTIMATOR = BAYES; 
  POINT = MEAN;", "\n",
        "THIN = ", THIN,";\n",
        "CHAINS = " , CHAINS,";\n",
        "PROCESS = 2;
  ALGORITHM = GIBBS; 
  BITERATIONS = ", BITERATIONS,";\n",
        collapse = "" ),
      
      MODEL = paste0(mpluscodeModel_M0_Z,
                     "\n",
                     "C by;\n",
                     "[C@0];\n",
                     "C@1;\n",
                     "IM on C (aIM_C);\n",
                     "SM on C (aSM_C);\n",
                     "IY on C (bIY_C);\n",
                     "SY on C (bSY_C);\n"
      ),
      MODELPRIORS = paste0( priors_muC ),
      
      MODELCONSTRAINT = paste0(
        # specific
        "NEW (X_IM_IY, X_SM_IY, X_IM_SY, X_SM_SY  );\n", 
        "X_IM_IY = aIM_X*bIY_IM; \n",
        "X_SM_IY = aSM_X*bIY_SM; \n",
        "X_IM_SY = aIM_X*bSY_IM; \n",
        "X_SM_SY = aSM_X*bSY_SM; \n"
        
      ),
      
      OUTPUT = "TECH1 TECH3 TECH4 TECH8;\n",
      
      usevariables = colnames(dat), autov = F,
      rdata = dat
    )
  }
  
  
  
  ##---fit Model M0----
  set.seed(12345)
  naive.mplusOut.Bayes = mplusModeler(naive.mplusModel.Bayes, run = 1, 
                          writeData = "ifmissing" ,modelout = "naiveBayes.inp")
  
  naive.mplusOut.Bayes$results$tech8$psr-> psr_naive
  if ( psr_naive[nrow(psr_naive), 'psr']<1.1 ) {
    randC_naive.params = as.data.frame(naive.mplusOut.Bayes[["results"]][["parameters"]][["unstandardized"]])  
    abmed_naiveBayes_raw = randC_naive.params[
      c((randC_naive.params$paramHeader=="IM.ON"&randC_naive.params$param=='X') |(randC_naive.params$paramHeader=="SM.ON"&randC_naive.params$param=='X') |
      (randC_naive.params$paramHeader=="IY.ON"&randC_naive.params$param=='IM') | (randC_naive.params$paramHeader=="IY.ON"&randC_naive.params$param=='SM') | 
      (randC_naive.params$paramHeader=="SY.ON"&randC_naive.params$param=='IM') | (randC_naive.params$paramHeader=="SY.ON"&randC_naive.params$param=='SM') |
        randC_naive.params$paramHeader=="New.Additional.Parameters") , ]
    
    abmed_naiveBayes = abmed_naiveBayes_raw[ c(5,6,1,3,2,4,7,9,8,10) ,-1 ]
    abmed_naiveBayes$param[1:6] = paste(abmed_naiveBayes_raw$paramHeader, abmed_naiveBayes_raw$param, sep="_")[c(5,6,1,3,2,4,7,9,8,10)][1:6]
    
    #datpreplist = mget( ls(), envir = environment())
  }

  if (psr_naive[nrow(psr_naive), 'psr'] >= 1.1  ){
    print(naive.mplusOut.Bayes$results$errors)
    print( tail(psr_naive) )
    abmed_naiveBayes='Bayesian estimation of the original model (Model M0) did not coverge (PSR>1.1).'
  }
  
  
  ##--fit Model M1----
  set.seed(23456)
  sensi.mplusOut.Bayes = mplusModeler(sensiC.mplusModel.Bayes, run = 1, 
                                      writeData = "ifmissing" ,modelout = "sensiBayes.inp")
  
  sensi.mplusOut.Bayes$results$tech8$psr-> psr_sensi
  if ( psr_sensi[nrow(psr_sensi), 'psr']<1.1 ) {
    randC_sensi.params = as.data.frame(sensi.mplusOut.Bayes[["results"]][["parameters"]][["unstandardized"]])  
    abmed_sensiBayes_raw = randC_sensi.params[
      c((randC_sensi.params$paramHeader=="IM.ON"&randC_sensi.params$param=='X') |(randC_sensi.params$paramHeader=="SM.ON"&randC_sensi.params$param=='X') |
          (randC_sensi.params$paramHeader=="IY.ON"&randC_sensi.params$param=='IM') | (randC_sensi.params$paramHeader=="IY.ON"&randC_sensi.params$param=='SM') | 
          (randC_sensi.params$paramHeader=="SY.ON"&randC_sensi.params$param=='IM') | (randC_sensi.params$paramHeader=="SY.ON"&randC_sensi.params$param=='SM') |
          randC_sensi.params$paramHeader=="New.Additional.Parameters") , ]
    
    abmed_sensiBayes = abmed_sensiBayes_raw[ c(5,6,1,3,2,4,7,9,8,10) ,-1 ]
    abmed_sensiBayes$param[1:6] = paste(abmed_sensiBayes_raw$paramHeader, abmed_sensiBayes_raw$param, sep="_")[c(5,6,1,3,2,4,7,9,8,10)][1:6]
    
    #datpreplist = mget( ls(), envir = environment())
  }
  
  if (psr_sensi[nrow(psr_sensi), 'psr'] >= 1.1  ){
    print(sensi.mplusOut.Bayes$results$errors)
    print( tail(psr_sensi) )
    abmed_sensiBayes='Bayesian estimation of the sensitivity analysis model (Model M1) did not coverge (PSR>1.1).'
  }
  
  BSAout = data.frame(cbind(abmed_naiveBayes, abmed_sensiBayes[ ,-1]))
  if(ncol(BSAout) > 2){
    colnames(BSAout) = c("param", paste0("Bayes_Original_",colnames(abmed_naiveBayes[,-1]) ), paste0("BSA_",colnames(abmed_sensiBayes[,-1]) ) )
    
  }
  
  BSAoutlist = mget( ls(), envir = environment())
  
 return(BSAoutlist)
}

