
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
  within_person_residual_cov = ("eM1-eM5 (sigma_eM);
 eY1-eY5 (sigma_eY);")
  Xcol = 11 # column for the independent variable 
  Z1cols = c(12) # columns for the level-2 covariates in both the level-2 models of mediator and outcome  
  Z2cols = NULL # columns for the level-2 covariates in the level-2 model of mediator only
  Z3cols = NULL # columns for the level-2 covariates in the level-2 model of outcome only
  alph=0.05 # alpha level per effect
}

fun.fixC.datprep_Z <- function(
  dat = datN200T5, # data frame containing variables in the LGCM model
  Ycols = c(1:5), # columns for repeated measures of outcome from the first to the last time points
  SY.loadings = c(-4,-3,-2,-1,0), # loadings of the outcome slope
  Mcols = c(6:10), # columns for repeated measures of mediator from the first to the last time points
  SM.loadings = c(0,1,2,3,4), # loadings of the mediator slope
  Xcol = 11, # column for the independent variable
  Z1cols = c(12), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
  Z2cols = NULL, # columns for the level-2 covariates in the level-2 model of mediator only
  Z3cols = NULL, # columns for the level-2 covariates in the level-2 model of outcome only
  within_person_residual_cov = ("eM1-eM5 (sigma_eM);
 eY1-eY5 (sigma_eY);"),
  alph=0.05 # alpha level per effect
)
{
  
  YTT = length(Ycols); MTT = length(Mcols)
  colnames(dat)[Ycols] = paste0('Y', 1:length(Ycols) )
  colnames(dat)[Mcols] = paste0('M', 1:length(Mcols) )
  colnames(dat)[Xcol] = 'X'
  Z1names=colnames(dat)[Z1cols]
  Z2names=colnames(dat)[Z2cols]
  Z3names=colnames(dat)[Z3cols]
  
  # original model M0 mplusOut ----
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
    "IM on X (a1_IM);\n",
    "SM on X (a1_SM);\n",
    "IY on X IM SM (b1_IY b2_IY b3_IY);\n",
    "SY on X IM SM (b1_SY b2_SY b3_SY);\n",
    
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
  if( !is.null(Z1cols) ){
    mpluscodeModel_M0_Z = paste0(mpluscodeModel_M0, 
                                 paste0('IM SM on ', Z1names, ";\n"),
                                 paste0('IY SY on ', Z1names, ";\n"),
                                 collapse = ""
                                 )
  }
  if( !is.null(Z2cols) ){
    mpluscodeModel_M0_Z = paste0(mpluscodeModel_M0_Z, 
                                    paste0('IM SM on ', Z2names, ";\n"),
                                 collapse = ""
    )
  }
  if( !is.null(Z3cols) ){
    mpluscodeIYSYonZ3 = paste0('IY SY on ', Z3names, ";\n")
    mpluscodeModel_M0_Z = paste0(mpluscodeModel_M0_Z, 
                                 paste0('IY SY on ', Z3names, ";\n"),
                                 collapse = ""
    )
  }
  
  # original model M0
  naive.mplusModel = mplusObject(
    TITLE = "lgcm",
    ANALYSIS = " 
  MODEL = NOCOVARIANCES;
  ESTIMATOR = ML; 
  ",#ML,a,b uncorrelated #  Starting with Mplus Version 5, TYPE=MEANSTRUCTURE is the default for all analyses.  
    MODEL = mpluscodeModel_M0_Z,
    MODELCONSTRAINT = paste0(
      # for calculate fixC se
      "NEW (detISM);\n", 
      "detISM = revarIM*revarSM - (IMwSM)^2 ;\n",
      "NEW (inISM11, inISM22, inISM12);\n", 
      "inISM11 = revarSM/detISM;\n",
      "inISM22 = revarIM/detISM;\n",
      "inISM12 = (-1)*IMwSM/detISM;\n",
      # specific
      "NEW (X_IM_IY, X_SM_IY, X_IM_SY, X_SM_SY  );\n", 
      "X_IM_IY = a1_IM*b2_IY; \n",
      "X_SM_IY = a1_SM*b3_IY; \n",
      "X_IM_SY = a1_IM*b2_SY; \n",
      "X_SM_SY = a1_SM*b3_SY; \n",
      # total
      "NEW (X_to_IY, X_to_SY  );\n", 
      "X_to_IY = X_IM_IY+X_SM_IY; \n",
      "X_to_SY = X_IM_SY+X_SM_SY; \n"
      
    ),
    
    OUTPUT = "TECH1 TECH3 TECH4;\n",
    
    usevariables = colnames(dat), autov = F,
    rdata = dat
  )
  
  
  mplusOut = mplusModeler(naive.mplusModel, run = 1, 
                          writeData = "ifmissing" ,modelout = "naive.inp")
  
  
  if ( length(mplusOut$results$errors)+length(mplusOut$results$warnings) == 0 ) {
    sample_naive.Parameters = readModels("naive.out",what = c("parameters")
                                         ,quiet = T)$parameters$unstandardized
    
    # for range limits  ----
    #ESTIMATED CORRELATION MATRIX FOR THE LATENT VARIABLES
    latCorEst=mplusOut$results$tech4$latCorEst
    latCorEst[upper.tri(latCorEst)] =
      t(latCorEst)[upper.tri( t(latCorEst) )] 
    # order: c('IM','SM','IY','SY','X','Z')
    Rc = latCorEst[c('IM','SM','IY','SY','X',Z1names,Z2names,Z3names),c('IM','SM','IY','SY','X',Z1names,Z2names,Z3names)]
    # for transform ----
    latCovEst=mplusOut$results$tech4$latCovEst
    latCovEst[upper.tri(latCovEst)] =
      t(latCovEst)[upper.tri( t(latCovEst) )] 
    sd_ISM = sqrt(diag(latCovEst[c('IM','SM'), c('IM','SM')]))
    sd_ISY = sqrt(diag(latCovEst[c('IY','SY'), c('IY','SY')]))
    COVc = latCovEst[c('IM','SM','IY','SY','X',Z1names,Z2names,Z3names),c('IM','SM','IY','SY','X',Z1names,Z2names,Z3names)]
    
    # for se(hat.b_fixC) 
    #ESTIMATED COVARIANCE MATRIX FOR PARAMETER ESTIMATES
    # the joint distribution of naive a,b,inISM ----
    numapath = mplusOut$results$tech1$parameterSpecification$X$beta[c('IM','SM'),c('X')]
    numbpath = mplusOut$results$tech1$parameterSpecification$X$beta[c('IY','SY'),c('IM','SM')]
    numbpath = c(numbpath[1,], numbpath[2,])
    numinISM = max(mplusOut$results$tech1$parameterSpecification$X$psi,na.rm = T)+c(2:4)
    #V.a_b_inISM
    V.a_b_inISM=mplusOut$results$tech3$paramCov[c(numapath,numbpath,numinISM),c(numapath,numbpath,numinISM)]
    rownames(V.a_b_inISM)=c("IM.ON_X","SM.ON_X",
                            "IY.ON_IM","IY.ON_SM","SY.ON_IM","SY.ON_SM",
                            "inISM11", "inISM22", "inISM12")
    colnames(V.a_b_inISM)=rownames(V.a_b_inISM)
    V.a_b_inISM[upper.tri(V.a_b_inISM)] =t(V.a_b_inISM)[upper.tri( t(V.a_b_inISM) )] 
    #E.a_b_inISM
    wsample_naive.Parameters = as_tibble(sample_naive.Parameters) %>%
      filter(!str_detect(paramHeader,".BY")) %>%
      filter(!str_detect(param,"E")) %>%
      unite("paramName",paramHeader,param) %>% select(paramName,est) %>%
      spread(paramName,est)
    E.a_b_inISM=wsample_naive.Parameters[,c("IM.ON_X","SM.ON_X",
                                            "IY.ON_IM","IY.ON_SM","SY.ON_IM","SY.ON_SM",
                                            "New.Additional.Parameters_INISM11", "New.Additional.Parameters_INISM22", "New.Additional.Parameters_INISM12")]
    E.a_b_inISM = unlist(E.a_b_inISM)
    names(E.a_b_inISM)=c("IM.ON_X","SM.ON_X",
                         "IY.ON_IM","IY.ON_SM","SY.ON_IM","SY.ON_SM",
                         "inISM11", "inISM22", "inISM12")
    est_b = matrix(E.a_b_inISM[c("IY.ON_IM","IY.ON_SM","SY.ON_IM","SY.ON_SM")], 2,2,byrow = T)
    estVarCov.b_inISM=V.a_b_inISM[-c(1:2),-c(1:2)]
    estE.b_inISM=E.a_b_inISM[-c(1:2)]
    estEE.b_inISM = estE.b_inISM %*% t(estE.b_inISM)
    names(estE.b_inISM)=c("IY.ON_IM","IY.ON_SM","SY.ON_IM","SY.ON_SM",
                          "inISM11", "inISM22", "inISM12")
    rownames(estEE.b_inISM)=c("IY.ON_IM","IY.ON_SM","SY.ON_IM","SY.ON_SM",
                              "inISM11", "inISM22", "inISM12")
    colnames(estEE.b_inISM)=rownames(estEE.b_inISM)
    
    
    # two-tailed wald test critical value for z-statistic est_se
    # alph=0.05 # if multiple testing, may need Bonferroni correction (0.05/#tests)
    alph=alph
    qchi = (qnorm(alph/2))^2 # or qchisq(1-alph, 1)
    
    critA = estEE.b_inISM - qchi*estVarCov.b_inISM
    
    # muaC in the order: muaC_IM muaC_SM
    #IY.ON
    critA.IY.ON_IM_fixC = critA[
      c("IY.ON_IM","inISM11","inISM12"), c("IY.ON_IM","inISM11","inISM12") ]
    estE.IY.ON_IM_fixC = estE.b_inISM[ c("IY.ON_IM","inISM11","inISM12") ] # solvebC.bsig( critA.IY.ON_IM_fixC, muaC = c(.3,.2)) # use this bC in fixC, Mplus (est_se)=1.96, pval=0.05
    
    critA.IY.ON_SM_fixC = critA[
      c("IY.ON_SM","inISM12","inISM22"), c("IY.ON_SM","inISM12","inISM22") ]
    estE.IY.ON_SM_fixC = estE.b_inISM[ c("IY.ON_SM","inISM12","inISM22") ]
    
    # SY.ON
    critA.SY.ON_IM_fixC = critA[
      c("SY.ON_IM","inISM11","inISM12"), c("SY.ON_IM","inISM11","inISM12") ]
    estE.SY.ON_IM_fixC = estE.b_inISM[ c("SY.ON_IM","inISM11","inISM12") ]
    
    critA.SY.ON_SM_fixC = critA[
      c("SY.ON_SM","inISM12","inISM22"), c("SY.ON_SM","inISM12","inISM22") ]
    estE.SY.ON_SM_fixC = estE.b_inISM[ c("SY.ON_SM","inISM12","inISM22") ]
    datpreplist = mget( ls(), envir = environment())
  }

  if ( length(mplusOut$results$errors)+length(mplusOut$results$warnings) != 0 ){
    print(mplusOut$results$errors)
    print(mplusOut$results$warnings)
    datpreplist='original model M0 returns errors and/or warnings.'
  }
  
 return(datpreplist)
}

