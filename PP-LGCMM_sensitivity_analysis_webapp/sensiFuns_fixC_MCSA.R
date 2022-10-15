
defaultW <- getOption("warn")
options(warn = -1)

library(mvtnorm)
library(ggplot2)
library(ggpubr)
library(tidyverse)
options(tibble.print_max = Inf, tibble.width = Inf)
library(MplusAutomation)

# source('fun_fixC_datprep_Z.R')
# source('sourcefun_fixC_sr_Z.R')
# source('sourcefun_plot_fixC_sr_Z.R')

# datprep_list = fun.fixC.datprep_Z( dat = read.csv('LGCMdata.csv',header = T)  )


######## Run Sensitivity Analysis ##########

######### FixC_sr_Z (JS test, based on analytical Solution, transform rC to muC) ########
# in the same order as Rc!!!: rCIm,rCSm,rCIy,rCSy,X,Z1,Z2,Z3
fun.sensifixC_JS_Z <-function( datpreplist, rCvalues=c(0.2,0.2,0.2,0.2) ){ # rCvalues: rCIm,rCSm,rCIy,rCSy used in the panel names 
  list2env(datpreplist, envir = environment())
  
  # plausible rC values
  smallfixC_s.onSm.all_rC = expand.grid(
    rCIm = c(0, rCvalues[1]) # corr(C, Im) # subplot labels
    , rCSm = seq(-1,1,by=0.01) # need a finer grid around zero to plot the vertical "asymptote" of y=1/x (x approaches 0, y approaches +-Inf)
    , rCIy = c(0, rCvalues[3]) # subplot labels
    , rCSy = c(0, rCvalues[4]) # no need, because it will be solved; but rCSy still influences range limits of rCIy; not sig or zero lines
    , rCX = c(0) # in the same order as Rc
    # , rCZ = c(0) # in the same order as Rc
  )
  smallfixC_s.onIm.all_rC = expand.grid(
    rCIm = seq(-1,1,by=0.01) # x-axis 
    #c(0, 0.3, 0.5)
    , rCSm = c(0,rCvalues[2]) # residual variance of SM in model M0 is not significant, so unlikely to be highly correlated with C 
    , rCIy = c(0,rCvalues[3]) 
    , rCSy = c(0,rCvalues[4])
    , rCX = c(0) # in the same order as Rc
    # , rCZ = c(0) # in the same order as Rc
  )
  # medium large
  medilargfixC_s.onSm.all_rC = expand.grid(
    rCIm = c(0.3, 0.5) # corr(C, Im) # subplot labels
    , rCSm = seq(-1,1,by=0.01) # need a finer grid around zero to plot the vertical "asymptote" of y=1/x (x approaches 0, y approaches +-Inf)
    , rCIy = c(0.3, 0.5) # subplot labels
    , rCSy = c(0.3, 0.5) # no need, because it will be solved; but rCSy still influences range limits of rCIy; not sig or zero lines
    , rCX = c(0) # in the same order as Rc
    # , rCZ = c(0) # in the same order as Rc
  )
  
  medilargfixC_s.onIm.all_rC = expand.grid(
    rCIm = seq(-1,1,by=0.01) # x-axis
    #c(0, 0.3, 0.5)
    , rCSm = c(0.3, 0.5) # residual variance of SM in model M0 is not significant, so unlikely to be highly correlated with C
    , rCIy = c(0.3, 0.5)
    , rCSy = c(0.3, 0.5)
    , rCX = c(0) # in the same order as Rc
    # , rCZ = c(0) # in the same order as Rc
  )
  
  if ( (length(Z1cols)+length(Z2cols)+length(Z3cols))>=1 ){
    rCZ123=matrix(0, nrow = nrow(smallfixC_s.onSm.all_rC), ncol = (length(Z1cols)+length(Z2cols)+length(Z3cols)) )
    colnames(rCZ123)=paste0('rC',c(Z1names,Z2names,Z3names))
    smallfixC_s.onSm.all_rCZ=cbind(smallfixC_s.onSm.all_rC, rCZ123)
    smallfixC_s.onIm.all_rCZ=cbind(smallfixC_s.onIm.all_rC,rCZ123)
    medilargfixC_s.onSm.all_rCZ=cbind(medilargfixC_s.onSm.all_rC, rCZ123)
    medilargfixC_s.onIm.all_rCZ=cbind(medilargfixC_s.onIm.all_rC, rCZ123)
    
  }
  if ( (length(Z1cols)+length(Z2cols)+length(Z3cols)) < 1 ){
    smallfixC_s.onSm.all_rCZ=smallfixC_s.onSm.all_rC
    smallfixC_s.onIm.all_rCZ=smallfixC_s.onIm.all_rC
    medilargfixC_s.onSm.all_rCZ=medilargfixC_s.onSm.all_rC
    medilargfixC_s.onIm.all_rCZ=medilargfixC_s.onIm.all_rC
  }
  

  smallfixC_s.onSm.SensiRes = fixC_s.SensiRes(
    all_rC = smallfixC_s.onSm.all_rCZ, envlist = datpreplist )
  
  smallfixC_s.onIm.SensiRes = fixC_s.SensiRes(
    all_rC = smallfixC_s.onIm.all_rCZ, envlist = datpreplist)

  
  medilargfixC_s.onSm.SensiRes = fixC_s.SensiRes(
    all_rC = medilargfixC_s.onSm.all_rCZ, envlist = datpreplist)
  
  medilargfixC_s.onIm.SensiRes = fixC_s.SensiRes(
    all_rC = medilargfixC_s.onIm.all_rCZ, envlist = datpreplist)
  
  # output results ----
  out.sensifixC_JS_Z=list(
    # original model ML estimation results
    mplusOut, sample_naive.Parameters, 
    # on Im sensi results
    smallfixC_s.onIm.SensiRes, medilargfixC_s.onIm.SensiRes,
    # on Sm sensi results
    smallfixC_s.onSm.SensiRes, medilargfixC_s.onSm.SensiRes
  )
  names(out.sensifixC_JS_Z) = c(
    'mplusOut', 'sample_naive.Parameters', 
    # on Im sensi results
    'smallfixC_s.onIm.SensiRes', 'medilargfixC_s.onIm.SensiRes',
    # on Sm sensi results
    'smallfixC_s.onSm.SensiRes', 'medilargfixC_s.onSm.SensiRes'
  )
  
  out = mget( ls(), envir = environment())
  return(out)
}

# out_sensifixC_JS_Z=fun.sensifixC_JS_Z( datpreplist = datprep_list )


# plot results ----
fun.plotsensifixC_JS_Z <- function(outsensifixC_JS_Z, alph=0.05 ){
  
  # list2env(datpreplist, envir = environment())
  list2env(outsensifixC_JS_Z, envir = environment())
  
  # Sy on Im ----
  pval_apath = sample_naive.Parameters[sample_naive.Parameters$paramHeader=="IM.ON"&sample_naive.Parameters$param=="X", c("pval")]
  pval_bpath = sample_naive.Parameters[sample_naive.Parameters$paramHeader=="SY.ON"&sample_naive.Parameters$param=="IM", c("pval")]
  
  if( max(pval_apath, pval_bpath)>=alph ){ifsig_original=FALSE}
  if( max(pval_apath, pval_bpath)<alph ){ifsig_original=TRUE}
  est_bpath=sample_naive.Parameters[sample_naive.Parameters$paramHeader=="SY.ON"&sample_naive.Parameters$param=="IM", c("est")]
  dire_bpath= (est_bpath>0)
  
  # zero small
  zerosmallplot=T
  if( zerosmallplot ){
    plotdat_all = as_tibble(smallfixC_s.onIm.SensiRes) %>%
      select( rCIm,rCSm,rCIy,
              range.lo_rCSy, range.up_rCSy,
              contains(".SyonIm")) %>%
      mutate( noconfounding.SyonIm = -(1/sd_ISM[1])*(rCSm*sd_ISM[2])*estE.b_inISM['inISM12']/estE.b_inISM['inISM11'] )
  }
  colnames(plotdat_all)=gsub("_","",gsub("rCSy", "", colnames(plotdat_all)))
  colnames(plotdat_all)=gsub(".SyonIm", "", colnames(plotdat_all))
  xaxis="rCIm"
  yaxis="rCSy"
  facetgrid="rCSm ~ rCIy"
  if( ifsig_original ){
    sensiplot01_SyonIm =sensiplot.sig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                      xaxis="rCIm", yaxis="rCSy",facetgrid="rCSm ~ rCIy")
  }
  if( !ifsig_original ){
    sensiplot01_SyonIm =sensiplot.nonsig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                         xaxis="rCIm", yaxis="rCSy",facetgrid="rCSm ~ rCIy")
  }
  
  # Sy on Im ----
  pval_apath = sample_naive.Parameters[sample_naive.Parameters$paramHeader=="IM.ON"&sample_naive.Parameters$param=="X", c("pval")]
  pval_bpath = sample_naive.Parameters[sample_naive.Parameters$paramHeader=="SY.ON"&sample_naive.Parameters$param=="IM", c("pval")]
  
  if( max(pval_apath, pval_bpath)>=alph ){ifsig_original=FALSE}
  if( max(pval_apath, pval_bpath)<alph ){ifsig_original=TRUE}
  est_bpath=sample_naive.Parameters[sample_naive.Parameters$paramHeader=="SY.ON"&sample_naive.Parameters$param=="IM", c("est")]
  dire_bpath= (est_bpath>0)
  
  # zero small
  zerosmallplot=T
  if( zerosmallplot ){
    plotdat_all = as_tibble(smallfixC_s.onIm.SensiRes) %>%
      select( rCIm,rCSm,rCIy,
              range.lo_rCSy, range.up_rCSy,
              contains(".SyonIm")) %>%
      mutate( noconfounding.SyonIm = -(1/sd_ISM[1])*(rCSm*sd_ISM[2])*estE.b_inISM['inISM12']/estE.b_inISM['inISM11'] )
  }
  colnames(plotdat_all)=gsub("_","",gsub("rCSy", "", colnames(plotdat_all)))
  colnames(plotdat_all)=gsub(".SyonIm", "", colnames(plotdat_all))
  xaxis="rCIm"
  yaxis="rCSy"
  facetgrid="rCSm ~ rCIy"
  if( ifsig_original ){
    sensiplot01_SyonIm =sensiplot.sig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                      xaxis="rCIm", yaxis="rCSy",facetgrid="rCSm ~ rCIy")
  }
  if( !ifsig_original ){
    sensiplot01_SyonIm =sensiplot.nonsig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                         xaxis="rCIm", yaxis="rCSy",facetgrid="rCSm ~ rCIy")
  }
  # medi large
  medilargplot=T
  if( medilargplot ){
    plotdat_all = as_tibble(medilargfixC_s.onIm.SensiRes) %>%
      select( rCIm,rCSm,rCIy,
              range.lo_rCSy, range.up_rCSy,
              contains(".SyonIm")) %>%
      mutate( noconfounding.SyonIm = -(1/sd_ISM[1])*(rCSm*sd_ISM[2])*estE.b_inISM['inISM12']/estE.b_inISM['inISM11'] )
    
  }
  colnames(plotdat_all)=gsub("_","",gsub("rCSy", "", colnames(plotdat_all)))
  colnames(plotdat_all)=gsub(".SyonIm", "", colnames(plotdat_all))
  xaxis="rCIm"
  yaxis="rCSy"
  facetgrid="rCSm ~ rCIy"
  if( ifsig_original ){
    sensiplot35_SyonIm =sensiplot.sig(plotdat_all=plotdat_all, dire_bpath=dire_bpath, 
                                      xaxis="rCIm", yaxis="rCSy",facetgrid="rCSm ~ rCIy")
  }
  if( !ifsig_original ){
    sensiplot35_SyonIm =sensiplot.nonsig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                         xaxis="rCIm", yaxis="rCSy",facetgrid="rCSm ~ rCIy")
  }
  
  
  # Iy on Im ----
  pval_apath = sample_naive.Parameters[sample_naive.Parameters$paramHeader=="IM.ON"&sample_naive.Parameters$param=="X", c("pval")]
  pval_bpath = sample_naive.Parameters[sample_naive.Parameters$paramHeader=="IY.ON"&sample_naive.Parameters$param=="IM", c("pval")]
  if( max(pval_apath, pval_bpath)>=alph ){ifsig_original=FALSE}
  if( max(pval_apath, pval_bpath)<alph ){ifsig_original=TRUE}
  
  est_bpath=sample_naive.Parameters[sample_naive.Parameters$paramHeader=="IY.ON"&sample_naive.Parameters$param=="IM", c("est")] 
  dire_bpath = (est_bpath>=0)
  # zero small
  if( zerosmallplot ){
    plotdat_all = as_tibble(smallfixC_s.onIm.SensiRes) %>%
      select( rCIm,rCSm,rCSy,
              range.lo_rCIy, range.up_rCIy,
              contains(".IyonIm")) %>%
      mutate( noconfounding.IyonIm = -(1/sd_ISM[1])*(rCSm*sd_ISM[2])*estE.b_inISM['inISM12']/estE.b_inISM['inISM11'] )
  }
  colnames(plotdat_all)=gsub("_","",gsub("rCIy", "", colnames(plotdat_all)))
  colnames(plotdat_all)=gsub(".IyonIm", "", colnames(plotdat_all))
  xaxis="rCIm"
  yaxis="rCIy"
  facetgrid="rCSm ~ rCSy"
  if( ifsig_original ){
    sensiplot01_IyonIm= sensiplot.sig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                      xaxis="rCIm", yaxis="rCIy",facetgrid="rCSm ~ rCSy")
    
  }
  if( !ifsig_original ){
    sensiplot01_IyonIm= sensiplot.nonsig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                         xaxis="rCIm", yaxis="rCIy",facetgrid="rCSm ~ rCSy")
  }
  # medi large
  if( medilargplot ){
    plotdat_all = as_tibble(medilargfixC_s.onIm.SensiRes) %>%
      select( rCIm,rCSm,rCSy,
              range.lo_rCIy, range.up_rCIy,
              contains(".IyonIm")) %>%
      mutate( noconfounding.IyonIm = -(1/sd_ISM[1])*(rCSm*sd_ISM[2])*estE.b_inISM['inISM12']/estE.b_inISM['inISM11'] )
  }
  colnames(plotdat_all)=gsub("_","",gsub("rCIy", "", colnames(plotdat_all)))
  colnames(plotdat_all)=gsub(".IyonIm", "", colnames(plotdat_all))
  xaxis="rCIm"
  yaxis="rCIy"
  facetgrid="rCSm ~ rCSy"
  if( ifsig_original ){
    sensiplot35_IyonIm= sensiplot.sig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                      xaxis="rCIm", yaxis="rCIy",facetgrid="rCSm ~ rCSy")
    
  }
  if( !ifsig_original ){
    sensiplot35_IyonIm= sensiplot.nonsig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                         xaxis="rCIm", yaxis="rCIy",facetgrid="rCSm ~ rCSy")
    
  }
  
  
  # Sy on Sm ----
  pval_apath = sample_naive.Parameters[sample_naive.Parameters$paramHeader=="SM.ON"&sample_naive.Parameters$param=="X", c("pval")]
  pval_bpath = sample_naive.Parameters[sample_naive.Parameters$paramHeader=="SY.ON"&sample_naive.Parameters$param=="SM", c("pval")]
  
  if( max(pval_apath, pval_bpath)>=alph ){ifsig_original=FALSE}
  if( max(pval_apath, pval_bpath)<alph ){ifsig_original=TRUE}
  est_b=sample_naive.Parameters[sample_naive.Parameters$paramHeader=="SY.ON"&sample_naive.Parameters$param=="SM", c("est")]
  dire_bpath=(est_b>0)
  
  if(zerosmallplot){
    plotdat_all = as_tibble(smallfixC_s.onSm.SensiRes) %>%
      select( rCIm,rCSm,rCIy
              , range.lo_rCSy, range.up_rCSy,
              contains(".SyonSm")) %>%
      mutate( noconfounding.SyonSm = -(1/sd_ISM[2])*(rCIm*sd_ISM[1])*estE.b_inISM['inISM12']/estE.b_inISM['inISM22'] )
  }
  
  colnames(plotdat_all)=gsub("_","",gsub("rCSy", "", colnames(plotdat_all)))
  colnames(plotdat_all)=gsub(".SyonSm", "", colnames(plotdat_all))
  xaxis="rCSm"
  yaxis="rCSy"
  facetgrid="rCIm ~ rCIy"
  if( ifsig_original ) {
    sensiplot01_SyonSm=sensiplot.sig( plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                      xaxis="rCSm", yaxis="rCSy",facetgrid="rCIm ~ rCIy" )
  }
  if( !ifsig_original ) {
    sensiplot01_SyonSm=sensiplot.nonsig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                         xaxis="rCSm", yaxis="rCSy",facetgrid="rCIm ~ rCIy" )
  }
  # medi large
  if(medilargplot){
    plotdat_all = as_tibble(medilargfixC_s.onSm.SensiRes) %>%
      select( rCIm,rCSm,rCIy
              , range.lo_rCSy, range.up_rCSy,
              contains(".SyonSm")) %>%
      mutate( noconfounding.SyonSm = -(1/sd_ISM[2])*(rCIm*sd_ISM[1])*estE.b_inISM['inISM12']/estE.b_inISM['inISM22'] )
  }
  
  colnames(plotdat_all)=gsub("_","",gsub("rCSy", "", colnames(plotdat_all)))
  colnames(plotdat_all)=gsub(".SyonSm", "", colnames(plotdat_all))
  xaxis="rCSm"
  yaxis="rCSy"
  facetgrid="rCIm ~ rCIy"
  
  if( ifsig_original ) {
    sensiplot35_SyonSm=sensiplot.sig( plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                      xaxis="rCSm", yaxis="rCSy",facetgrid="rCIm ~ rCIy" )
  }
  if( !ifsig_original ) {
    sensiplot35_SyonSm=sensiplot.nonsig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                         xaxis="rCSm", yaxis="rCSy",facetgrid="rCIm ~ rCIy" )
  }
  
  
  # Iy on Sm ----
  pval_apath = sample_naive.Parameters[sample_naive.Parameters$paramHeader=="SM.ON"&sample_naive.Parameters$param=="X", c("pval")]
  pval_bpath = sample_naive.Parameters[sample_naive.Parameters$paramHeader=="IY.ON"&sample_naive.Parameters$param=="SM", c("pval")]
  
  if( max(pval_apath, pval_bpath)>=alph ){ifsig_original=FALSE}
  if( max(pval_apath, pval_bpath)<alph ){ifsig_original=TRUE}
  est_b=sample_naive.Parameters[sample_naive.Parameters$paramHeader=="IY.ON"&sample_naive.Parameters$param=="SM", c("est")]
  dire_bpath=(est_b>0)
  
  xaxis="rCSm"
  yaxis="rCIy"
  facetgrid="rCIm ~ rCSy"
  
  if(zerosmallplot){
    plotdat_all = as_tibble(smallfixC_s.onSm.SensiRes) %>%
      select( rCIm,rCSm,rCSy
              , range.lo_rCIy, range.up_rCIy,
              contains(".IyonSm")) %>%
      mutate( noconfounding.IyonSm = -(1/sd_ISM[2])*(rCIm*sd_ISM[1])*estE.b_inISM['inISM12']/estE.b_inISM['inISM22'] ) 
  }
  colnames(plotdat_all)=gsub("_","",gsub("rCIy", "", colnames(plotdat_all)))
  colnames(plotdat_all)=gsub(".IyonSm", "", colnames(plotdat_all))
  if( ifsig_original ){
    sensiplot01_IyonSm= sensiplot.sig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                      xaxis="rCSm", yaxis="rCIy",facetgrid="rCIm ~ rCSy")
    
  }
  if( !ifsig_original ){
    sensiplot01_IyonSm= sensiplot.nonsig(plotdat_all=plotdat_all, dire_bpath=dire_bpath, 
                                         xaxis="rCSm", yaxis="rCIy",facetgrid="rCIm ~ rCSy")
    
  }
  # medi large
  if(medilargplot){
    plotdat_all = as_tibble(medilargfixC_s.onSm.SensiRes) %>%
      select( rCIm,rCSm,rCSy
              , range.lo_rCIy, range.up_rCIy,
              contains(".IyonSm")) %>%
      mutate( noconfounding.IyonSm = -(1/sd_ISM[2])*(rCIm*sd_ISM[1])*estE.b_inISM['inISM12']/estE.b_inISM['inISM22'] ) 
  }
  colnames(plotdat_all)=gsub("_","",gsub("rCIy", "", colnames(plotdat_all)))
  colnames(plotdat_all)=gsub(".IyonSm", "", colnames(plotdat_all))
  if( ifsig_original ){
    sensiplot35_IyonSm= sensiplot.sig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                      xaxis="rCSm", yaxis="rCIy",facetgrid="rCIm ~ rCSy")
    
  }
  if( !ifsig_original ){
    sensiplot35_IyonSm= sensiplot.nonsig(plotdat_all=plotdat_all, dire_bpath=dire_bpath,
                                         xaxis="rCSm", yaxis="rCIy",facetgrid="rCIm ~ rCSy")
    
  }
  
  outplots=list(
    sensiplot01_IyonIm, sensiplot35_IyonIm,
    sensiplot01_SyonIm, sensiplot35_SyonIm,
    sensiplot01_IyonSm, sensiplot35_IyonSm,
    sensiplot01_SyonSm, sensiplot35_SyonSm
  )
  names(outplots)=c(
    'sensiplot01_IyonIm', 'sensiplot35_IyonIm',
    'sensiplot01_SyonIm', 'sensiplot35_SyonIm',
    'sensiplot01_IyonSm', 'sensiplot35_IyonSm',
    'sensiplot01_SyonSm', 'sensiplot35_SyonSm'
  )
  
  sensiplots_IyonIm = annotate_figure(ggarrange(
    sensiplot01_IyonIm, sensiplot35_IyonIm, common.legend = T , nrow=1)
    ,top = text_grob('Sensitivity plots for the indirect effect of X on the outcome intercept (Iy) through the mediator intercept (Im)'
                     ,face = 'bold', size = 14) )
  sensiplots_SyonIm = annotate_figure(ggarrange(
    sensiplot01_SyonIm, sensiplot35_SyonIm, common.legend = T , nrow=1)
    ,top = text_grob('Sensitivity plots for the indirect effect of X on the outcome slope (Sy) through the mediator intercept (Im)'
                     ,face = 'bold', size = 14) )
  sensiplots_IyonSm = annotate_figure(ggarrange(
    sensiplot01_IyonSm, sensiplot35_IyonSm, common.legend = T , nrow=1)
    ,top = text_grob('Sensitivity plots for the indirect effect of X on the outcome intercept (Iy) through the mediator slope (Sm)'
                     ,face = 'bold', size = 14) )
  sensiplots_SyonSm = annotate_figure(ggarrange(
    sensiplot01_SyonSm, sensiplot35_SyonSm, common.legend = T , nrow=1)
    ,top = text_grob('Sensitivity plots for the indirect effect of X on the outcome intercept (Sy) through the mediator intercept (Sy)'
                     ,face = 'bold', size = 14) )
  
  sensiplots_all=ggarrange(sensiplots_IyonIm,sensiplots_SyonIm,
                           sensiplots_IyonSm,sensiplots_SyonSm, nrow = 4)
  
  out = mget( ls(), envir = environment())
  
  return(out)
}

# outsensi_ResPlot=fun.plotsensifixC_JS_Z(outsensifixC_JS_Z = out_sensifixC_JS_Z )

###################### HybridC_sr #########
fun.hybridmuC_sr <- function(datpreplist, Mean_C=c(0.39, 0.39, -0.39, -0.39), 
         Sigma_C=diag(rep(0.5^2), nrow = 4),
         K=1000, nboot=1, seed=12
         ){
  list2env(datpreplist, envir = environment())
  
  # (Step.1) specify a joint informative priors ----
  #[e.g., Multivariate Normal for aC,bC; Uniform for rC] for the sensitivity parameters to account for uncertainty. 
 
  # Mean_C=c(0.39, 0.39, -0.39, -0.39)
  names(Mean_C)=c('muaC.Im','muaC.Sm','mubC.Iy','mubC.Sy')
  # Sigma_C=diag(rep(0.5^2), nrow = 4)
  
  # (Step.2) draw K (e.g., K=1000) sets of sensitivity parameters from the specified prior. ----
  # K=100
  # set.seed(seed = 12)
  set.seed(seed = seed)
  # prior of muC
  muC.Kdraws=rmvnorm(K, Mean_C, Sigma_C)
  
  # (Step.3) for each set of drawn sensitivity parameters, ----
  #obtain frequentist MLE estimates (and SEs, z statistics, p-values) for the indirect effects in the sensitivity analysis model
  muC1.hybrid <- function( i=1 , nboot=1 ){
    muaC=muC.Kdraws[i, 1:2]
    mubC_IY=muC.Kdraws[i,3]
    mubC_SY=muC.Kdraws[i,4]
    
    comb=matrix(0, 6, ncol(V.a_b_inISM),
                dimnames = list(paste0("fixC_",colnames(V.a_b_inISM)[1:6]),
                                paste0("naive_",colnames(V.a_b_inISM))
                ) )
    # path a
    comb[c("fixC_IM.ON_X"),c("naive_IM.ON_X")]=1 
    comb[c("fixC_SM.ON_X"),c("naive_SM.ON_X")]=1
    # path b
    #Iy
    comb[c("fixC_IY.ON_IM"),c("naive_IY.ON_IM","naive_inISM11","naive_inISM12")]=
      c(1, -mubC_IY*muaC)
    comb[c("fixC_IY.ON_SM"),c("naive_IY.ON_SM","naive_inISM12","naive_inISM22")]=
      c(1, -mubC_IY*muaC)
    #Sy
    comb[c("fixC_SY.ON_IM"),c("naive_SY.ON_IM","naive_inISM11","naive_inISM12")]=
      c(1, -mubC_SY*muaC)
    comb[c("fixC_SY.ON_SM"),c("naive_SY.ON_SM","naive_inISM12","naive_inISM22")]=
      c(1, -mubC_SY*muaC)
    # est, se of paths a b ----
    fixC_E.a_b = comb %*% E.a_b_inISM
    fixC_V.a_b = comb %*% V.a_b_inISM %*% t(comb)
    fixC_SE.a_b = sqrt(diag(fixC_V.a_b))
    # fixC_Corr.a_b = diag(1/fixC_se.a_b,length(fixC_se.a_b))%*%fixC_V.a_b%*%diag(1/fixC_se.a_b,length(fixC_se.a_b))
    fixC_E.med = as_tibble( t(fixC_E.a_b) ) %>% 
      mutate( fixC_X_IM_IY = fixC_IM.ON_X * fixC_IY.ON_IM ) %>% 
      mutate( fixC_X_SM_IY = fixC_SM.ON_X * fixC_IY.ON_SM ) %>%
      mutate(fixC_X_IM_SY = fixC_IM.ON_X * fixC_SY.ON_IM ) %>%
      mutate(fixC_X_SM_SY = fixC_SM.ON_X * fixC_SY.ON_SM ) %>%
      # total
      mutate(fixC_X_TO_IY = fixC_X_IM_IY + fixC_X_SM_IY) %>%
      mutate(fixC_X_TO_SY = fixC_X_IM_SY + fixC_X_SM_SY) %>% 
      gather(key = "paramName","est")
    
    # # mcci ----
    nboot=1
    onebootfixC_a_b = rmvnorm(nboot, mean = fixC_E.a_b, sigma = fixC_V.a_b)
    colnames(onebootfixC_a_b)=colnames(fixC_V.a_b)
    
    onebootfixC_params = as_tibble(onebootfixC_a_b) %>% 
      mutate( fixC_X_IM_IY = fixC_IM.ON_X * fixC_IY.ON_IM ) %>% 
      mutate( fixC_X_SM_IY = fixC_SM.ON_X * fixC_IY.ON_SM ) %>%
      mutate(fixC_X_IM_SY = fixC_IM.ON_X * fixC_SY.ON_IM ) %>%
      mutate(fixC_X_SM_SY = fixC_SM.ON_X * fixC_SY.ON_SM ) %>%
      # total
      mutate(fixC_X_TO_IY = fixC_X_IM_IY + fixC_X_SM_IY) %>%
      mutate(fixC_X_TO_SY = fixC_X_IM_SY + fixC_X_SM_SY) 
    
    # transform to rC to see if in the range limit ----
    est_b=matrix(fixC_E.a_b[3:6], 2,2,byrow = T)
    rownames(est_b)=c('IY.ON','SY.ON')
    colnames(est_b)=c('_IM','_SM')
    rCISm = muaC/sd_ISM
    rCIy = (t(est_b[1, ]) %*% muaC + mubC_IY) / (sd_ISY[1])
    rCSy = (t(est_b[2, ]) %*% muaC + mubC_SY) / (sd_ISY[2])
    Rlatent=diag(1, 1+nrow(Rc),1+ncol(Rc))
    Rlatent[1,2:5]=c(rCISm,rCIy,rCSy)
    Rlatent[2:5,1]=c(rCISm,rCIy,rCSy)
    Rlatent[-1,-1]=Rc
    colnames(Rlatent)=c('C',colnames(Rc))
    rownames(Rlatent)=colnames(Rlatent)
    ifrangelim = (det(Rlatent)>0)
    
    # # gather output
    # onesr_hybrid.a_b = data.frame(paramName = rownames(fixC_V.a_b), est = fixC_E.a_b, se = fixC_SE.a_b, ifrangelim = ifrangelim)
    # onesr_hybrid.med = full_join(fixC_E.med[7:12, ], fixC_MCCI[7:12, ], by = 'paramName') %>%
    #   mutate( ifrangelim = ifrangelim )
    # 
    # onesr_hybrid.out = list( onesr_hybrid.a_b, onesr_hybrid.med )
    # names(onesr_hybrid.out) = c( 'a_b', 'med' )
    
    # gather output
    onesr_hybrid.out=onebootfixC_params %>% mutate( ifrangelim = ifrangelim )
    
    return(onesr_hybrid.out)
  }
  
  # K=1000
  # set.seed(seed = 12)
  muC.Kdraws=rmvnorm(K, Mean_C, Sigma_C)
  # nboot=1000
  hybridres = lapply(1:K, muC1.hybrid, nboot=1)
  mhybrid_res=data.frame(do.call(rbind, hybridres))
  # Mean_C2rC = c(Mean_C[1:2]/sd_ISM,
  #               bC_IY.to.rCIy(mubC_IY = Mean_C[3], mubC_SY = Mean_C[4], muaC = Mean_C[1:2]),
  #               bC_SY.to.rCSy(mubC_SY = Mean_C[4], mubC_IY = Mean_C[3], muaC = Mean_C[1:2])
  # )
  out=mget(ls(), envir = environment())
  return(out)
}


fun.hybridrC_sr <- function(
  datpreplist, 
  K=1000,Z2names=NULL,
  Min_rC=c(0.3,0.3,0.3,0.3),
  Max_rC=c(0.5,0.5,0.5,0.5)
  , nboot=1, seed=12
){
  list2env(datpreplist, envir = environment())
  
  # draw rCIm rCSm rCIy rCSy. Note that the rC used in 'fun_fixC_sr_Z.R' also contain rCX=rCZ=0 
  set.seed(seed = seed)
  # rC4.Kdraws=data.frame(rCIm=runif(K, 0.3, 0.5),rCSm=runif(K, 0.3, 0.5),
  #                       rCIy=runif(K, 0.3, 0.5),rCSy=runif(K, 0.3, 0.5))
  rC1.hybrid <- function( i=1 , nboot=1 ){
    rC4.Kdraws=as.matrix(rC4.Kdraws)
    rCISm=rC4.Kdraws[i, 1:2]
    rCIy=rC4.Kdraws[i, 3]
    rCSy=rC4.Kdraws[i, 4]
    # see if in the range limit ----
    Rlatent=diag(1, 1+nrow(Rc),1+ncol(Rc))
    colnames(Rlatent)=c('C',colnames(Rc))
    rownames(Rlatent)=colnames(Rlatent)
    Rlatent[1,2:5]=c(rCISm,rCIy,rCSy)
    Rlatent[2:5,1]=c(rCISm,rCIy,rCSy)
    Rlatent[-1,-1]=Rc
    ifrangelim = (det(Rlatent)>0)
    
    # transform to muC ----
    COVlatent=diag(c(1,sqrt(diag(COVc))),nrow = 1+nrow(Rc))%*%Rlatent%*%diag(c(1,sqrt(diag(COVc))),nrow = 1+nrow(Rc))
    colnames(COVlatent)=c('C',colnames(Rc))
    rownames(COVlatent)=colnames(COVlatent)
    # OLSb=solve(COVlatent[c('C','IM','SM','X','Z'),c('C','IM','SM','X','Z')])%*%COVlatent[c('C','IM','SM','X','Z'),c('IY','SY')]
    # OLSxcols = which(!colnames(COVlatent)%in%c('IY','SY'))
    OLSxcols = which(!colnames(COVlatent)%in%c('IY','SY', Z2names)) # Z2 are covariates in the models of ISM only 
    OLSb=solve(COVlatent[OLSxcols,OLSxcols])%*%COVlatent[OLSxcols,c('IY','SY')]
    # COVlatent[c('IY','SY'),c('IY','SY')] - COVlatent[c('IY','SY'),OLSxcols]%*%solve(COVlatent[OLSxcols,OLSxcols])%*%COVlatent[OLSxcols,c('IY','SY')]
    mubC=OLSb[1,]
    names(mubC)=c("mubC.Iy", "mubC.Sy")
   
    #  # COVlatent
    # COVlatent=diag(c(1,sqrt(diag(COVc))),nrow = 1+nrow(Rc))%*%Rlatent%*%diag(c(1,sqrt(diag(COVc))),nrow = 1+nrow(Rc))
    # colnames(COVlatent)=c('C',colnames(Rc))
    # rownames(COVlatent)=colnames(COVlatent)
    # OLSb=solve(COVlatent[c('C','IM','SM','X','Z'),c('C','IM','SM','X','Z')])%*%COVlatent[c('C','IM','SM','X','Z'),c('IY','SY')]
    # mubC=OLSb[1,]
    
    muaC = rCISm*sd_ISM
    mubC_IY=mubC[1]
    mubC_SY=mubC[2]
    
    comb=matrix(0, 6, ncol(V.a_b_inISM),
                dimnames = list(paste0("fixC_",colnames(V.a_b_inISM)[1:6]),
                                paste0("naive_",colnames(V.a_b_inISM))
                ) )
    # path a
    comb[c("fixC_IM.ON_X"),c("naive_IM.ON_X")]=1 
    comb[c("fixC_SM.ON_X"),c("naive_SM.ON_X")]=1
    # path b
    #Iy
    comb[c("fixC_IY.ON_IM"),c("naive_IY.ON_IM","naive_inISM11","naive_inISM12")]=
      c(1, -mubC_IY*muaC)
    comb[c("fixC_IY.ON_SM"),c("naive_IY.ON_SM","naive_inISM12","naive_inISM22")]=
      c(1, -mubC_IY*muaC)
    #Sy
    comb[c("fixC_SY.ON_IM"),c("naive_SY.ON_IM","naive_inISM11","naive_inISM12")]=
      c(1, -mubC_SY*muaC)
    comb[c("fixC_SY.ON_SM"),c("naive_SY.ON_SM","naive_inISM12","naive_inISM22")]=
      c(1, -mubC_SY*muaC)
    # est, se of paths a b ----
    fixC_E.a_b = comb %*% E.a_b_inISM
    fixC_V.a_b = comb %*% V.a_b_inISM %*% t(comb)
    fixC_SE.a_b = sqrt(diag(fixC_V.a_b))
    # fixC_Corr.a_b = diag(1/fixC_se.a_b,length(fixC_se.a_b))%*%fixC_V.a_b%*%diag(1/fixC_se.a_b,length(fixC_se.a_b))
    fixC_E.med = as_tibble( t(fixC_E.a_b) ) %>% 
      mutate( fixC_X_IM_IY = fixC_IM.ON_X * fixC_IY.ON_IM ) %>% 
      mutate( fixC_X_SM_IY = fixC_SM.ON_X * fixC_IY.ON_SM ) %>%
      mutate(fixC_X_IM_SY = fixC_IM.ON_X * fixC_SY.ON_IM ) %>%
      mutate(fixC_X_SM_SY = fixC_SM.ON_X * fixC_SY.ON_SM ) %>%
      # total
      mutate(fixC_X_TO_IY = fixC_X_IM_IY + fixC_X_SM_IY) %>%
      mutate(fixC_X_TO_SY = fixC_X_IM_SY + fixC_X_SM_SY) %>% 
      gather(key = "paramName","est")
    
    # # mcci ----
    # nbootfixC_a_b = rmvnorm(nboot, mean = fixC_E.a_b, sigma = fixC_V.a_b)
    # colnames(nbootfixC_a_b)=colnames(fixC_V.a_b)
    # 
    # nbootfixC_params = as_tibble(nbootfixC_a_b) %>% 
    #   mutate( fixC_X_IM_IY = fixC_IM.ON_X * fixC_IY.ON_IM ) %>% 
    #   mutate( fixC_X_SM_IY = fixC_SM.ON_X * fixC_IY.ON_SM ) %>%
    #   mutate(fixC_X_IM_SY = fixC_IM.ON_X * fixC_SY.ON_IM ) %>%
    #   mutate(fixC_X_SM_SY = fixC_SM.ON_X * fixC_SY.ON_SM ) %>%
    #   # total
    #   mutate(fixC_X_TO_IY = fixC_X_IM_IY + fixC_X_SM_IY) %>%
    #   mutate(fixC_X_TO_SY = fixC_X_IM_SY + fixC_X_SM_SY) 
    # fixC_CIlo = nbootfixC_params %>%
    #   summarise_all( quantile, alph/2 ) %>% gather(key = "paramName","mcci.lo")
    # fixC_CIup = nbootfixC_params %>%
    #   summarise_all( quantile, 1-alph/2 ) %>% gather(key = "paramName","mcci.up")
    # 
    # fixC_MCCI = full_join(fixC_CIlo, fixC_CIup, by='paramName') %>%
    #   mutate(mcci.ifsig = (mcci.lo*mcci.up)>0 )
    
    # the modified hybrid approach does not obtain test result for each of the K sets of sensi param; rather, it use the K sets of random samples from the adjusted sampling distribution for inference (e.g., 95% CI)
    ##3[ii]draw a random sample from the adjusted sampling distribution
    nboot=1
    onebootfixC_a_b = rmvnorm(nboot, mean = fixC_E.a_b, sigma = fixC_V.a_b)
    colnames(onebootfixC_a_b)=colnames(fixC_V.a_b)
    
    onebootfixC_params = as_tibble(onebootfixC_a_b) %>% 
      mutate( fixC_X_IM_IY = fixC_IM.ON_X * fixC_IY.ON_IM ) %>% 
      mutate( fixC_X_SM_IY = fixC_SM.ON_X * fixC_IY.ON_SM ) %>%
      mutate(fixC_X_IM_SY = fixC_IM.ON_X * fixC_SY.ON_IM ) %>%
      mutate(fixC_X_SM_SY = fixC_SM.ON_X * fixC_SY.ON_SM ) %>%
      # total
      mutate(fixC_X_TO_IY = fixC_X_IM_IY + fixC_X_SM_IY) %>%
      mutate(fixC_X_TO_SY = fixC_X_IM_SY + fixC_X_SM_SY) 
    
    # # gather output
    # onesr_hybrid.a_b = data.frame(paramName = rownames(fixC_V.a_b), est = fixC_E.a_b, se = fixC_SE.a_b, ifrangelim = ifrangelim)
    # onesr_hybrid.med = full_join(fixC_E.med[7:12, ], fixC_MCCI[7:12, ], by = 'paramName') %>%
    #   mutate( ifrangelim = ifrangelim )
    # 
    # onesr_hybrid.out = list( onesr_hybrid.a_b, onesr_hybrid.med )
    # names(onesr_hybrid.out) = c( 'a_b', 'med' )
    
    # gather output
    onesr_hybrid.out=onebootfixC_params %>% mutate( ifrangelim = ifrangelim )
    
    return(onesr_hybrid.out)
  }
  
  # set.seed(seed = 12)
  # draw rCIm rCSm rCIy rCSy. Note that the rC used in 'fun_fixC_sr_Z.R' also contain rCX=rCZ=0 
  # K=1000
  # Min_rC=c(0.3,0.3,0.3,0.3)
  # Max_rC=c(0.5,0.5,0.5,0.5)
  rC4.Kdraws=data.frame(rCIm=runif(K, Min_rC[1], Max_rC[1]),rCSm=runif(K, Min_rC[2], Max_rC[2]),
                        rCIy=runif(K, Min_rC[3], Max_rC[3]),rCSy=runif(K, Min_rC[4], Max_rC[4]))
  hybridres = lapply(1:K, rC1.hybrid, nboot=1)
  
  mhybrid_res=data.frame(do.call(rbind, hybridres))
  
  out=mget(ls(), envir = environment())
  return(out)
}


fun.summhybridC<-function(hybridreslist){
  list2env(hybridreslist, envir = environment())
  
  prop_rangelim =mean(mhybrid_res$ifrangelim)
  mhybrid_res1=mhybrid_res[ mhybrid_res$ifrangelim==1, 1:12] %>% 
    pivot_longer(1:12, names_to = 'param') %>% 
    mutate( param=str_replace(param,'fixC_','') )
  
  abmedsensi_raw = mhybrid_res1 %>% group_by(param) %>%
    summarise_all( list(mean= ~mean(.), sd=~sd(.), ci.lo=~quantile(. ,0.025), ci.up=~quantile(. ,0.975)) )
  
  abmedsensi = abmedsensi_raw[c(1,4,2,5,3,6,7:10) , ]
  
  abmednaive_raw = sample_naive.Parameters[
    c((sample_naive.Parameters$paramHeader=="IM.ON"&sample_naive.Parameters$param=='X') |(sample_naive.Parameters$paramHeader=="SM.ON"&sample_naive.Parameters$param=='X') |
        (sample_naive.Parameters$paramHeader=="IY.ON"&sample_naive.Parameters$param=='IM') | (sample_naive.Parameters$paramHeader=="IY.ON"&sample_naive.Parameters$param=='SM') | 
        (sample_naive.Parameters$paramHeader=="SY.ON"&sample_naive.Parameters$param=='IM') | (sample_naive.Parameters$paramHeader=="SY.ON"&sample_naive.Parameters$param=='SM') |
        (sample_naive.Parameters$paramHeader=="New.Additional.Parameters")&(sample_naive.Parameters$param%in%c('X_IM_IY','X_SM_IY','X_IM_SY','X_SM_SY')) ) , ]
  
  abmednaive = abmednaive_raw[ c(5,6,1,3,2,4,7,9,8,10) ,-1 ]
  abmednaive$param[1:6] = abmedsensi$param[1:6]
  MCSAout = cbind(abmednaive, abmedsensi[,-1], prop_rangelim)
  colnames(MCSAout)=c("param", paste0("ML_Original_",colnames(abmednaive[,-1]) ), paste0("MCSA_",colnames(abmedsensi[,-1]) ), "Proportion_of_draws_within_admissible_ranges" )
  
  
  out=mget(ls(), envir = environment())
  return(out)
}

# aaa=fun.hybridrC_sr(datpreplist = datprep_list,K=4)
# aaasumm=fun.summhybridC( hybridreslist = aaa )


