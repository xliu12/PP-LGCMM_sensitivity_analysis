

# including Z: only solverC.range() is changed # rC's order: C with c('IM','SM','IY','SY','X','Z') #rC in the same order as Rc!!!

solve.quadeq <- function(A=1,B=2,C=1) { # Ax^2+Bx+C=0
  if(A==0) {
    s1=s2= -C/B
  }
  if(A!=0) {
    delta=B^2-4*A*C
    s1=( -B + sqrt(delta) ) / (2*A) # sqrt(-1)=NaN
    s2=( -B - sqrt(delta) ) / (2*A)
  }
  
  return( c(s1,s2) )
}
# solve range limits of rC ----
# rC's order: C with c('IM','SM','IY','SY','X','Z')


#rC in the same order as Rc!!!
solverC.range <- function( whichrC=3,R.c=Rc, rC ) { 
  # solved from rC'*(Rc^-1)*rC <=1
  invRc = solve(R.c) # require R_C posi-def
  A = invRc[whichrC, whichrC] 
  B = 2 * invRc[whichrC, -whichrC] %*% rC[-whichrC]
  C = -1 + t(rC[-whichrC]) %*% invRc[-whichrC, -whichrC] %*% rC[-whichrC]
  # if( B^2-4*A*C < 0 ): NaNs for range; non-convergence
  s12.rCrange = solve.quadeq(A=A,B=B,C=C)
  names(s12.rCrange)=c("s1.rCrange","s2.rCrange")
  s12.rCrange=c(s12.rCrange)
  return(s12.rCrange)
}  

# solve bC needed given muaC ----
# if !is.nan(range) but bsig bzero solutions outside range limits: 
#convergence but warning for non-positive definite latent variable covariance matrix 
# if is.nan(range): non-convergence; so the solutions are meaningless
solvebC.bsig <- function( critA_bpath, muaC=c(0,0) ) { 
  # critA_bpath: 3by3, bpath, inISM**,inISM**   
  # e.g., solve mubC_IY for (est/se)^2 = qchi #xEEx'-qchi*xVx'=x(EE-qchi*V)x'=0. let critA=EE-qchi*V
  
  A = t(muaC) %*% critA_bpath[-1,-1] %*% muaC
  B = -2*t(muaC) %*% critA_bpath[-1,1]
  C = critA_bpath[1,1]
  s12.bsig = solve.quadeq(A=A,B=B,C=C)
  names(s12.bsig)=c("s1.bsig","s2.bsig")
  s12.bsig=c(s12.bsig)
  return(s12.bsig)
}  

solvebC.bzero <- function( estE_bpath, muaC=c(0,0) ) {
  s.bzero = estE_bpath[1] / ( t(muaC) %*% estE_bpath[-1] ) # 1/0=Inf 
  s.bzero=c(s.bzero)
  if( is.infinite(s.bzero) ) {
    s.bzero = NA
  }
  names(s.bzero)="s.bzero"
  
  return(s.bzero)
}



# transform rC to aC,bC ----
# rCISm.to.aC <- function(COVlatent){
#   OLSa=solve(COVlatent[c('C','X','Z'),c('C','X','Z')])%*%COVlatent[c('C','X','Z'),c('IM','SM')]
#   muaC=OLSa[1,]
#   names(muaC)=c("muaC.Im", "muaC.Sm")
#   return(muaC)
# }
rCISy.to.bC <- function(rC, Rc, COVc, Z2names){ # Z2 are covariates in the models of ISM only 
  # Rlatent
  Rlatent=diag(1, 1+nrow(Rc),1+ncol(Rc))
  Rlatent[1, -1]=rC
  Rlatent[-1,1]=rC
  Rlatent[-1,-1]=Rc
  colnames(Rlatent)=c('C',colnames(Rc))
  rownames(Rlatent)=colnames(Rlatent)
  # det(Rlatent)
  # COVlatent
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
  return(mubC)
}
rCISm.to.aC <- function(rCISm, sd.ISM=sd_ISM ){
  muaC = rCISm*sd.ISM
  return(muaC)
}
# rCISy.to.bC <- function(rCISy=c(0.2,0.2), sd_ISY,
#                         est_b, rCISm=c(0.2,0.2), sd_ISM){
#   mubC = (rCISy*sd_ISY) - est_b %*% (rCISm*sd_ISM)
#   return(mubC)
# }

# transform aC bC to rC ----
aC.to.rCISm <- function(muaC=c(0,0), sd_ISM){
  rCISm = muaC/sd_ISM
  return(rCISm)
}
bC_IY.to.rCIy <- function(mubC_IY, mubC_SY, muaC=c(0,0), V.a_b_inISM, E.a_b_inISM, sd_ISY ) {
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
  # est, se of paths a b 
  fixC_E.a_b = comb %*% E.a_b_inISM
  est_b=matrix(fixC_E.a_b[3:6], 2,2,byrow = T)
  rownames(est_b)=c('IY.ON','SY.ON')
  colnames(est_b)=c('_IM','_SM')
  rCIy = (t(est_b[1, ]) %*% muaC + mubC_IY) / (sd_ISY[1])
  return(rCIy)
}
bC_SY.to.rCSy <- function(mubC_SY, mubC_IY, muaC=c(0,0), V.a_b_inISM ,E.a_b_inISM, sd_ISY) {
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
  fixC_E.a_b = comb %*% E.a_b_inISM
  est_b=matrix(fixC_E.a_b[3:6], 2,2,byrow = T)
  rownames(est_b)=c('IY.ON','SY.ON')
  colnames(est_b)=c('_IM','_SM')
  rCSy = (t(est_b[2, ]) %*% muaC + mubC_SY) / (sd_ISY[2])
  return(rCSy)
}
# bC1.to.rC1 <- function(bC1=0, estb.onISM, sd1, muaC=c(0,0) ) {
#   rC1 = (t(estb.onISM) %*% muaC + bC1) / (sd1)
#   return(rC1)
# }


# sensitivity analysis ----
solve.sensi = function(rC=c(rep(0.2,4),0, 0), envlist ) {
  list2env(envlist, envir = environment())
  
  rC = as.numeric(rC) # rC col-vec
  
  # range limits rCIy rCSy ----
  range12_rCIy = solverC.range(whichrC = 3,R.c=Rc, rC=rC ) 
  range.lo_rCIy = min(range12_rCIy)
  range.up_rCIy = max(range12_rCIy)
  range12_rCSy = solverC.range(whichrC = 4, R.c=Rc, rC=rC ) 
  range.lo_rCSy = min(range12_rCSy)
  range.up_rCSy = max(range12_rCSy)
  
  # # Rlatent
  # Rlatent=diag(1, 1+nrow(Rc),1+ncol(Rc))
  # Rlatent[1, -1]=rC
  # Rlatent[-1,1]=rC
  # Rlatent[-1,-1]=Rc
  # colnames(Rlatent)=c('C',colnames(Rc))
  # rownames(Rlatent)=colnames(Rlatent)
  # # COVlatent
  # COVlatent=diag(c(1,sqrt(diag(COVc))),nrow = 1+nrow(Rc))%*%Rlatent%*%diag(c(1,sqrt(diag(COVc))),nrow = 1+nrow(Rc))
  # colnames(COVlatent)=c('C',colnames(Rc))
  # rownames(COVlatent)=colnames(COVlatent)

  # transform to muaC -----
  # muaC = rCISm.to.aC(COVlatent)
  muaC = rCISm.to.aC(rCISm = rC[1:2], sd.ISM = sd_ISM)
  # mubC = rCISy.to.bC(rCISy = rC[3:4],sd_ISY,est_b,rCISm = rC[1:2],sd_ISM)
  mubC = rCISy.to.bC(rC=rC, Rc=Rc, COVc = COVc, Z2names = Z2names)
  # bsig: solve bC given muaC ----
  bC_IYs12_sigb.IyonIm= solvebC.bsig(
    critA_bpath = critA.IY.ON_IM_fixC, muaC = muaC) # bC_IY influence IY on IM SM
  bC_IYs12_sigb.IyonSm= solvebC.bsig(
    critA_bpath = critA.IY.ON_SM_fixC, muaC = muaC)
  bC_SYs12_sigb.SyonIm= solvebC.bsig(
    critA_bpath = critA.SY.ON_IM_fixC, muaC = muaC) # bC_SY influence SY on IM SM
  bC_SYs12_sigb.SyonSm= solvebC.bsig(
    critA_bpath = critA.SY.ON_SM_fixC, muaC = muaC)
  
  # bzero: solve bC given muaC ----
  bC_IY_zerob.IyonIm= solvebC.bzero(
    estE_bpath = estE.IY.ON_IM_fixC, muaC = muaC) # bC_IY influence IY on IM SM
  bC_IY_zerob.IyonSm= solvebC.bzero(
    estE_bpath = estE.IY.ON_SM_fixC, muaC = muaC)
  bC_SY_zerob.SyonIm= solvebC.bzero(
    estE_bpath = estE.SY.ON_IM_fixC, muaC = muaC) # bC_SY influence SY on IM SM
  bC_SY_zerob.SyonSm= solvebC.bzero(
    estE_bpath = estE.SY.ON_SM_fixC, muaC = muaC)
  
  # transform bC_IY bC_SY to rCIy rCSy ----
  
  # rCIy
  # s12rCIy_sigb.IyonIm =
  #   sapply(bC_IYs12_sigb.IyonIm[1:2], FUN = bC1.to.rC1,
  #          estb.onISM = est_b[1,], sd1 = sd_ISY[1],muaC = muaC)
  s12rCIy_sigb.IyonIm =
    sapply(bC_IYs12_sigb.IyonIm[1:2], FUN = bC_IY.to.rCIy, mubC_SY=mubC[2],muaC = muaC, V.a_b_inISM=V.a_b_inISM,  E.a_b_inISM=E.a_b_inISM, sd_ISY=sd_ISY  )
  # s12rCIy_sigb.IyonIm = s12rCIy_sigb.IyonIm[ 
  #   (s12rCIy_sigb.IyonIm < 1) & (s12rCIy_sigb.IyonIm > -1)
  #   ]
  bound.lo_rCIy_sigb.IyonIm = ifelse(length(s12rCIy_sigb.IyonIm)==0, 
                                     NA, min(s12rCIy_sigb.IyonIm) )
  bound.up_rCIy_sigb.IyonIm = ifelse(length(s12rCIy_sigb.IyonIm)==0, 
                                     NA, max(s12rCIy_sigb.IyonIm) )
  # rCIy_zerob.IyonIm = bC1.to.rC1(
  #   bC1 = bC_IY_zerob.IyonIm,
  #   estb.onISM = est_b[1,], sd1 = sd_ISY[1],muaC = muaC)
  rCIy_zerob.IyonIm = bC_IY.to.rCIy(mubC_IY = bC_IY_zerob.IyonIm, mubC_SY = mubC[2], muaC = muaC, V.a_b_inISM=V.a_b_inISM,  E.a_b_inISM=E.a_b_inISM, sd_ISY=sd_ISY )
  # rCIy_zerob.IyonIm = rCIy_zerob.IyonIm [
  #   (rCIy_zerob.IyonIm < 1) & (rCIy_zerob.IyonIm > -1)
  # ]
  rCIy_zerob.IyonIm = ifelse(length(rCIy_zerob.IyonIm)==0,NA,rCIy_zerob.IyonIm)
  
  # s12rCIy_sigb.IyonSm =
  #   sapply(bC_IYs12_sigb.IyonSm[1:2], FUN = bC1.to.rC1,
  #          estb.onISM = est_b[1,], sd1 = sd_ISY[1],muaC = muaC)
  s12rCIy_sigb.IyonSm =
    sapply(bC_IYs12_sigb.IyonSm[1:2], FUN = bC_IY.to.rCIy,
           mubC_SY = mubC[2],muaC = muaC, V.a_b_inISM=V.a_b_inISM,  E.a_b_inISM=E.a_b_inISM , sd_ISY=sd_ISY)
  # s12rCIy_sigb.IyonSm = s12rCIy_sigb.IyonSm[
  #   (s12rCIy_sigb.IyonSm < 1) & (s12rCIy_sigb.IyonSm > -1)
  # ]
  bound.lo_rCIy_sigb.IyonSm = ifelse(length(s12rCIy_sigb.IyonSm)==0,
                                     NA, min(s12rCIy_sigb.IyonSm) )
  bound.up_rCIy_sigb.IyonSm = ifelse(length(s12rCIy_sigb.IyonSm)==0,
                                     NA, max(s12rCIy_sigb.IyonSm) )
  # rCIy_zerob.IyonSm = bC1.to.rC1(
  #   bC1 = bC_IY_zerob.IyonSm,
  #   estb.onISM = est_b[1,], sd1 = sd_ISY[1],muaC = muaC)
  rCIy_zerob.IyonSm = bC_IY.to.rCIy(
    mubC_IY = bC_IY_zerob.IyonSm, mubC_SY = mubC[2],muaC = muaC, V.a_b_inISM=V.a_b_inISM,  E.a_b_inISM=E.a_b_inISM , sd_ISY=sd_ISY)
  # rCIy_zerob.IyonSm = rCIy_zerob.IyonSm [
  #   (rCIy_zerob.IyonSm < 1) & (rCIy_zerob.IyonSm > -1)
  #   ]
  rCIy_zerob.IyonSm = ifelse(length(rCIy_zerob.IyonSm)==0,NA,rCIy_zerob.IyonSm)
  
  # # rCSy
  # s12rCSy_sigb.SyonIm =
  #   sapply(bC_SYs12_sigb.SyonIm[1:2], FUN = bC1.to.rC1,
  #          estb.onISM = est_b[2,], sd1 = sd_ISY[2],muaC = muaC)
  s12rCSy_sigb.SyonIm =
    sapply(bC_SYs12_sigb.SyonIm[1:2], FUN = bC_SY.to.rCSy,
           mubC_IY=mubC[1],muaC = muaC, V.a_b_inISM=V.a_b_inISM,  E.a_b_inISM=E.a_b_inISM, sd_ISY=sd_ISY )
  # s12rCSy_sigb.SyonIm = s12rCSy_sigb.SyonIm[
  #   (s12rCSy_sigb.SyonIm < 1) & (s12rCSy_sigb.SyonIm > -1)
  # ]
  bound.lo_rCSy_sigb.SyonIm = ifelse(length(s12rCSy_sigb.SyonIm)==0,
                                     NA, min(s12rCSy_sigb.SyonIm))
  bound.up_rCSy_sigb.SyonIm = ifelse(length(s12rCSy_sigb.SyonIm)==0,
                                     NA, max(s12rCSy_sigb.SyonIm))
  # rCSy_zerob.SyonIm = bC1.to.rC1(
  #   bC1 = bC_SY_zerob.SyonIm,
  #   estb.onISM = est_b[2,], sd1 = sd_ISY[2],muaC = muaC )
  rCSy_zerob.SyonIm = bC_SY.to.rCSy(
    mubC_SY = bC_SY_zerob.SyonIm,mubC_IY = mubC[1], muaC = muaC, V.a_b_inISM=V.a_b_inISM ,  E.a_b_inISM=E.a_b_inISM, sd_ISY=sd_ISY )
  # rCSy_zerob.SyonIm = rCSy_zerob.SyonIm[
  #   (rCSy_zerob.SyonIm < 1) & (rCSy_zerob.SyonIm > -1)
  # ]
  rCSy_zerob.SyonIm = ifelse(length(rCSy_zerob.SyonIm)==0,NA,rCSy_zerob.SyonIm)
  
  # s12rCSy_sigb.SyonSm =
  #   sapply(bC_SYs12_sigb.SyonSm[1:2], FUN = bC1.to.rC1,
  #          estb.onISM = est_b[2,], sd1 = sd_ISY[2],muaC = muaC)
  s12rCSy_sigb.SyonSm =
    sapply(bC_SYs12_sigb.SyonSm[1:2], FUN = bC_SY.to.rCSy,
           mubC_IY = mubC[1], muaC = muaC, V.a_b_inISM=V.a_b_inISM,  E.a_b_inISM=E.a_b_inISM , sd_ISY=sd_ISY)
  # s12rCSy_sigb.SyonSm = s12rCSy_sigb.SyonSm[
  #   (s12rCSy_sigb.SyonSm < 1) & (s12rCSy_sigb.SyonSm > -1)
  # ]
  bound.lo_rCSy_sigb.SyonSm = ifelse(length(s12rCSy_sigb.SyonSm)==0,
                                     NA, min(s12rCSy_sigb.SyonSm) )
  bound.up_rCSy_sigb.SyonSm = ifelse(length(s12rCSy_sigb.SyonSm)==0,
                                     NA, max(s12rCSy_sigb.SyonSm) )
  # rCSy_zerob.SyonSm = bC1.to.rC1(
  #   bC1 = bC_SY_zerob.SyonSm,
  #   estb.onISM = est_b[2,], sd1 = sd_ISY[2],muaC = muaC )
  rCSy_zerob.SyonSm = bC_SY.to.rCSy(
    mubC_SY = bC_SY_zerob.SyonSm, mubC_IY = mubC[1],muaC = muaC, V.a_b_inISM=V.a_b_inISM,  E.a_b_inISM=E.a_b_inISM , sd_ISY=sd_ISY )
  # rCSy_zerob.SyonSm = rCSy_zerob.SyonSm[
  #   (rCSy_zerob.SyonSm < 1) & (rCSy_zerob.SyonSm > -1)
  # ]
  rCSy_zerob.SyonSm = ifelse(length(rCSy_zerob.SyonSm)==0,NA,rCSy_zerob.SyonSm)
  # gather results ----
  lsnames = ls(envir = environment())
  outnames = c(grep("mu",lsnames,value = T),
    grep('bound\\.',lsnames,value = T),
    grep('^rC[A-Za-z]*\\_zerob',lsnames,value = T),
    grep('range\\.',lsnames,value = T))
  
  out = mget(outnames, envir = environment())
  
  return(out)
}



fixC_s.SensiRes <- function(all_rC, envlist) { # all_rC: # in the same order as Rc!!!: rCIm,rCSm,rCIy,rCSy,X
  # list2env(envlist, envir = environment())
  SensiRes = NULL
  for(i in 1:nrow(all_rC)){
    rC=all_rC[i, ]
    outi = solve.sensi(rC = rC, envlist=envlist )
    SensiRes = rbind(SensiRes, unlist(c(all_rC[i, ], outi)) )
  }
  
  SensiRes[is.infinite(SensiRes)] = NA
  SensiRes1 = as_tibble(SensiRes)
  return(SensiRes1)
}






