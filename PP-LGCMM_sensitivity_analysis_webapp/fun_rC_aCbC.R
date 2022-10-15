## transformation between rC (zero-order confounder correlations) and aCbC (confounder path coeffs)


# transform rC to aC,bC ----

rC.to.aCbC <- function(datpreplist, 
                        rC_ImSmIySy=c(.1,.1,.1,.1) ){ # Z2 are covariates in the models of ISM only 
  list2env(datpreplist, envir = environment())
  if ( length(mplusOut$results$errors)  == 0 ) {
    rC = c(rC_ImSmIySy, 0, rep(0, length(Z1cols)+length(Z2cols)+length(Z3cols) ))
    # Rlatent
    Rlatent=diag(1, 1+nrow(Rc),1+ncol(Rc))
    Rlatent[1, -1]=rC
    Rlatent[-1,1]=rC
    Rlatent[-1,-1]=Rc
    colnames(Rlatent)=c('C',colnames(Rc))
    rownames(Rlatent)=colnames(Rlatent)
    # det(Rlatent)
    ifrangelim = (det(Rlatent)>0)
    
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
    names(mubC)=c("bIY_C", "bSY_C")
    
    ## rCISm to aC
    rCISm = rC_ImSmIySy[1:2]
    muaC = rCISm*sd_ISM
    
    aCbC = c(muaC, mubC)
    names(aCbC) = c("aIM_C",'aSM_C',"bIY_C", "bSY_C")
    out=data.frame(t(c(aCbC, ifrangelim)))
    colnames(out) = c("aIM_C",'aSM_C',"bIY_C", "bSY_C","if_within_admissible_ranges")
  }
  if ( length(mplusOut$results$errors)  > 0 ){
    out=data.frame( error_originalmodel=c('Mplus fitting original model M0 returns errors and/or warnings.',
                          unlist(mplusOut$results$errors), unlist(mplusOut$results$warnings) ) 
                    )
  }
  
  return(out)
}


# transform aC bC to rC ----

aCbC.to.rC <- function( datpreplist, aCbC=c(0.14, 0.14, 0.14, 0.14) ) {
  list2env(datpreplist, envir = environment())
  
  ## bC.to.rC
  mubC_IY=aCbC[3]; mubC_SY=aCbC[4]; muaC=aCbC[1:2]
  
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
  rCSy = (t(est_b[2, ]) %*% muaC + mubC_SY) / (sd_ISY[2])
  
  ## aC to rC
  rCISm = aCbC[1:2]/sd_ISM

  rC_ImSmIySy = c(rCISm, rCIy, rCSy)
  names(rC_ImSmIySy) = c('rIM_C', 'rSM_C','rIY_C', 'rSY_C')
  
  # if range limits
  rC = c(rC_ImSmIySy, 0, rep(0, length(Z1cols)+length(Z2cols)+length(Z3cols) ))
  # Rlatent
  Rlatent=diag(1, 1+nrow(Rc),1+ncol(Rc))
  Rlatent[1, -1]=rC
  Rlatent[-1,1]=rC
  Rlatent[-1,-1]=Rc
  colnames(Rlatent)=c('C',colnames(Rc))
  rownames(Rlatent)=colnames(Rlatent)
  # det(Rlatent)
  ifrangelim = (det(Rlatent)>0)
  
  out= data.frame(t(c(rC_ImSmIySy, ifrangelim)))
  colnames(out) = c('rIM_C', 'rSM_C','rIY_C', 'rSY_C','if_within_admissible_ranges')
  
  if ( length(mplusOut$results$errors)  > 0 ){
    out=data.frame( error_originalmodel=c('Mplus fitting original model M0 returns errors and/or warnings.',
                                          unlist(mplusOut$results$errors), unlist(mplusOut$results$warnings) ) 
    )
  }
  
  return(out)
}




