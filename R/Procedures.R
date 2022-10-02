### Main procedures

#' Global SSR minimizers procedure
#'
#' Function obtains global minimizers using the recursive algorithm to compute
#' and minimize SSR over all possible segments. The procedure is required to
#' conduct supF, UDMax, WDMax and supF(l+1|l) test
#'
#'@aliases doglob
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@return A list containing the following components:
#'\itemize{
#'\item{glb} {Minimum global SSR}
#'\item{datevec} {Vector of dates (optimal minimizers)}
#'\item{bigvec} {Associated SSRs}}
#'@export
#'
doglob = function (y,z,x,m,eps,eps1,maxi,fixb,betaini,printd){

#check if model is pure or partial change
if (is.null(x)) {p = 0}
else {p = dim(x)[2]}

q = dim(z)[2]
bigT = dim(y)[1]
h = round(eps1*bigT)


if(p == 0) {
  if (printd == 1){
  print('This is a pure structural change model with the following specifications:')
  print(paste(q,'regressors z with allowed to change coefficients'))
  print(paste('maximum number of breaks:',m)) }
  out = dating(y,z,h,m,q,bigT)
  }
else{
  if (printd == 1){
  print('This is a partial structural change model with the following specifications:')
  print(paste(q,'regressors z with allowed to change coefficients'))
  print(paste(p,'regressors x with fixed coefficients'))
  print(paste('maximum number of breaks:',m))
  print(paste('initial values of β option (1)=TRUE/(0)=FALSE',fixb))
  if(fixb == 1) {print('initial values β')
    print(betaini)}
  print(paste('convergence criterion:',eps))
  print(paste('print iteration option (1)=TRUE/(0)=FALSE',printd))}
  out = nldat(y,z,x,h,m,p,q,bigT,fixb,eps,maxi,betaini,printd)
}
glb = out$glb
datevec = out$datevec
bigvec = out$bigvec
#printing results
if (printd == 1) {
for (i in 1:m){
  print(paste('Model with',i,'breaks has SSR:',glb[i,1]))
  print('The dates of breaks are:')
  print(datevec[1:i,i])
}
}


date = out$datevec


  return (out)
}

#' SupF, UDMax & WDMax testing procedure
#'@import 
#' The procedure calculate the test statistics and print results of the 2 main tests:
#' \itemize{
#' \item{SupF test} {F test of 0 vs m breaks}
#' \item{Double Max test} {UDMax: the unweighted version
#' and WDMax: the weighted version}
#'}
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param prewhit option to use AR(1) for prewhitening
#'@param robust,hetdat,hetvar options on error terms assumptions
#'@return A list that contains following:
#'\itemize{
#'\item{ftest}{SupF test statistics}
#'\item{cv}{Critical values for Sup F test }
#'\item{wftest}{Weighted SupF test}
#'\item{Dmax}{Double Max test statistics}
#'\item{cv_Dmax}{Critical values for Double Max test}
#'\item{WDmax}{WDmax test statistics}
#'
#'}
#'@export
#'
#'
dotest = function(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,robust,
                  hetdat,hetvar){
  siglev=matrix(c(10,5,2.5,1),4,1)
  #printd = 0 #suppress output printing
  #print('Output from testing procedure')
  out = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
  datevec = out$datevec

  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  bigT = dim(y)[1]
  h = round(eps1*bigT)
  #procedure for F test
  #print('a) supF tests against a fixed number of breaks')

  ftest = matrix(0L, nrow = m, ncol = 1)
  wftest = matrix(0L, nrow = m, ncol = 1)

  for (i in 1:m){
    ftest[i,1] = pftest(y,z,i,q,bigT,datevec,prewhit,robust,x,p,hetdat,hetvar)
    if(printd==1){
    print(paste('supF test for 0 versus',i,'breaks (scaled by q):',ftest[i,1]))}
  }

  cv_supF = matrix(0L,4,m)
  for (c in 1:4){
    #critical values for supF test
    cv = getcv1(c,eps1)
    cv_supF[c,] = cv[q,1:m,drop=FALSE]
    if (printd==1){
    print(paste('The critical values at the',siglev[c,1],'% level are (for k = 1 to',m,'):'))
    print(cv[q,1:m,drop=FALSE])}
  }

  #procedure for Dmax and UDmax test

  #print('b) Dmax test against an unknown number of breaks')
  #print(paste('The UDmax test is:',max(ftest)))
  cv_Dmax = matrix(0L,4,1)
  for (c in 1:4) {
    #critical values for Dmax test
    cvm = getdmax(c,eps1)
    cv_Dmax[c,1] = cvm[q,1]
    if(printd==1){
    print(paste('The critical values at the',siglev[c,1],'% level is:',
                cvm[q,1]))}
  }


  for (c in 1:4){
    #computation of WDmax test
    cv = getcv1(c,eps1)
    for( i in 1:m){
      wftest[i,1] = cv[q,1] * ftest[i,1] / cv[q,i]
    }
    if (printd==1){
    print(paste('WDmax test at the',siglev[c,1],'% level is:',max(wftest)))}
  }
  rownames(cv_supF) = siglev
  rownames(cv_Dmax) = siglev
  
 
  out = list('ftest' = ftest, 'cv_supF' = cv_supF,
             'cv_Dmax' = cv_Dmax, 'UDmax' = max(ftest))
  out$mbreak = m
  class(out) = 'Wtest'
  
  out = compile.Wtest(out)
  return(out)
}



#' SupF(l+1|l) test
#'
#' Function computes the procedure of SupF(l+1|l) test. The function returns
#' the test statistics of supF(l+1|l) test
#' with null hypothesis is maximum number of break is l
#' and alternative hypothesis is l+1.
#' The l breaks under the null hypothesis are taken from the
#' global minimization. Also, new date (if available) and critical values based on
#' significant levels are returned for plotting and inference
#'
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param prewhit option to use AR(1) for prewhitening
#'@param robust,hetdat,hetvar options on error terms assumptions
#' @return A list that contains following:
#' \itemize{
#'\item{supfl}{SupF(l+1|l) test statistics}
#'\item{cv}{Critical values for SupF(l+1|l) test}
#'\item{ndat}{New date (if available)} }
#'
#'@export

dospflp1 = function(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
                    robust,hetdat,hetvar) {
  siglev=matrix(c(10,5,2.5,1),4,1)
  #print('supF(l+1|l) tests using global optimizers under the null')
  out = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
  datevec = out$datevec
  bigvec = out$bigvec


  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  bigT = dim(y)[1]
  h = round(eps1*bigT)

  supfl = matrix (0L,m-1,1)
  ndat = matrix (0L,m-1,1)
  for (i in seq(1,m-1,1)){
    out1 = spflp1(bigvec,datevec[1:i,i,drop=FALSE],i+1,y,z,h,q,prewhit,robust,x,p,hetdat,hetvar)
    supfl[i,1] = out1$maxf
    ndat[i,1] = out1$newd
    #print(paste('The supF(',i+1,'|',i,') test is',supfl[i,1]))
    #print(paste('It corresponds to a new break at:',ndat[i,1]))
  }

  cv_supFl = matrix(0L,4,m)

  for (c in 1:4){
    #critical values for supF(l+1|l) test
    cv = getcv2(c,eps1)
    cv_supFl[c,] = cv[q,1:m,drop=FALSE]
    #print(paste('The critical values at the',siglev[c,1],'% level are (for k = 1 to',m,'):'))
    #print(cv[q,1:m,drop=FALSE])
  }
  rownames(cv_supFl) = siglev

  
  out = list('supfl' = supfl, 'ndat' = ndat, 'cv' = cv_supFl)
  out$mbreak = m
  class(out) = 'supfltest'
  out = compile.supfltest(out)
  print(out)
  return(out)
  
}

#' Order estimation procedure
#'
#' The function carry out the procedure to estimate order
#'  using BIC and the criterion of Liu, Wu and Zidek
#'
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param bic indicator which criterion is used in selecting number of breaks
#'@return A list that contains following:
#'\item{mBIC}{number of breaks selected by BIC}
#'\item{mLWZ}{number of breaks selected by LWZ}
#'@export
#'@references
#
doorder = function(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,bic_opt) {


  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  bigT = dim(y)[1]
  h = round(eps1*bigT)

  if (p == 0){zz = z}
  else{zz = cbind(z,x)}
  temp = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
  glb = temp$glb
  bigT = dim(y)[1]
  ssr0 = nssr(y,zz)
  delta0 = 0.1 #optimal parameters in LWZ paper
  c0 = 0.299
  glob= matrix(0L, nrow = m+1, ncol=1)
  glob[1,1] = ssr0
  glob[seq(2,m+1),1] = glb

  bic = matrix(0L,nrow = m+1, ncol = 1)
  lwz = matrix(0L,nrow = m+1, ncol = 1)

  for (i in seq(1,m+1)){
    bic [i,1] = log(glob[i,1]/bigT) + log(bigT)*(i-1)*(q+1)/bigT
    lwz[i,1] = log(glob[i,1]/(bigT-i*q-i+1)) +
      ((i-1)*(q+1)*c0*(log(bigT))^(2+delta0))/bigT

   
  }

  mBIC = which.min(bic) - 1
  mLWZ = which.min(lwz) - 1
  
  if (bic_opt == 1){
    m=mBIC
    p_name = 'BIC'
  }else{
    m=mLWZ
    p_name = 'LWZ'
  }
  date = temp$datevec[seq(1,m,1),m,drop=FALSE]
  if (m == 0){
    print('There are no breaks in this model and estimation is skipped')
    return(NULL)
  }
  else{
    hetq=0
    hetomega=0
    hetdat=1
    hetomega=1
    hetvar=1
    robust=1
    prewhit=1
    out = estim(m,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = p_name
    out$nbreak = m
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out = compile.model(out)
    return(out)
  }
  
}


#sequential procedure
sequa = function(m,signif,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1){


  dv = matrix(0L, nrow = m+2, ncol = 1)
  dv2 = matrix(0L, nrow = m+2, ncol = 1)
  ftestv = matrix(0L, nrow = m+1,ncol = 1)

  cv = getcv2(signif,eps1)
  dv[1,1] = 0

  if (p == 0){
    y_rev = rot90(rot90(y))
    z_rev = rot90(rot90(z))
    vssrev = ssr(1,y_rev,z_rev,h,bigT)
    vssr = ssr(1,y,z,h,bigT)
    out = partione(h,bigT-h,bigT,vssr,vssrev)
    datx = out$dx
    ssrmin = out$ssrmin
  }
  else{
    out = onebp(y,z,x,h,1,bigT)
    datx = out$bd
    ssrmin = out$ssrind
  }

  dv[2,1] = datx

  ftest=pftest(y,z,1,q,bigT,dv[2,1,drop=FALSE],prewhit,robust,x,p,hetdat,hetvar)

  if (ftest < cv[q,1]) {
    nbreak = 0
    dv[2,1] = 0
    #dv0 = dv[seq(2,nbreak+1,1),1]
    nseg = 1
  }
  else{
    #print(paste('First break found at:',datx))
    nbreak = 1
    nseg = 2
    dv[nseg+1,1] = bigT
  }

  while(nseg <= m){
    ds = matrix(0L,nseg+1,1)
    ftestv = matrix(0L,nseg+1,1)

    i_s = 1

    while(i_s <= nseg){
      length = dv[i_s+1,1] - dv[i_s,1]

      if(length >= 2*h){
        if(p==0){
          y_temp = y[seq(dv[i_s,1]+1,dv[i_s+1,1]),1,drop=FALSE]
          z_temp = z[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          vssr = ssr(1,y_temp,z_temp,h,length)
          y_temp_rev = rot90(rot90(y_temp))
          z_temp_rev = rot90(rot90(z_temp))
          vssrev = ssr(1,y_temp_rev,z_temp_rev,h,length)
          out = partione(h,length-h,length,vssr,vssrev)
          ds[i_s,1] = out$dx
          ftestv[i_s,1] = pftest(y_temp,z_temp,1,q,length,ds[i_s,1,drop=FALSE],prewhit,
                                 robust,0,p,hetdat,hetvar)
        }
        else{
          y_temp = y[seq(dv[i_s,1]+1,dv[i_s+1,1]),1,drop=FALSE]
          z_temp = z[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          x_temp = x[seq(dv[i_s,1]+1,dv[i_s+1,1]),,drop=FALSE]
          out = onebp(y,z,x,h,dv[i_s,1]+1,dv[i_s+1,1])
          ds[i_s,1] = out$bd - dv[i_s,1]
          ftestv[i_s,1] = pftest(y_temp,z_temp,1,q,length,ds[i_s,1],
                                 prewhit,robust,x_temp,p,hetdat,hetvar)
        }
      }
      else{
        ftestv[i_s,1] = 0.0
      }
      i_s = i_s + 1
    }

    maxf = max(ftestv[seq(1,nseg,1),1])

    if (maxf < cv[q,nseg]){
      #print(nbreak)
      #dv0 = dv[seq(2,nbreak+1,1),1]
    }
    else {
      newseg = which.max(ftestv[seq(1,nseg),1])
      dv[nseg+2,1] = ds[newseg,1] + dv[newseg,1]
      nbreak = nbreak + 1
      #check this sort
      dv2 = sort(dv[seq(2,nseg+2,1),1])
      dv2 = matrix(dv2, ncol = 1)
      dv[1,1] = 0
      dv[seq(2,nseg+2,1),1] = dv2

    }
    nseg = nseg + 1
  }

  #print('The sequential procedure has reached the upper limit')
  if (nbreak < 1) {dv0 = c()}
  else{
    dv0 = dv[seq(2,nbreak+1,1),1]}
  out = list('nbreak' = nbreak, 'dv0' = dv0)
}

#' Sequential procedure
#'
#'function to apply sequential procedure to obtain number of breaks and break
#'dates. Current version only allows pure structural changes. This will be
#'generalized
#'
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param prewhit option to use AR(1) for prewhitening process
#'@param robust,hetdat,hetvar options on error terms assumptions
#' @return A list that contains following:
#' \itemize{
#'\item{nbreak}{Number of breaks}
#'\item{dateseq}{Sequence of break dates}}
#'@export
dosequa = function(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
                   robust,hetdat,hetvar) {
  nbreak = matrix(0L, nrow = 4, ncol = 1)
  dateseq = matrix(0L,nrow = 4, ncol = m)
  siglev=matrix(c(10,5,2.5,1),4,1)

  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  bigT = dim(y)[1]
  h = round(eps1*bigT)

  for (j in 1:4){
   # print(paste('Output from the sequential procedure at significance level',
   #              siglev[j,1],'%'))
    out = sequa(m,j,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
    nbr = out$nbreak
    datese = as.matrix(out$dv0)
    #print(paste('The sequential procedure estimated the number of breaks at:',nbr))
    if (nbr > 0) {
    #  print('The break dates are:')
    #  print(datese)
    }
    nbreak[j,1] =nbr

    if (nbr!=0){
      dateseq[j,seq(1,nbreak[j,1])] = t(datese)
    }
  }
  m = nbreak[2,1]
  date = as.matrix(dateseq[2,seq(1,nbreak[2,1],1),drop=FALSE])
  date = t(date)
  if (m == 0){
    print('There are no breaks in this model and estimation is skipped')
    return(NULL)
    }
  else{
    hetq=0
    hetomega=0
    out = estim(m,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = 'dosequa'
    out$nbreak = m
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out = compile.model(out)
    return(out)
  }
  
}

#preparti
preparti = function(y,z,nbreak,dateseq,h,x,p) {
  bigT = dim(z)[1]
  q = dim(z)[2]

  #care if nbreak is matrix or scalar
  dv = matrix(0L, nrow = nbreak+2, ncol = 1)
  dv[1,1] = 0
  dv[seq(2,nbreak+1,1),1] = dateseq

  dv[nbreak+2,1] = bigT
  ds = matrix(0L,nrow = nbreak, ncol = 1)
  dr = matrix(0L,nrow = nbreak, ncol = 1)

  for (is in 1:nbreak){
    length = dv[is+2,1] - dv[is,1]
    if (p == 0){
      index = seq(dv[is,1]+1,dv[is+2,1],1)
      y_temp = y[index,1,drop=FALSE]
      z_temp = z[index,,drop=FALSE]
      vssr = ssr(1,y_temp,z_temp,h,length)
      y_temp_rev = rot90(rot90(y_temp))
      z_temp_rev = rot90(rot90(z_temp))
      vssrev = ssr(1,y_temp_rev,z_temp_rev,h,length)
      out = partione(h,length-h,length,vssr,vssrev)
      ds[is,1] = out$dx
      dr[is,1] = ds[is,1] + dv[is,1]
    }
    else{
      out = onebp(y,z,x,h,dv[is,1]+1,dv[is+2,1])
      ds[is,1] = out$bd
      dr[is,1] = ds[is,1]
    }

  }
  return(dr)
}

#'Repartition procedure
#'
#'The following procedure constructs the so-called repartition
#'estimates of the breaks obtained by the sequential method (see Bai
#'(1995), Estimating Breaks one at a time, Econometric Theory, 13,
#'315-352. It alows estimates that have the same asymptotic
#'distribution as those obtained by global minimization. Otherwise, the
#'output from the procedure "estim" below do not deliver asymptotically
#'correct confidence intervals for the break dates.
#'
#'@param y dependent variables in matrix form
#'@param z matrix of independent variables with coefficients are allowed to change across
#'regimes
#'@param x matrix of independent variables with coefficients constant across regimes
#'@param m maximum number of breaks
#'@param eps1 trimming level
#'@param eps convergence criterion for iterative recursive computation
#'@param maxi maximum number of iterations
#'@param fixb option to use fixed initial input \eqn{\beta}. If \code{1},
#'the model will use values given in \code{betaini}. If \code{0}, betaini is skipped
#'@param betaini Initial \eqn{beta_0} to use in estimation
#'@param printd option to print results of iterations for partial change model
#'@param prewhit option to use AR(1) for prewhitening process
#'@param robust,hetdat,hetvar options on error terms assumptions
#'@return reparv Repartition method estimation of break dates
#'
#'@export
dorepart = function(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
                    robust,hetdat,hetvar){
  reparv = matrix (0L,4,m)
  siglev=matrix(c(10,5,2.5,1),4,1)

  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  bigT = dim(y)[1]
  h = round(eps1*bigT)
  nbreak = matrix(0L,4,1)
  
  
  for (j in 1:4){
    temp = sequa(m,j,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
    nbreak[j,1] = temp$nbreak
    if (temp$nbreak == 0){
     # print(('The sequential procedure found no break and the repartition procedure is skipped.'))
    }
    else {
      repartda = preparti(y,z,temp$nbreak,
                          temp$dv0,
                          h,x,p)
      reparv[j,seq(1:temp$nbreak)] = repartda
    }
  }
  
  #estimate the date at 5% significant level
  m = nbreak[2,1]
  date = as.matrix(reparv[2,seq(1,nbreak[2,1],1),drop=FALSE])
  date = t(date)
  
  if (m == 0){
    print('There are no breaks in this model and estimation is skipped')
    return(NULL)
  }
  else{
    hetq=0
    hetomega=0
    out = estim(m,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
    out$p_name = 'dorepart'
    out$nbreak = m
    class(out) = 'model'
    out$numz = q
    out$numx = p
    out = compile.model(out)
    return(out)
  }
  
  
  return (reparv)
}


#Format output of n break model
#'@export
compile.model = function (x,digits = -1,...){
  if (digits==-1){digits = 3}
  if (x$nbreak == 0){return(NULL)}
  else {
    
    #format date estimation
    CI_95 = c()
    CI_90 = c()
    coln = c()
    for (i in 1:x$nbreak){
      coln = c(coln, paste('Break',i,sep=''))
      CI_95 = c(CI_95,paste('(',x$CI[i,1],',',x$CI[i,2],')',sep=''))
      CI_90 = c(CI_90,paste('(',x$CI[i,3],',',x$CI[i,4],')',sep=''))
    }
    date_tab = data.frame(date = t(x$date),stringsAsFactors = FALSE)
    date_tab = rbind(date_tab,CI_95)
    date_tab = rbind(date_tab,CI_90)
    colnames(date_tab) = coln
    rownames(date_tab) = c('Date','95% CI','90% CI')
    
    
    #format regime-specific coefficients
    rnameRS = c()
      for (i in 1:x$numz){
        rnameRS = cbind(rnameRS,paste('z',i,sep=''))
      }
    
    cnameRS = c()
    coefRS = c()
    for (i in 1:(x$nbreak+1)){
    cnameRS = cbind(cnameRS,paste('Regime',i))
    }
    for (j in 1:x$numz){
      coefRSj = c()
    for (i in 1:(x$nbreak+1)){
      coefRSj = cbind(coefRSj,paste(format(round(x$beta[(j-1)*(x$nbreak+1)+i,1],digits),nsmall=digits),
                              paste('(',format(round(x$SE[(j-1)*(x$nbreak+1)+i,1],digits),nsmall=digits),')',sep='')))
      }
      coefRS = rbind(coefRS,coefRSj)
    }
    
  }
    rnameRS = paste(rnameRS,'(SE)',sep=' ')
    coef_tabRS=data.frame(co = coefRS)
    
    
    rownames(coef_tabRS) = rnameRS
    colnames(coef_tabRS) = cnameRS
    
    
    #format regime-wise coefficients
    
    if(x$numx == 0){coef_tabRW = NULL}
    else{ 
      rnameRW = c()
      for (i in 1:x$numx){
        rnameRW = cbind(rnameRW,paste('x',i,sep=''))
      }
    
    cnameRW = 'All regimes'
    coefRW = c()
    for (j in 1:x$numx){
      coefRW = cbind(coefRW,paste(format(round(x$beta[x$numz*(x$nbreak+1)+j,1],digits),nsmall=digits),
                                  paste('(',format(round(x$SE[x$numz*x$nbreak+j,1],digits),nsmall=digits),')',sep='')))}
    
    rnameRW = paste(rnameRW,'(SE)',sep=' ')
    coef_tabRW=data.frame(co = coefRW)
    
    rownames(coef_tabRW) = rnameRW
    colnames(coef_tabRW) = cnameRW}
    
    table = list('date_tab' = date_tab,'RS_tab' = coef_tabRS, 'RW_tab' = coef_tabRW)
    x$tab = table
    return(x)
  }
  
  




#'Summary output of a n breaks model
#'
#'Function to format the output of the n-break model
#'@param nbreak number of breaks in the model
#'@param glb global minimum SSR of the model
#'@param datevec estimated break date
#'@param y dependent variables
#'@param z independent variables with regime-specific coefficients
#'@param x independent variables with regime-wise coefficients
#'@return tbl Tables containings:
#'i) global minimum SSR
#'ii) estimated date
#'iii) estimated coefficients

print.model <- function(x,digits = -1,...)
{
  #print procedure to select number of breaks
  proc = switch(x$p_name,'dosequa'='sequential procedure', 'BIC' = 'BIC', 'LWZ' = 'LWZ',
                'dorepart' = 'repartition procedure', 'fix' = 'specified number of breaks')
  cat(paste('\nThe number of breaks is estimated by',proc,'\n'))
  
  if (digits==-1){digits = 3}
  if (x$nbreak == 0){print('No break found due to estimated number of break is 0')
    break}
  if(x$numx == 0){
    cat(paste('Pure change model with',x$nbreak,'estimated breaks.',sep=' '))
  }else if(x$numx > 0){
    cat(paste('Partial change model with', x$nbreak,'estimated breaks.',sep=' '))
  }
  cat('\nMinimum SSR =',
              format(round(x$SSR,3),nsmall=3),'\n')
  
  cat('\nEstimated date:\n')
  print(x$tab$date_tab,quote=FALSE)
  
  cat('\nEstimated regime-specific coefficients:\n')
  print(x$tab$RS_tab,quote=FALSE)
  
  
  if(x$numx == 0) {cat('\nNo regime-wise regressors\n')}
  else{
    cat('\nEstimated regime-wise coefficients:\n\n')
    print(x$tab$RW_tab,quote='FALSE')}
    

  invisible(x)
}


#'Summary output of Sup Wald test
#'
#'Function to format the output of the Sup F test
#'@param nbreak number of breaks in the model
#'@param glb global minimum SSR of the model
#'@param datevec estimated break date
#'@param y dependent variables
#'@param z independent variables with regime-specific coefficients
#'@param x independent variables with regime-wise coefficients
#'@return tbl Tables containings:
#'i) global minimum SSR
#'ii) estimated date
#'iii) estimated coefficients
#'


compile.Wtest <- function(x,digits = -1,...)
{
  cnames1 = c()
  for (i in 1:x$mbreak){
    if (i == 1){
    cnames1 = cbind(cnames1, paste(i,'break'))}
    else
    {
      cnames1 = cbind(cnames1,paste(i,'breaks'))
    }
  }
  ftest = t(x$ftest)
  
  
  supF1 = data.frame(ftest = format(round(ftest,3),nsmall=3)) 
  cv_supF = format(round(x$cv_supF,3),nsmall = 3)
  colnames(cv_supF) = colnames(supF1)
  supF1 = rbind(supF1,cv_supF)
 
  rownames(supF1) = c('Sup F','10% CV','5% CV','2.5% CV','1% CV')
  colnames(supF1) = cnames1
  
  UDmax = data.frame(UDmax = format(round(x$UDmax,3),nsmall = 3), 
                     cv = format(round(t(x$cv_Dmax),3),nsmall = 3))

  colnames(UDmax) = c('UDMax','10% CV','5% CV','2.5% CV','1% CV')
 
  
  x$supF1 = supF1
  x$UDmax = UDmax
  return(x)
}

#'Summary output of Sup(l+1|l) test
#'
#'Function to format the output of the Sup F test
#'@param nbreak number of breaks in the model
#'@param glb global minimum SSR of the model
#'@param datevec estimated break date
#'@param y dependent variables
#'@param z independent variables with regime-specific coefficients
#'@param x independent variables with regime-wise coefficients
#'@return tbl Tables containings:
#'i) global minimum SSR
#'ii) estimated date
#'iii) estimated coefficients
#'


print.Wtest <- function(x)
{
  cat('\na) SupF tests against a fixed number of breaks\n\n')
  print(x$supF1,quote=FALSE)
  cat('\nb) UDmax tests against an unknown number of breaks\n\n')
  print(x$UDmax,quote=FALSE)
  invisible(x)
}

compile.supfltest = function(x){
  nbreak = sum(!is.na(x$ndat))
  cnames = c()
  for (i in 1:nbreak){
    cnames = cbind(cnames, paste('supF(',i+1,'|',i,')',sep=''))
  }
  x$supfl = format(round(x$supfl,3),nsmall = 3)
  x$ndat = format(x$ndat,nsmall=0)
  x$cv = format(round(x$cv,3),nsmall = 3)
  
  sfl =data.frame(supfl = t(x$supfl[seq(1,nbreak,1),1,drop=FALSE]))
  ndat = t(x$ndat[seq(1,nbreak,1),1,drop=FALSE])
  cv = x$cv[,seq(1,nbreak,1)]
  colnames(ndat) = colnames(sfl)
  colnames(cv) = colnames(sfl)
  sfl = rbind(sfl,ndat)
  sfl = rbind(sfl,cv)
  colnames(sfl) = cnames
  rownames(sfl) = c('Seq supF','Estimated date','10% CV','5% CV', '2.5% CV', '1% CV')
  x$sfl = sfl
  return (x)
}

print.supfltest = function(x){
  cat('\nsupF(l+1|l) tests using global optimizers under the null\n\n')
  print(x$sfl,quote=FALSE)
  invisible(x)
}


