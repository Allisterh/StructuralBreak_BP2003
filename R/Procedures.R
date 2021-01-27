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
#'\item{glb}{minimum global SSR}
#'\item{datevec}{Vector of dates (optimal minimizers)}
#'\item{bigvec}{Associated SSRs}
#'@export
#'
doglob = function (y,z,x,m,eps,eps1,maxi,fixb,betaini,printd){

#check if model is pure or partial change
if (is.null(x)) {p = 0}
else {p = dim(x)[2]}

q = dim(z)[2]
h = round(eps1*bigT)
bigT = dim(y)[1]

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

#printing results
if (printd == 1) {
for (i in 1:m){
  print(paste('Model with',i,'breaks has SSR:',glb[i,1]))
  print('The dates of breaks are:')
  print(datevec[1:i,i])
}
}
  return (out)
}

#' SupF, UDMax & WDMax testing procedure
#'
#' The procedure calculate the test statistics and print results of the 2 main tests:
#' \item{SupF test} {F test of 0 vs m breaks}
#' \item{Double Max test} {UDMax: the unweighted version
#' and WDMax: the weighted version}
#'
#'
#'@return A list that contains following:
#'\item{ftest}{SupF test statistics}
#'\item{cv}{Critical values for Sup F test }
#'\item{wftest}{Weighted Double Max test}
#'\item{cvm}{Critical values for Dmax test}
#'@export
#'
#'
dotest = function(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,robust,
                  hetdat,hetvar){
  siglev=matrix(c(10,5,2.5,1),4,1)
  printd = 0 #suppress output printing
  print('Output from testing procedure')
  out = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
  datevec = out$datevec

  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  h = round(eps1*bigT)
  bigT = dim(y)[1]

  #procedure for F test
  print('a) supF tests against a fixed number of breaks')

  ftest = matrix(0L, nrow = m, ncol = 1)
  wftest = matrix(0L, nrow = m, ncol = 1)

  for (i in 1:m){
    ftest[i,1] = pftest(y,z,i,q,bigT,datevec,prewhit,robust,x,p,hetdat,hetvar)
    print(paste('supF test for 0 versus',i,'breaks (scaled by q):',ftest[i,1]))
  }

  cv_supF = matrix(0L,4,m)
  for (c in 1:4){
    #critical values for supF test
    cv = getcv1(c,eps1)
    cv_supF[c,] = cv[q,1:m,drop=FALSE]
    print(paste('The critical values at the',siglev[c,1],'% level are (for k = 1 to',m,'):'))
    print(cv[q,1:m,drop=FALSE])
  }

  #procedure for Dmax and UDmax test

  print('b) Dmax test against an unknown number of breaks')
  print(paste('The UDmax test is:',max(ftest)))
  cv_Dmax = matrix(0L,4,1)
  for (c in 1:4) {
    #critical values for Dmax test
    cvm = getdmax(c,eps1)
    cv_Dmax[c,1] = cvm[q,1]
    print(paste('The critical values at the',siglev[c,1],'% level is:',
                cvm[q,1]))
  }


  for (c in 1:4){
    #computation of WDmax test
    cv = getcv1(c,eps1)
    for( i in 1:m){
      wftest[i,1] = cv[q,1] * ftest[i,1] / cv[q,1]
    }
    print(paste('WDmax test at the',siglev[c,1],'% level is:',max(wftest)))
  }


  rownames(cv_supF) = siglev
  rownames(cv_Dmax) = siglev
  out = list('ftest' = ftest, 'wftest' = wftest, 'cv_supF' = cv_supF,
             'cv_Dmax' = cv_Dmax)
  return(out)
}



#'
#'

dospflp1 = function(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
                    robust,hetdat,hetvar) {
  siglev=matrix(c(10,5,2.5,1),4,1)
  print('supF(l+1|l) tests using global optimizers under the null')
  out = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
  datevec = out$datevec
  bigvec = out$bigvec

  if (is.null(x)) {p = 0}
  else {p = dim(x)[2]}

  q = dim(z)[2]
  h = round(eps1*bigT)
  bigT = dim(y)[1]

  supfl = matrix (0L,m-1,1)
  ndat = matrix (0L,m-1,1)
  for (i in seq(1,m-1,1)){
    out1 = spflp1(bigvec,datevec[1:i,i,drop=FALSE],i+1,y,z,h,q,prewhit,robust,x,p,hetdat,hetvar)
    supfl[i,1] = out1$maxf
    ndat[i,1] = out1$newd
    print(paste('The supF(',i+1,'|',i,') test is',supfl[i,1]))
    print(paste('It corresponds to a new break at:',ndat[i,1]))
  }

  cv_supFl = matrix(0L,4,m)

  for (c in 1:4){
    #critical values for supF(l+1|l) test
    cv = getcv2(c,eps1)
    cv_supFl[c,] = cv[q,1:m,drop=FALSE]
    print(paste('The critical values at the',siglev[c,1],'% level are (for k = 1 to',m,'):'))
    print(cv[q,1:m,drop=FALSE])
  }
  rownames(cv_supFl) = siglev


  result = list('supfl' = supfl, 'ndat' = ndat, 'cv' = cv_supFl)
}

#'
#'
#'
doorder = function(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd) {
  if (p == 0){zz = z}
  else{zz = cbind(z,x)}
  out = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
  glb = out$glb
  bigT = dim(y)[1]
  ssr0 = nssr(y,zz)
  out = order(ssr0,glb,bigT,m,q)
  #store result
  mbic = out$mBIC
  mlwz = out$mLWZ
}


#
dosequa = function() {
  nbreak = matrix(0L, nrow = 4, ncol = 1)
  dateseq = matrix(0L,nrow = 4, ncol = m)

  for (j in 1:4){
    print(paste('Output from the sequential procedure at significance level',
                siglev[j,1],'%'))

    out = sequa(m,j,q,h,bigT,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
    nbr = out$nbreak
    datese = out$dv0
    print(paste('The sequential procedure estimated the number of breaks at:',nbr))
    if (nbr > 0) {
      print('The break dates are:')
      print(datese)
    }
    nbreak[j,1] =nbr

    if (nbr!=0){
      dateseq[j,seq(1,nbreak[j,1])] = t(datese)
    }

  }
}

#'
dorepart = function(){}

#
