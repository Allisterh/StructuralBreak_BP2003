#procedure and sub-procedures to obtain corrected cov-variance estimators
#based on assumptions on the error terms (e.g: heteroskedasticity/homoskedasticity,
# serially correlated/serially uncorrelated, ...)


#' Long-run covariance matrix computation
#'
#' The procedure to compute the long run covariance matrix \code{jhat} based
#' on variance matrix \code{vmat}
#'
#' @param vmat variance matrix
#' @return jhat Long run covariance matrix
#' @export
jhatpr = function(vmat) {
  nt = dim(vmat)[1]
  d = dim(vmat)[2]
  #automatic bandwidth selection
  st = bandw(vmat)
  #lag 0 covariance
  jhat = t(vmat) %*% vmat

  #forward sum
  for (j in seq(1,nt-1,1)){
    vl = vmat[(j+1):nt,,drop=FALSE]
    vr = vmat[1:(nt-j),,drop=FALSE]
    vmat_s = t(vl) %*% vr
    jhat = jhat + as.vector(kern(j/st)) * vmat_s
  }

  #backward sum
  for (j in seq(1,nt-1,1)) {
    vl = vmat[seq(1,nt-j,1),,drop=FALSE]
    vr = vmat[seq(j+1,nt,1),,drop=FALSE]
    vmat_s = t(vl) %*% vr
    jhat = jhat + as.vector(kern(j/st)) * vmat_s
  }

  #small sample correction
  jhat = jhat/(nt-d)

  return(jhat)
}

#' Heteroskedasticy and autocorrelation consistency correction
#'
#'function corrects the standard errors based on options of prewhitening
#'using a AR(1) process estimation of error terms to obtain HAC errors/
#'
#'@param reg matrix of regressors
#'@param res residuals
#'@param prewhit Option of using prewhitening process. If \code{1}, an AR(1)
#'process will be used to filter. If \code{0}, skipped the filtering process
#'@return hac Heteroskedasticy and autocorrelation consistent errors
####
correct = function (reg,res,prewhit){
  #initialize storage
  nt = dim(reg)[1]
  d = dim(reg)[2]
  b = matrix(0L, nrow = d, ncol = 1)
  bmat = matrix(0L, nrow = d, ncol = d)
  vstar = matrix(0L, nrow = nt-1, ncol = d)
  vmat = matrix(0L, nrow = nt, ncol = d)

  #Construct matrix z_t * u_t
  for (i in 1:d){
    vmat[,i] = reg[,i,drop=FALSE] * res
  }

  #Prewhitening to matrix vmat by filtering with a AR(1). If prewhit = 0, skip
  if (prewhit == 1){
    for(i in 1:d){
      #check carefully if errors
      b = OLS(vmat[seq(2,nt,1),i,drop=FALSE],vmat[seq(1,nt-1,1),,drop=FALSE])
      bmat[i,] = t(b)
      vstar[,i] = vmat[seq(2,nt,1),i,drop=FALSE] - vmat[seq(1,nt-1,1),,drop=FALSE] %*% b
    }
    #kernel on the residuals
    jh = jhatpr(vstar)

    #recolor
    inv = solve( diag(1,d) - bmat)
    hac = inv %*% jh %*% t(inv)
  }
  else {

    hac = jhatpr(vmat)
  }

  return (hac)
}
#' Covariance matrix of estimator delta construction
#'
#'@param y dependent variable
#'@param z matrix of independent variables with coefficients allowed to change
#'across regimes
#'@param x atrix of independent variables with constant coefficients
#'across regimes
#'@param q number of regressors \code{z}
#'@param p number x regressors
#'@param i maximum number of breaks
#'@param bigT sample period T
#'@param b vector of estimated dates of breaks
#'@param prewhit Option of using prewhitening process. If \code{1}, an AR(1)
#'process will be used to filter. If \code{0}, skipped the filtering process
#'@param robust,withb,hetdat,hetvar options for assumptions of error terms
#'structure
#'@return vdel Covariance matrix of delta
#'@export

pvdel = function(y,z,i,q,bigT,b,prewhit,robust,x,p,withb,hetdat,hetvar) {

  ev = matrix(0L , nrow = i+1, ncol = 1)
  zbar = diag_par(z,i,b)

  if (p == 0){
    delv = OLS(y,zbar)
    res = y - zbar %*% delv
    reg = zbar
  }
  else{
    delv = OLS(y, cbind(zbar,x))
    res = y - cbind(zbar,x) %*% delv

    if (withb == 0) {reg = zbar - x %*% solve(t(x) %*% x) %*% t(x) %*% zbar}
    else {reg = cbind(x,zbar)}
  }


  if (robust == 0) {
    #testing with no serial correlation in errors
    if(p==0) {
      if (hetdat==1 && hetvar == 0){
        sig = t(res) %*% res / bigT
        vdel = drop(sig) * solve(t(reg) %*% reg)

      }
      if (hetdat == 1 && hetvar == 1){
        sig = psigmq(res,b,q,i,bigT)
        vdel = kron(sig,diag(1,q)) %*% solve(t(reg) %*% reg)
      }
      if (hetdat == 0 && hetvar == 0){
        lambda = plambda(b,i,bigT)
        sig = t(res) %*% res / bigT
        vdel = drop(sig) * solve(kron(lambda,t(z) %*% z))
      }
      if (hetdat == 0 && hetvar == 1) {
        lambda = plambda(b,i,bigT)
        sig = psigmq(res,b,q,i,bigT)
        vdel = kron(sig,diag(1,q)) %*% solve(kron(lambda, t(z) %*% z))
      }
    }
    else {
      if (hetdat == 0) {
        print(paste('hetdat == 0 is not allowed','vdel is returned zeros'))
        vdel = matrix (0L, nrow = q*(i+1), ncol = q*(i+1))
      }
      if (hetdat == 1 && hetvar == 0) {
        sig = t(res) %*% res / bigT
        vdel = drop(sig) * solve(t(reg) %*% reg)
      }
      if (hetdat == 1 && hetvar == 1) {
        wbar = diag_par(reg,i,b)
        ww = t(wbar) %*% wbar
        sig = psigmq(res,b,q,i,bigT)
        gg = matrix (0L, nrow = (i+1)*q + p*withb, ncol = (i+1)*q + p*withb)
        ie = 1
        while(ie <= i + 1){
          index = seq((ie-1)*((i+1)*q+p*withb)+1,ie*((i+1)*q+p*withb),1)
          increment = sig[ie,ie] * ww[index,index]
          gg = gg + increment
          ie = ie + 1
        }
        vdel = solve(t(reg) %*% reg) %*% gg %*% solve(t(reg) %*% reg)
      }
    }
  }
  else {
    #testing with serial correlation in errors
    if(hetdat == 0) {
      print(paste('hetdat = 0 is not allowed','vdel is returned zeros'))
      vdel = matrix(0L,nrow = q*(i+1),ncol = q*(i+1))
    }

    if (p==0){
      if (hetvar == 1){
        hac = matrix(0L, nrow = q*(i+1), ncol = q*(i+1))
        vdel = matrix(0L, nrow = q*(i+1), ncol = q*(i+1))
        ind_temp = seq(1,b[1,1],1)
        temp = correct(z[ind_temp,,drop=FALSE], res[ind_temp,1,drop=FALSE] , prewhit)
        hac[1:q,1:q] = b[1,1] * temp
        if(dim(b)[1] > 1){
          for (j in 2:i) {
            ind_hac = seq((j-1)*q+1,j*q,1)
            ind_temp = seq(b[j-1,1]+1,b[j,1],1)
            temp = correct(z[ind_temp,,drop=FALSE], res[ind_temp,1,drop=FALSE], prewhit)
            hac[ind_hac,ind_hac] = (b[j,1] - b[j-1,1]) * temp
          }
        }
        ind_hac = seq(i*q+1,(i+1)*q,1)
        ind_temp = seq(b[i,1]+1,bigT,1)
        temp = correct(z[ind_temp,,drop=FALSE],res[ind_temp,1,drop=FALSE],prewhit)
        hac[ind_hac,ind_hac] = (bigT - b[i,1]) * temp
        vdel = solve(t(reg) %*% reg) %*% hac %*% solve(t(reg) %*% reg)

      }
      else {
        hac = correct(z,res,prewhit)
        lambda = plambda(b,i,bigT)
        vdel = bigT * solve(t(reg) %*% reg) %*% kron(lambda,hac) %*% solve(t(reg) %*% reg)
      }
    }
    else{
      hac = correct(reg,res,prewhit)
      vdel = bigT * solve(t(reg) %*% reg) %*% hac %*% solve(t(reg) %*% reg)
    }
  }

  return(vdel)
}

