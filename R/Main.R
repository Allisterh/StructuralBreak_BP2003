#Main program to handle users' input

#' Robust structural change estimation
#'
#' Function executes main procedures described in Bai and Perron
#'
#' Function executes the main procedures selected by users. There are 7 main
#' procedures including:\itemize{
#' \item{\code{\link{doglob}}}{ A recursive procedure to obtain global minimizers of SSR,
#' corresponding to the break date.}
#' \item{\code{\link{dotest}}} {Procedure to conduct SupF test of 0 versus m breaks and
#'  Double Max test.}
#' \item{\code{\link{dosupflp1}}} {Procedure to conduct SupF(l+1|l).}
#' \item{\code{\link{doorder}}} {Procedure to find number of break by criteria selection}
#' \item{\code{\link{dosequa}}} {Procedure to obtain break dates via sequential method}
#' \item{\code{\link{dorepart}}} {Procedure to obtain break dates via repartition method}
#'}
#' To obtain the confidence intervals for the break dates and corrected errors for
#' estimates of the model, set \code{est} = \code{1}. If \code{est} = \code{1}, specify
#' one of the following methods for estimation procedure:
#' \itemize{
#' \item{\code{"BIC"}} {Estimate and compute confidence intervals following
#' order selection method to pick number of breaks by BIC. For more details, see
#' \code{\link{doorder}},\code{\link{estim}}}
#' \item{\code{"LWZ"}} {Estimate and compute confidence intervals following
#' order selection method to pick number of breaks by LWZ. For more details, see
#' \code{\link{doorder}},\code{\link{estim}}}
#' \item{\code{"seq"}} {Estimate and compute confidence intervals following
#' sequential method estimation of the break date. For more details, see
#' \code{\link{dosequa}},\code{\link{estim}}}
#' \item{\code{"rep"}} {Estimate and compute confidence intervals following
#' repartition method estimation of the break date. For more details, see
#' \code{\link{dorepart}},\code{\link{estim}}}
#' \item{\code{"fix"}} {Estimate and compute confidence intervals following
#' pre-specified number of breaks, \code{i}. For more details, see
#' \code{\link{estim}}}}
#' \emph{All \code{default} values} of error assumptions (\code{robust},
#' \code{hetdat}, \code{hetvar}, \code{hetq}) will be set to 1
#'
#' @param y name of dependent variable
#' @param z name of independent variables which coefficients are allowed to change
#' across regimes. \code{default} is vector of 1 (Mean-shift model)
#' @param x name of independent variables which coefficients are constant across
#' regimes. \code{default} is NULL
#' @param data the data set for estimation
#' @param eps1 value of trimming (in percentage) for the construction
#' and critical values of \itemize{
#' \item {supF type tests} {supF test, the Dmax (\code{\link{pftest}})
#' and the supF(l+1|l) (\code{\link{spflp1}})
#' \item  {sequential procedure}. }
#' If these tests are used, minimal segment length will be set
#' at \code{default} = int(\code{eps1}*T) (T is total sample size).
#' But if the tests are not required,
#' estimation can be done with an arbitrary h. There are five
#' options: \itemize{
#' \item{\code{eps1}=.05} {Maximal value of \code{m} = 10}
#' \item{\code{eps1}=.10} {Maximal value of \code{m} = 8}
#' \item{\code{eps1}=.15} {Maximal value of \code{m} = 5}
#' \item{\code{eps1}=.20} {Maximal value of \code{m} = 3}
#' \item{\code{eps1}=.25} {Maximal value of \code{m} = 2}
#' } \code{default} = 0.15
#' @param m Maximum number of structural changes allowed. If not specify,
#' m will be set to \code{default} value from \code{eps1}
#' @param prewhit set to \code{1} if want to apply AR(1) prewhitening prior to estimating
#' the long run covariance matrix.
#' @param robust set to \code{1} if want to allow for heterogeneity
#' and autocorrelation in the residuals, \code{0} otherwise.
#' The method used is \emph{Andrews(1991)} automatic bandwidth with AR(1) approximation and the quadratic
#' kernel. Note: Do not set to \code{1} if lagged dependent variables are
#' included as regressors.
#' @param hetdat option for the construction of the F tests. Set to 1 if want to
#' allow different moment matrices of the regressors across segments.
#' If \code{hetdat} = \code{0}, the same moment matrices are assumed for each segment
#' and estimated from the ful sample. It is recommended to set
#' \code{hetdat}=\code{1} if number of regressors \code{x} > \code{0}.
#' @param hetvar option for the construction of the F tests.Set to \code{1}
#' if want to allow
#' for the variance of the residuals to be different across segments.
#' If \code{hetvar}=\code{0}, the variance of the residuals is assumed constant
#' across segments and constructed from the full sample. This option
#' is not available when \code{robust}=\code{1}
#' @param hetomega used in the construction of the confidence intervals for the break
#' dates. If \code{hetomega}=\code{0}, the long run covariance matrix of zu is
#' assumed identical across segments
#' (the variance of the errors u if \code{robust}={0})
#' @param hetq used in the construction of the confidence intervals for the break
#' dates. If \code{hetq}=\code{0}, the moment matrix of the data is assumed identical
#' across segments
#' @param procedure Selectable procedures to execute:
#' \itemize{
#' \item{\code{\link{doglob}}}{ A recursive procedure to obtain global minimizers of SSR,
#' corresponding to the break date.}
#' \item{\code{\link{dotest}}} {Procedure to conduct SupF test of 0 versus m breaks and
#'  Double Max test.}
#' \item{\code{\link{dospflp1}}} {Procedure to conduct SupF(l+1|l).}
#' \item{\code{\link{doorder}}} {Procedure to find number of break by criteria selection}
#' \item{\code{\link{dosequa}}} {Procedure to obtain break dates via sequential method}
#' \item{\code{\link{dorepart}}} {Procedure to obtain break dates via repartition method}
#'}
#' @param CI Option to estimate dates, model parameters and their confidence intervals.
#' If \code{1}, specify method to estimate in \code{method} below.
#' @param method \itemize{
#' \item{\code{"BIC"}} {Estimate and compute confidence intervals following
#' order selection method to pick number of breaks by BIC. For more details, see
#' \code{\link{doorder}},\code{\link{estim}}}
#' \item{\code{"LWZ"}} {Estimate and compute confidence intervals following
#' order selection method to pick number of breaks by LWZ. For more details, see
#' \code{\link{doorder}},\code{\link{estim}}}
#' \item{\code{"seq"}} {Estimate and compute confidence intervals following
#' sequential method estimation of the break date. For more details, see
#' \code{\link{dosequa}},\code{\link{estim}}}
#' \item{\code{"rep"}} {Estimate and compute confidence intervals following
#' repartition method estimation of the break date. For more details, see
#' \code{\link{dorepart}},\code{\link{estim}}}
#' \item{\code{"fix"}} {Estimate and compute confidence intervals following
#' pre-specified number of breaks, \code{i}. For more details, see
#' \code{\link{estim}}}}
#' @param maxi number of maximum iterations for recursive calculations of finding
#' global minimizers.\code{default} = 10
#' @param fixb Option to use initial values of beta. Set to 1 if initial values
#' are entered by users.
#' @param eps convergence criterion for recursive calculations
#' @param betaini values of initial betas. Only used when \code{fixb} = 1.
#' @param fixn number of pre-specified breaks. \code{default} = -1. It will be replaced
#' automatically to 2 if no specification is given
#' @param printd Print option for model estimation. \code{default} = 0, to suppress outputs
#' @export

model_estimation = function(y_name,z_name = NULL,x_name = NULL,data,eps1 = 0.15,m = -1,prewhit = 1,
                           robust = 1,hetdat = 1,hetvar = 1,hetomega = 1,hetq = 1,
    procedure = c('doglob','dotest','dospflp1','doorder','dosequa','dorepart'),
    CI = 1,method = c('BIC','LWZ','seq','rep','fix'),maxi = 10,fixb = 0,
    eps = 0.00001,betaini = 0,fixn=-1,printd = 0){

  #handle data
  y_ind = match(y_name,colnames(data))
  x_ind = match(x_name,colnames(data))
  z_ind = match(z_name,colnames(data))

  if(is.na(y_ind)){
    print('No dependent variable found. Please try again')
    return(NULL)}

  y = data[,y_ind]
  y = data.matrix(y)
  T = dim(y)[1]

  if (is.null(x_name)) {x = c()}
  else{
    if(anyNA(x_ind)){print('No x regressors found. Please try again')}
    else{x = data.matrix(data[,x_ind,drop=FALSE])}}
  if (is.null(z_name)) {z = matrix(1L,T,1)}
  else{
    if(anyNA(z_ind)){print('No z regressors found. Please try again')}
    else{z = data.matrix(data[,z_ind,drop=FALSE])}}

  #set maximum breaks
  v_eps1 = c(0.05,0.10,0.15,0.20,0.25)
  v_m = c(10,8,5,3,2)
  index = match(eps1,v_eps1)

  #set significant level
  siglev=matrix(c(10,5,2.5,1),4,1)

  if(is.na(index)) {
    print('Invalid trimming level, please choose 1 of the following values')
    print(v_eps1)
    return(NULL)
  }
  if (m == -1) {
    m = v_m[index]
  }
  if (fixn == -1){
    fixn = 2
  }

  if(length(dim(x))==0){p = 0}
  else{p = dim(x)[2]}
  q = dim(z)[2]


  procedure = unlist(procedure,',')


  out = list('procedure' =  procedure)

  if('doglob' %in% procedure){
    t = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
    out$glob = t[c(1:2)]
  }
  if('dotest' %in% procedure){
   out$test = dotest(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,robust,
           hetdat,hetvar)

  }
  if('dospflp1' %in% procedure){
    out$flp = dospflp1(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
             robust,hetdat,hetvar)

  }
  if('doorder' %in% procedure){
    out$ord = doorder(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)

  }
  if('dosequa' %in% procedure){
   out$sequ = dosequa(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
            robust,hetdat,hetvar)

  }
  if('dorepart' %in% procedure){
    out$repart = dorepart(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
             robust,hetdat,hetvar)
  }

  if(CI == 1){
    #method = unlist(method,',')
    if('BIC' %in% method){
      t = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
      t_out = doorder(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
      mbic = t_out$BIC
      datevec = t$datevec
      if (mbic == 0) {
        print('No break selected by BIC')
        }
      else{
      out$est_BIC = estim(mbic,q,z,y,datevec[,mbic,drop=FALSE],robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
      }
    }
    if('LWZ' %in% method){
      t = doglob(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
      t_out = doorder(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd)
      mlwz = t_out$LWZ
      datevec = t$datevec
      if (mlwz == 0) {
        print('No break selected by BIC')
        }
      else{
      out$est_LWZ = estim(mlwz,q,z,y,datevec[,mlwz,drop=FALSE],robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)}
    }
    if('seq' %in% method){
      t_out = dosequa(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
                      robust,hetdat,hetvar)
      nbreak = t_out$nbreak
      dateseq = t_out$dateseq
      ii = 0
      j = 1
      while(j<=4){
        if (ii == 0){
          if (nbreak[j,1] != 0){
            print(paste('Output from the estimation of the model selected by the sequential method at significance level',
                        siglev[j,1], '%'))
            out$est_seq = estim(nbreak[j,1],q,z,y,t(dateseq[j,seq(1,nbreak[j,1]),drop=FALSE]),robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
          }
        }
        j = j+1

        if(j <= 4){
          if(nbreak[j,1] == nbreak[j-1,1]){
            if (nbreak[j,1] == 0){
              print(paste('For the',siglev[j,1], '% level, the model is the same as for the',
                          siglev[j-1,1], '% level.'))
              print('The estimation is not repeated')
              ii = 1
            }
            else{
              if (identical(dateseq[j,seq(1,nbreak[j,1]),drop = FALSE],dateseq[j-1,seq(1,nbreak[j-1,1]),drop = FALSE])){
                print(paste('For the',siglev[j,1], '% level, the model is the same as for the',
                            siglev[j-1,1], '% level.'))
                print('The estimation is not repeated')
                ii = 1
              }
            }
          }
        }
        else{ii = 0}
      }
    }

    if('rep' %in% method){
      t_out = dosequa(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
                      robust,hetdat,hetvar)
      t_out1 = dorepart(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
                        robust,hetdat,hetvar)
      nbreak = t_out$nbreak
      reparv = t_out1
      ii = 0
      while (j <= 4){
        if (ii==0){
          if (nbreak[j,1] == 0){
            print(paste('The sequential procedure at the significance level',
                        siglev[j,1],'% found no break and'))
            print('the repartition procedure was skipped')
          }
          else{
            print(paste('Output from the estimation of the model selected by the repartition method from the sequential procedure at the significance level',
                        siglev[j,1],'%'))
            out$est_repart = estim(nbreak[j,1],q,z,y,t(reparv[j,seq(1,nbreak[j,1]),drop=FALSE]),
                  robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)
          }
        }
        j = j+1
        if (j <= 4) {
          if (nbreak[j,1] == nbreak[j-1,1]){
            if (nbreak[j,1] == 0){
              print(paste('For the',siglev[j,1], '% level, the model is the same as for the',
                          siglev[j-1,1], '% level.'))
              print('The estimation is not repeated')
              ii = 1
            }
            else {
              if (identical(dateseq[j,seq(1,nbreak[j,1]),drop = FALSE],dateseq[j-1,seq(1,nbreak[j-1,1]),drop = FALSE])){
                print(paste('For the',siglev[j,1], '% level, the model is the same as for the',
                            siglev[j-1,1], '% level.'))
                print('The estimation is not repeated')
                ii = 1
              }
            }
          }
          else {
            ii = 0
          }
        }
      }
    }

    if('fix' %in% method){
      t_out = doglob(y,z,x,fixn,eps,eps1,maxi,fixb,betaini,printd)
      datevec = t_out$datevec
      if(length(datevec) == 0){
        print('No break is found')
      }else{
      print(paste('Output from the estimation of the model with',fixn,'breaks'))
      out$est_fix = estim(fixn,q,z,y,datevec[,fixn,drop=FALSE],robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)}
      }
  }
  else{print('estimation methods are skipped')}
  return(out)
}

