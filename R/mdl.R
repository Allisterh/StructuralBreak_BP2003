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

mdl <- function(y_name,z_name = NULL,x_name = NULL,data,eps1 = 0.15,m = -1,prewhit = 1,
                           robust = 1,hetdat = 1,hetvar = 1,hetomega = 1,hetq = 1,
    procedure = c('dotest'),
    CI = 1,maxi = 10,fixb = 0,
    eps = 0.00001,betaini = 0,fixn=-1,printd = 0){
  #handle data
  y_ind = match(y_name,colnames(data))
  x_ind = match(x_name,colnames(data))
  z_ind = match(z_name,colnames(data))
  
  cat('The options chosen are:\n')
  cat(paste(' i) hetdat = ',hetdat),'\n',paste('ii) hetvar = ',hetvar),'\n',
      paste('iii) hetomega = ',hetomega),'\n',paste('iv) hetq = ',hetq),'\n',
      paste('v) robust = ',robust),'\n',paste('vi) prewhite = ',prewhit),'\n')
  
  mdl = list()

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


  

  if('dotest' %in% procedure){
   mdl$Wtest = dotest(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,robust,
           hetdat,hetvar)
  }

  if('dospflp1' %in% procedure){
    mdl$spflp1 = dospflp1(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
             robust,hetdat,hetvar)
  }
  
  if('doorder' %in% procedure){
    mdl$BIC = doorder(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,1)
    mdl$LWZ = doorder(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,0)
  }
  
  if('dosequa' %in% procedure){
   mdl$sequa = dosequa(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
            robust,hetdat,hetvar)
  }
  
  if('dorepart' %in% procedure){
    mdl$repart = dorepart(y,z,x,m,eps,eps1,maxi,fixb,betaini,printd,prewhit,
             robust,hetdat,hetvar)
  }
  
  if('fix' %in% procedure){
    t_out = doglob(y,z,x,fixn,eps,eps1,maxi,fixb,betaini,printd)
    datevec = t_out$datevec
    if(length(datevec) == 0){
      print('No break is found')
      return(NULL)
    }else{
      #  print(paste('Output from the estimation of the model with',fixn,'breaks'))
    date = datevec[,fixn,drop=FALSE]
    fix_mdl = estim(fixn,q,z,y,date,robust,prewhit,hetomega,hetq,x,p,hetdat,hetvar)}
    mdl$fix = fix_mdl
    mdl$fix$p_name = 'fix'
    mdl$fix$nbreak = fixn
    class(mdl$fix) = 'model'
    mdl$fix$numz = q
    mdl$fix$numx = p
    mdl$fix = compile.model(mdl$fix)
  }



  #reorganize the results into the list
  class(mdl) <- 'mdl'
  mdl$maxb = m
  mdl$procedure = procedure
  
  return(mdl)
  }



print.mdl <- function(x,digits = -1,...)
{
  proc = x$procedure
  cat(paste('\nProcedures invoked for maximum',x$maxb,'breaks:\n\n'))
  for (p_name in proc){
    cat(paste(p_name,'\n'))
  }
  cat(paste('\nTo obtain information about specific procedure, 
  type stored variable name + \'$\' + procedure name'))
  
  invisible(x)
}






