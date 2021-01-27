#Main program to handle users' input

#' Robust structural break date estimation
#'
#' Function executes the main procedures selected by users. There are 7 main
#' procedures including:
#' \code{\link{doglob}} A recursive procedure to obtain global minimizers of SSR,
#' corresponding to the break date.
#' \code{\link{dotest}} Procedure to conduct SupF test of 0 versus m breaks and
#'  Double Max test.
#' \code{\link{dosupflp1}} Procedure to conduct SupF(l+1|l).
#'
#' @param y name of dependent variable
#' @param z name of independent variables which coefficients are allowed to change
#' across regimes
#' @param x name of independent variables which coefficients are constant across
#' regimes
#' @param data the data set we estimate
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
#' @export

date_estimation = function(y,z,x,data,prewhit,robust,hetdat,hetvar,hetomega,hetq,
    procedure = c('doglobal','dotest','dospflp1','doorder','dosequa','dorepart'),
    estdate,estim_method = c('BIC','LWZ','seq','rep','fix')) {

return (0)
}
