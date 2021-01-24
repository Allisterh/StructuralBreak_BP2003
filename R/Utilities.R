#' Utilities functions
#'

###### Matrix auxiliary functions (only available in MATLAB) ######

#' Rotate a matrix
#'
#'Function rotate sa matrix 90 degree counterclockwise
#'@author Linh Nguyen
#'@param A matrix to be rotated 90 degree counterclockwise
#'@return B: rotated matrix A
#'@examples
#'A = matrix(c(1:10),5,2)
#'rot90(A)
#'@export
rot90 = function(A){
  rA = dim(A)[1]
  cA = dim(A)[2]
  B = matrix (0L, nrow = cA, ncol = rA)
  for (i in 1:rA){
    tempA = A[i,,drop=FALSE]
    tempA = tempA[,seq(dim(A)[2],1,-1)] #reverse the row elements
    tempB = t(tempA) #transpose the row
    B[,i] = tempB
  }
  return(B)
}

#'Kronecker Tensor Product
#'
#'Function computes Kronecker Tensor Product between matrix A and matrix B
#'
#'@author Linh Nguyen
#'@param A Matrix A
#'@param B Matrix B
#'@return out: Kronecker Tensor Product of A (X) B
#'@examples
#'A = matrix(c(1,2),1,2)
#'B = matrix(c(1:4),2,2)
#'kron(A,B)
#'@export
kron = function(A,B){
  ###
  # A,B = m x n matrix A & p x q matrix B
  # ****
  # out = A (x) B
  ###
  rA = dim(A)[1]
  cA = dim(A)[2]
  rB = dim(B)[1]
  cB = dim(B)[2]

  out = matrix(0L, nrow = rA * rB, ncol = cA * cB)

  for (i in 1:rA){
    for (j in 1:cA) {
      r_index = seq((i-1)*rB+1,i*rB)
      c_index = seq((j-1)*cB+1,j*cB)
      out[r_index,c_index] = A[i,j] * B
    }
  }
  return(out)
}

#'Diagonal partition given break dates
#'
#' Function constructs the matrix of regressors z which coefficients are changed
#' on the estimated break dates
#'
#'
#'@param input matrix of independent variables z with coefficients allowed to
#'change overtime
#'@param m number of breaks in the series
#'@param date vector of estimated break dates
#'@examples
#' z = matrix(c(1:100),50,2)
#' m = 2  #2 breaks
#' date = matrix(c(15,30),2,1) #first break at t = 15; second break at t = 30
#' diag_par(z,m,date)
#'
#'@return output: matrix of partitioned variables corresponds to break dates
#'@export
###
diag_par = function(input,m,date) {

  nt = dim(input)[1]
  q1 = dim(input)[2]
  #create output matrix of zeros with size: nt x (break+1)*q1
  output = matrix(0L,nrow = nt, ncol = (m+1)*q1)
  #copy values of 1st segment of z to the output matrix
  output[c(1:date[1,1]),c(1:q1)] = input[c(1:date[1,1]),,drop=FALSE]
  i = 2
  while (i <= m){
    #copy values of i-th segment of z to output matrix corresponding to date vector
    r_i = seq(date[i-1,1]+1,date[i,1],1) #indices of rows to copy input values
    c_i = seq((i-1)*q1+1,i*q1,1)     #indices of cols to copy input values
    output[r_i,c_i]=input[r_i,,drop=FALSE]
    i = i+1
  }
  rl_i = seq(date[m,1]+1,nt,1)
  c_i = seq(m*q1+1,(m+1)*q1,1)
  output[rl_i,c_i] = input[rl_i,,drop=FALSE]
  return (output)
}


######Computation auxiliary functions ############
#' OLS regression in matrix form
#'
#' Function computes OLS estimates of the regression y on x
#' in matrix form, \eqn{\beta = (X'X)^{-1} X'Y}
#'
#'@param y matrix of dependent variables
#'@param x matrix of independent variables
#'@return b: OLS estimates of the regression
#'
#'@export
OLS = function(y,x) {
  b = solve((t(x) %*% x)) %*% t(x) %*% y
  b = matrix(b)
  return (b)
}


#' SSR computation
#'
#' Function computes Sum of Squared Residuals based on OLS estimates, \eqn{\beta}
#'
#'@import olsqr
#'@param y matrix of dependent variables
#'@param zz matrix of independent variables
#'@return SSR: Sum of Squared Residuals
#'@export
ssr = function(y,zz) {
  delta = OLS(y,zz)
  resid = y - zz %*% delta
  ssr = t(resid) %*% resid
  return(ssr)
}

