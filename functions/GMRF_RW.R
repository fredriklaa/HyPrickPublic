
#'Creates the structure matrix for random walk order 1 (RW1) and RW2
#'The precision matrix is not scaled
#'
#'@param n gives the number of time steps
#'@param order gives the order of the random walk
#'@return a matrix with the precision matrix
#'@export
#'@name GMRF_RW
GMRF_RW <- function(n = 10,order=1) {
  precinc = 1
  stopifnot(order %in% c(1,2))

  if(is.null(n)) n<- nrow(mu)

  mu <- matrix(0,nrow=n,ncol=1)
  i = c(1:(n-1),1:(n-1))
  j = c(1:(n-1),2:n)
  x <- numeric(length=((n-1)*2))
  x[1:(n-1)] = -1
  x[n:((2*n)-2)] = 1
  Dmat = Matrix::sparseMatrix(i,j,x=x)
  R = Matrix::t(Dmat)%*%Dmat
  if (order == 1) {
    Q = R
  }
  if (order == 2) {
    R <- R %*% R
    R[1,(1:3)] = c(1,-2,1)
    R[2,(1:4)] = c(-2,5,-4,1)
    R[(n-1),(n-3):n] = c(1,-4,5,-2)
    R[(n),(n-2):n] = c(1,-2,1)
    Q = R
  }

  return(Q)
}
