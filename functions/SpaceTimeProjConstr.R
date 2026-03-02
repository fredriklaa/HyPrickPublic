#' Create the projector matrix and the corresponding constraints
#' 
#' @param 'ns' The number of spatial points
#' @param 'nt' The number of temporal points
#' @param 'type' Either 'GCF' in which case the Goicoa (GCF) constraints are used or 
#' 'SC' in which case the additional SchrC8dle (SC) constraints are used. Default is 'GCF'
#' @param 'dim' Either 'space' in which case the constraints with respect to space are used, or 
#' 'time' in which case the constraints with respect to time is used. If no specification is given, 
#' 'space' is given if ns>nt and 'time' in the opposite case.  
#' 
#' @return A list containing the projection matrix 'P' and the remaining constraints 'C', both of class 'dgCMatrix'.
#'
#' @export
SpaceTimeProjConstr = function(ns,nt,type="GCF",dim="DEFAULT")
  {
  require(Matrix)
  if(dim=="DEFAULT")
  {
    if(ns>nt)
      dim = "space"
    else
      dim = "time"
  }
  r = list()
  if(dim=="space")
  {
    A1 = kronecker(matrix(rep(1,nt),nrow=1),Diagonal(ns))
    A2=kronecker(Diagonal(nt),matrix(rep(1,ns),nrow=1))
    A2 = A2[-nrow(A2),]
    P = Diagonal(ns*nt)-Matrix::crossprod(A1)/nt
    P1 = P
    A = rbind(A1,A2)
    B1 = NULL
    if(type=="SC")
    {
      tstar = c(1:nt)-mean(c(1:nt))
      tstar = tstar/sqrt(sum(tstar^2))
      B1 = kronecker(matrix(tstar,nrow=1),Diagonal(ns))
      B1 = B1[-nrow(B1),]
      P = P - Matrix::crossprod(B1)
      A = rbind(A,B1)
    }
    
    return(list(P=P,A2=cbind(A2,0*A2),A=A,
                A1=A1,B1=B1,P1=P1))
  }
  else if(dim=="time")
  {
    A1 = kronecker(Diagonal(nt),matrix(rep(1,ns),nrow=1))
    A2=kronecker(matrix(rep(1,nt),nrow=1),Diagonal(ns))
    P = Diagonal(ns*nt)-Matrix::crossprod(A1)/ns
    A2 = A2[-nrow(A2),]
    A = rbind(A1,A2)
    
    if(type=="SC")
    {
      tstar = c(1:nt)
      tstar = tstar/sqrt(sum(tstar^2))
      B1 = kronecker(matrix(tstar,nrow=1),Diagonal(ns))
      B1 = B1[-nrow(B1),]
      A2=rbind(A2,B1)
      A = rbind(A,B1)
    }
    
    return(list(P=P,A2=cbind(A2,0*A2),A=A))
  }
  else
  {
    print("Unkown dim")
    return(NULL)
  }
}

