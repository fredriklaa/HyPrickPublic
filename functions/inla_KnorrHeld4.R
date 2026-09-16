##
# @Analysis of the Knorr-Held-IV model using inla
#
# @param  dataST  A data frame containing (at least) for variables
#  Observations (response), a numeric variable
#  A factor describing spatial locations
#  A factor describing temporal locations
#  A factor describing deltas
# @param  indd A vector of length four giving the columns for spatial, temporal, interaction and response
# @param  Qtemp Precision matrix for temporal component
# @param  Qspat Precision matrix for spatial component
# @param  family  The likelihood family. 
# @param  constr  Type of constraint, either "SC" or " "  
# @param  method  Method for performing computation. Either "standard", "hymik", "hybw" or "hyprick"
# @param  prec.intercept   Precision for intercept, default 0.001
# @param  diagonal  Diagonal element added to precision matrices to make them non-singular
# @param  kappa The kappa parameter in the HyMik1 and HyPrick methods
#
# @return Description of the return value.
#
# @details
# Longer description of what the function does, any assumptions,
# performance notes, references, etc.
#
# @note   Optional extra notes.
# @warning Optional warnings or caveats.
#
# @examples
##
inla_KnorrHeld4 = function(dataST,Qtemp,Qspat,indd=c(1,2,3,4),extracov="",family="gaussian",scale=FALSE,
                           constr="SC",method="standard",nthreads = "12:2",prec.intercept=0,diagonal=1e-05,kappa=1e06)
{
  #Convert names to "standard" names
  names(dataST)[indd[1]] = "alpha"        # Temporal component
  names(dataST)[indd[2]] = "gamma"        # Spatial component
  names(dataST)[indd[3]] = "delta"        # Interaction compoment
  names(dataST)[indd[4]] = "Y"            # Response


  
  ns=nrow(Qspat)
  nt=nrow(Qtemp)
  eps = diagonal
  
  d = c(1:nt)-mean(1:nt)
  dtilde = d/sqrt(sum(d^2))
  Qt11 = Qtemp[1,1]
  Qs11 = Qspat[1,1]
  kappast = kappa
  if(scale)
  {
    #Scale models
    #scf1= exp(mean(log(diag(MASS::ginv(as.matrix(Qtemp))))))
    #Qtemp = scf1*Qtemp
    Qtemp = inla.scale.model(Qtemp,list(A=rbind(rep(1/sqrt(nt),nt),dtilde),e=c(0,0)))
    #scf2 = exp(mean(log(diag(MASS::ginv(as.matrix(Qspat))))))
    #Qspat = scf2*Qspat
    Qspat= inla.scale.model(Qspat,constr=list(A=matrix(rep(1/sqrt(ns),ns),nrow=1),e=0))
    scf1= Qtemp[1,1]/Qt11
    scf2 = Qspat[1,1]/Qs11
    kappast = kappa*scf1*scf2
    print(paste0("Scaled Qtemp by ",sf1,"Scaled Qspat by ",sf2))
  }
  #Precision matrix for interaction term
  Q_st=kronecker(Qtemp,Qspat)
  prior.fixed=list(prec.intercept=prec.intercept)
  
  
  #Compute constraints
  Ins = Diagonal(ns)
  Int = Diagonal(nt)
  Ifull = Diagonal(ns*nt)
  Atime = kronecker(matrix(rep(1/sqrt(nt),nt),nrow=1),Ins)
  Atime = Atime[-nrow(Atime),]
  Atime2 = kronecker(matrix(dtilde,nrow=1),Ins)
  Atime2 = Atime2[-nrow(Atime2),]       ## Explain this
  Aspace = kronecker(Int,matrix(rep(1/sqrt(ns),ns),nrow=1))
  if(ns>nt)
  {
    A1 = Atime
    if(constr=="SC")
    {
      A1 = rbind(A1, Atime2)
    }
    A2=Aspace
  }
  else
  {
    A1 = Aspace
    A2 = Atime
    if(constr=="SC")
      A2 = rbind(A2,Atime2)
  }
  #Remove singular constraint from A2.
  #A2 = A2[-nrow(A2),]
  A = rbind(A1,A2)
  Z = Ifull-Matrix::crossprod(A1)

  #Formula including covariates and main effects
  extracov = paste0(extracov,'+')
  form = paste0('Y~',extracov,'f(gamma,model="generic0",diagonal=eps,Cmatrix=Qtemp,constr=T)+f(alpha,model="generic0",diagonal=eps,Cmatrix=Qspat,constr=T)')
  
  if(method=="standard")
  {
    form = as.formula(paste0(form,'+f(delta,model="generic0",Cmatrix=Q_st,diagonal=eps,constr=F,extraconstr = list(A=A,e=rep(0,nrow(A))))'))
    resinla =inla(form,data=dataST,verbose=T,family=family,control.fixed=prior.fixed,num.threads=nthreads)
  }
  if(method=="hymik")
  {
    AHyMiK = A2%*%Z
    form = as.formula(paste0(form,'+f(delta,model="z",diagonal=eps,precision=kappa,Z=Z,Cmatrix = Q_st,constr=F,extraconstr = list(A=cbind(AHyMiK,AHyMiK*0),e=rep(0,nrow(AHyMiK))))'))
    resinla =inla(form,data=dataST,verbose=T,family=family,control.fixed=prior.fixed,num.threads=nthreads)
  }
  if(method=="hyprick")
  {
   form = as.formula(paste0(form,'+f(delta,model="generic0",diagonal=eps,Cmatrix=Q_st+kappast*(Ifull-Z),constr=F,extraconstr=list(A=A2,e=rep(0,nrow(A2))))'))
    resinla = inla(form,data=dataST,verbose=T,family=family,control.fixed=prior.fixed,num.threads=nthreads)
  }
  if(method=="hybw")
  {
    A_Bolin=A1
    TMat=c_basis2(A_Bolin)
    
    QT=(Q_st)%*%(TMat$T)
    QT=(Q_st)%*%(TMat$T)
    QT=drop0(QT,tol=0)
    TQT=(t(TMat$T)%*%QT)
    TQT=drop0(TQT,tol=0)
    
    C=1:nrow(A_Bolin)
    U=(1+nrow(A_Bolin)):(ncol(A_Bolin))
    
    A2Z = A2%*%Z
    
    form = as.formula(paste0(form,'+f(delta,model="z",diagonal=eps,Z=TMat$T[,U],precision = kappa,Cmatrix=TQT[U,U],constr=F,extraconstr=list(A=cbind((A2Z),matrix(0,nrow=nrow(A2Z),ncol=length(U))),e=rep(0,nrow(A2Z))))'))
    resinla = inla(form,data=dataST,verbose=T,family=family,control.fixed=prior.fixed,num.threads=nthreads)
  }
  if(method=="bw")
  {
    #rcpp called.
    Rcpp::sourceCpp(paste0(dir,"cp.cpp"))
    A_Bolin=A
    TMat=c_basis2(A_Bolin)
    
    QT=(Q_st)%*%(TMat$T)
    QT=(Q_st)%*%(TMat$T)
    QT=drop0(QT,tol=0)
    TQT=(t(TMat$T)%*%QT)
    TQT=drop0(TQT,tol=0)
    
    C=1:nrow(A_Bolin)
    U=(1+nrow(A_Bolin)):(ncol(A_Bolin))
    
    form = as.formula(paste0(form,'+f(delta,model="z",diagonal=0,Z=TMat$T[,U],precision = kappa,Cmatrix=TQT[U,U]+eps*Diagonal(length(U)),constr=F)'))
    resinla = inla(form,data=dataST,verbose=T,family=family,control.fixed=prior.fixed,num.threads=nthreads)
  }
  if(method=="bw_A_formulation")
  {
    #rcpp called.
    Rcpp::sourceCpp(paste0(dir,"cp.cpp"))
    print("Note that extracov is *not* supported!")
    if(length(extracov)>1){
      print("Not supported!")
      return(0)
    }
    A_Bolin=A
    #browser()
    TMat=c_basis2(A_Bolin)
    
    QT=(Q_st)%*%(TMat$T)
    QT=drop0(QT,tol=0)
    TQT=(t(TMat$T)%*%QT)
    TQT=drop0(TQT,tol=0)
    
    C=1:nrow(A_Bolin)
    U=(1+nrow(A_Bolin)):(ncol(A_Bolin))
    
    Xalpha=sparse.model.matrix(~-1+as.factor(alpha),dataST)
    Xgamma=sparse.model.matrix(~-1+as.factor(gamma),dataST)
    Xdelta=sparse.model.matrix(~-1+as.factor(delta),dataST)
    Xtot=cbind(1,Xgamma,Xalpha,Xdelta%*%TMat$T[,U])
    
    form = paste0('Y~-1+Intercept+','f(gamma,model="generic0",diagonal=eps,Cmatrix=Qtemp,constr=T)+f(alpha,model="generic0",diagonal=eps,Cmatrix=Qspat,constr=T)')
    
    dataST_list=list(Y=dataST$Y,Intercept=c(1,rep(NA,ns+nt+length(U))),gamma=c(NA,1:nt,rep(NA,ns+length(U))),alpha=c(rep(NA,1+nt),1:ns,rep(NA,length(U))),delta=c(rep(NA,1+ns+nt),1:length(U)))
    form = as.formula(paste0(form,'+f(delta,model="generic0",diagonal=0,Cmatrix=TQT[U,U]+eps*Diagonal(length(U)),constr=F)'))
    resinla = inla(form,data=dataST_list,verbose=T,family=family,control.fixed=prior.fixed,num.threads=nthreads,control.predictor=list(A=Xtot))
  }
  resinla
}


BWAlg1 = function(A)
{
  n = ncol(A)
  k = nrow(A)
  Tmat = Diagonal(n)
  Asum = colSums(A!=0)
  D = c(1:n)[Asum>0]
  Asvd = svd(A)
}
