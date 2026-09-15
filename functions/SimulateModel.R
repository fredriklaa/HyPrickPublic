library(Matrix)
library(INLA)
library(sparseMVN)
#library(INLAconstraints)
library(data.table)
#Number of time points
dir = "functions/"
#Update:
source(paste0(dir,"SpaceTimeProjConstr.R"))
source(paste0(dir,"GMRF_RW.R"))



set.seed(20)
nt=10
graph=system.file("demodata/germany.graph", package="INLA")
I=INLA::inla.graph2matrix(graph)
diag(I)=0
I=-I
diag(I)=-rowSums(I)

I2=I
ns=nrow(I2)
#Small value, so that we get proper precision matrix
#We first sample from a proper distribution with too large variances in some directions, and later correct these variances. 
epsilon=1e-05
#Set up constraints
A_t = matrix(1,1,nt)
A_t=rbind(A_t,seq(1,nt))
A_s=matrix(1,ncol=ns,nrow=1)
A_st=cbind(kronecker(rep(1,nt),diag(ns)),kronecker(diag(nt),rep(1,ns))[,-1],kronecker(seq(1,nt),diag(ns))[,-1])
A_st=t(A_st)
#Set up precision matrix
Q_RW2_before=GMRF_RW(n=nt,order=2)
Q_RW2 = Q_RW2_before
#Scale model
#Q_RW2 = inla.scale.model(Q_RW2_before,list(A=A_t,e=rep(0,2)))

#Q_RW2_scaled=inla.scale.model(GMRF_RW(order=2,n=10),constr = list(A=rbind(rep(1,nt),seq(1,nt)),e=rep(0,2)))
Q_ICAR = I2
#Q_ICAR=inla.scale.model(I2,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
Q_st=kronecker(Q_RW2,Q_ICAR)

#Compute cholesky factorization
CH_RW2=Cholesky(Q_RW2+diag(nt)*epsilon)
CH_ICAR=Cholesky(Q_ICAR+diag(ns)*epsilon)
#CH_Qst = kronecker(CH_RW2,CH_ICAR)
CH_Qst=Cholesky(Q_st+diag(ns*nt)*epsilon)

#Simulate from field
x_t_main=sparseMVN::rmvn.sparse(1,mu=rep(0,nt),CH=CH_RW2,prec=TRUE)
x_s_main=sparseMVN::rmvn.sparse(1,mu=rep(0,ns),CH=CH_ICAR,prec=TRUE)
x_st_main=sparseMVN::rmvn.sparse(1,mu=rep(0,ns*nt),CH=CH_Qst,prec=TRUE)
#Correct for constraints: Conditioning by kriging: 
x_t_main=as.vector(x_t_main)-solve(Q_RW2+diag(nt)*epsilon)%*%t(A_t)%*%solve(A_t%*%solve(Q_RW2+diag(nt)*epsilon)%*%t(A_t))%*%A_t%*%as.vector(x_t_main)
x_s_main=as.vector(x_s_main)-solve(Q_ICAR+diag(ns)*epsilon)%*%t(A_s)%*%solve(A_s%*%solve(Q_ICAR+diag(ns)*epsilon)%*%t(A_s))%*%A_s%*%as.vector(x_s_main)
solveQst=solve(Q_st+diag(ns*nt)*epsilon)
x_st_main=as.vector(x_st_main)-solveQst%*%t(A_st)%*%solve(A_st%*%solveQst%*%t(A_st))%*%A_st%*%as.vector(x_st_main)




#NB: we must add the slope to get a simple sum to zero constraint. 
slope = rnorm(1,mean=0,sd=2)
slopes_st=slope + rnorm(ns,mean=0,sd=0.1)

main=seq(1,ns*nt)%%ns
main[main==0]=ns
temporal=ceiling(seq(1,ns*nt)/ns)

data=data.table(main_spatial=main,main_temporal=temporal,interaction=seq(1,ns*nt))
dataF=copy(data)
dataF$main_spatial=factor(dataF$main_spatial)
dataF$main_temporal=factor(dataF$main_temporal)
dataF$interaction=factor(dataF$interaction)

DesignMatrixMainSpatial=model.matrix(~main_spatial-1,dataF)
DesignMatrixMainTemporal=model.matrix(~main_temporal-1,dataF)
DesignMatrixMainInteraction=model.matrix(~interaction-1,dataF)

NewConstraints=kronecker(seq(1,nt),diag(ns))
NewConstraintsScaled=NewConstraints/sqrt(colSums(NewConstraints^2))
NewConstraintsScaledProjectorMatrix=diag(ns*nt)-(NewConstraintsScaled)%*%t(NewConstraintsScaled)


mean_pred=DesignMatrixMainTemporal%*%as.vector(x_t_main)+DesignMatrixMainSpatial%*%as.vector(x_s_main)+DesignMatrixMainInteraction%*%as.vector(x_st_main)
mean_pred=mean_pred+slope*data$main_temporal/nt
for(i in 1:(nt)){
  print(i)
  mean_pred=mean_pred+slopes_st[i]*DesignMatrixMainInteraction%*%NewConstraintsScaled[,i]
}



combined_pred=do.call(rbind,replicate(30,mean_pred,simplify=FALSE))
data=do.call(rbind,replicate(30,data,simplify = FALSE))

set.seed(20)
mu = 1
data$Y=rpois(n=nrow(combined_pred),lambda=exp(-1+combined_pred))
show(range(data$Y))

SpatioTemporalDataset = data
saveRDS(SpatioTemporalDataset,file="data/SpatioTemporalData.RDS")

#saveRDS(data,file="D:/SpatioTemposave(ralDataset.RDS")
m = nt
min_limit_y=min(c(x_t_main[,1]+slope*seq(1,m)/m,x_t_main[,1]))
max_limit_y=max(c(x_t_main[,1]+slope*seq(1,m)/m,x_t_main[,1]))


#Create plot, using base graphics. Black: simple sum to zero, while blue is "fully" constrained field. 
plot(x_t_main+slope*seq(1,m)/m,ylim=c(-1+min_limit_y,max_limit_y+1),type="l")
points(x_t_main,col="blue",type="l")



