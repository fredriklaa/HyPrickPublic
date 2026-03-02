library(Matrix)
library(INLA)
library(sparseMVN)
library(INLAconstraints)
library(data.table)
#Number of time points


dirdata=""
set.seed(20)
nt=20
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
A_st=cbind(kronecker(rep(1,nt),Diagonal(ns)),kronecker(Diagonal(nt),rep(1,ns))[,-1],kronecker(seq(1,nt),Diagonal(ns))[,-1])
A_st=t(A_st)
#Set up precision matrix
Q_RW2_before=GMRF_RW(n=nt,order=2)
#Scale model
Q_RW2 = inla.scale.model(Q_RW2_before+Diagonal(nt)*epsilon,list(A=A_t,e=rep(0,2)))
Q_ICAR=inla.scale.model(I2+Diagonal(ns)*epsilon,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
Q_st=kronecker(Q_RW2,Q_ICAR)

#Simulate from field
x_t_main_init=inla.qsample(n=1,Q=Q_RW2+Diagonal(nt)*epsilon,seed = 1)
x_s_main_init=inla.qsample(n=1,Q=Q_ICAR+Diagonal(ns)*epsilon,seed = 1)
x_st_main_init=inla.qsample(n=1,Q=Q_st+Diagonal(ns*nt)*epsilon,seed = 1)

#Correct for constraints: Conditioning by kriging: 
x_t_main=inla.qsample(n=1,Q=Q_RW2+Diagonal(nt)*epsilon,mu = x_t_main_init,constr = list(A=A_t,e=c(0,0)),seed = 1)
x_s_main_main=inla.qsample(n=1,Q=Q_ICAR+Diagonal(ns)*epsilon,mu = x_s_main_init,constr = list(A=A_s,e=c(0)),seed = 1)
x_st_main=inla.qsample(n=1,Q=Q_st+Diagonal(ns*nt)*epsilon,mu = x_st_main_init,constr = list(A=A_st,e=rep(0,nrow(A_st))),seed = 1)



main=seq(1,ns*nt)%%ns
main[main==0]=ns
temporal=ceiling(seq(1,ns*nt)/ns)

data_in=data.table(main_spatial=main,main_temporal=temporal,interaction=seq(1,ns*nt))
dataF=copy(data_in)
dataF$main_spatial=factor(dataF$main_spatial)
dataF$main_temporal=factor(dataF$main_temporal)
dataF$interaction=factor(dataF$interaction)

DesignMatrixMainSpatial=model.matrix(~main_spatial-1,dataF)
DesignMatrixMainTemporal=model.matrix(~main_temporal-1,dataF)
DesignMatrixMainInteraction=model.matrix(~interaction-1,dataF)

mean_pred=0.1*DesignMatrixMainTemporal%*%as.vector(x_t_main)+0.5*DesignMatrixMainSpatial%*%as.vector(x_s_main)+0.3*DesignMatrixMainInteraction%*%as.vector(x_st_main)
combined_pred_replicate=do.call(rbind,replicate(30,mean_pred,simplify = FALSE))
data=do.call(rbind,replicate(30,data_in,simplify = FALSE))

set.seed(20)
data$Y=rpois(n=length(combined_pred_replicate[,1]),lambda=exp(0.5+combined_pred_replicate[,1]))
show(range(data$Y))

SpatioTemporalData = data
saveRDS(SpatioTemporalData,file=paste0(dirdata,"SpatioTemporalData.RDS"))
