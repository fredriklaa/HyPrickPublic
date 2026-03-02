library(Matrix)
library(INLA)
#Can include pardiso if available
#INLA::inla.setOption("pardiso.license","xxx/sys/licences/pardiso.lic")
#INLA::inla.pardiso.check()
library(data.table)
library(ggplot2)
library(xtable)
library(Rcpp)

rm(list=ls())
dir=""
dirfun = paste0(dir,"/functions/")
dirdata=paste0(dir,"/data/")

source(paste0(dirfun,"SpaceTimeProjConstr.R"))
source(paste0(dirfun,"GMRF_RW.R"))

#rcpp called.
Rcpp::sourceCpp(paste0(dirfun,"cp.cpp"))

#Argument to INLA:
nthreads = "3:2"


graph=system.file("demodata/germany.graph", package="INLA")


#Read data:

df = readRDS(paste0(dirdata,"/SpatioTemporalData.RDS"))

Q_ICAR=INLA::inla.graph2matrix(graph)
diag(Q_ICAR)=0
Q_ICAR=-Q_ICAR
diag(Q_ICAR)=-rowSums(Q_ICAR)
ns=nrow(Q_ICAR)

nt=20

Q_RW2=GMRF_RW(n=nt,order=2)

Q_st=kronecker(Q_RW2,Q_ICAR)
prior.fixed=list(prec.intercept=0.001)

eps=1e-05
kap=1e06

A_Bolin=kronecker(rep(1,nt),Diagonal(ns))

tstar = c(1:nt)
tstar = tstar/sqrt(sum(tstar^2))
B1 = kronecker(tstar,Diagonal(ns))
A_Bolin=cbind(B1,A_Bolin)

TMat=c_basis2(t(A_Bolin))
eps=1e-05

QT=(Q_st+Diagonal(ns*nt)*eps)%*%(TMat$T)
QT=drop0(QT,tol=0)
TQT=(t(TMat$T)%*%QT)
TQT=drop0(TQT,tol=0)

C=1:nrow(t(A_Bolin))
U=(1+nrow(t(A_Bolin))):(ncol(t(A_Bolin)))

PC =SpaceTimeProjConstr(ns = ns,nt=nt,type ="SC",dim = "space")

A_BB=kronecker(Diagonal(nt),rep(1,ns))
A_BB2=t(PC$P%*%A_BB)
A_newHymiK=t(PC$P%*%(A_BB)[,-1])


A_newB=t(PC$P%*%t(PC$A2[,seq(1,ns*nt)]))
A_new_BOLIN=A_newB%*%t(TMat$T)
A_new_BOLIN_sub=A_new_BOLIN[,U]

A_old=t(A_BB)[-1,]
A_BB_alternative=t(t(A_BB)%*%TMat$T[,U])

library(INLA)

fil = "sim.Standard.RDS"
if(file.exists(fil))
  StandardINLA = readRDS(fil)
if(!file.exists(fil))
{
  StandardINLA =inla(Y~f(main_temporal,diagonal=0,model="generic0",Cmatrix=Q_RW2+Diagonal(nt)*eps,constr=T)+
                       f(main_spatial,diagonal=0,model="generic0",Cmatrix=Q_ICAR+Diagonal(ns)*eps,constr=T)+
                       f(interaction,model="generic0",Cmatrix=Q_st+Diagonal(ns*nt)*eps,constr=F,
                         extraconstr = list(A=PC$A,e=rep(0,nrow(PC$A)))),
                     data=df,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=nthreads)
  #control.predictor=list(compute=TRUE
  saveRDS(StandardINLA,file=fil)
}

fil = "sim.HyMiK.RDS" 
if(file.exists(fil))
  HyMiK = readRDS(fil)
if(!file.exists(fil))
{
  HyMiK =inla(Y~f(main_temporal,diagonal=0,model="generic0",Cmatrix=Q_RW2+Diagonal(nt)*eps,constr=T)+
                f(main_spatial,diagonal=0,model="generic0",Cmatrix=Q_ICAR+Diagonal(ns)*eps,constr=T)+
                f(interaction,model="z",diagonal=0,precision=1e5,Z=PC$P,Cmatrix = Q_st+Diagonal(ns*nt)*eps,constr=F,
                  extraconstr = list(A=cbind(A_newHymiK,A_newHymiK*0),e=rep(0,nrow(A_newHymiK)))),
              data=df,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=nthreads)
  #control.predictor=list(compute=TRUE
  saveRDS(HyMiK,file=fil)
}

fil = "sim.HyPrick.RDS"
if(file.exists(fil))
  HyPrick = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick = inla(Y~f(main_temporal,diagonal=0,model="generic0",Cmatrix=Q_RW2+Diagonal(nt)*eps,constr=T)+
                   f(main_spatial,diagonal=0,model="generic0",Cmatrix=Q_ICAR+Diagonal(ns)*eps,constr=T)+
                   f(interaction,model="generic0",diagonal=0,Cmatrix=Q_st+10^5*(Diagonal(nrow(Q_st))-PC$P)+Diagonal(nrow(Q_st))*eps,constr=F,extraconstr=list(A=A_old,e=rep(0,nrow(A_old)))),
                 data=df,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=nthreads)
  saveRDS(HyPrick,file=fil)
} 

fil = "sim.bolinWallin.RDS"
if(file.exists(fil))
  BolinWallin = readRDS(fil)
if(!file.exists(fil))
{
  BolinWallin = inla(Y~f(main_temporal,diagonal=0,model="generic0",Cmatrix=Q_RW2+Diagonal(nt)*eps,constr=T)+
                       f(main_spatial,diagonal=0,model="generic0",Cmatrix=Q_ICAR+Diagonal(ns)*eps,constr=T)+
                       f(interaction,Z=TMat$T[,U],precision = kap,model="z",diagonal = 0,Cmatrix=TQT[U,U]+eps*Diagonal(length(U)),
                         constr=F,extraconstr=list(A=cbind((A_newB),(A_new_BOLIN_sub)*0),e=rep(0,-1+ncol(A_BB)))),
                     data=df,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=nthreads)
  saveRDS(BolinWallin,file=fil)
}

tab =cbind(StandardINLA$summary.hyperpar$mean,HyMiK$summary.hyperpar$mean,
           BolinWallin$summary.hyperpar$mean,HyPrick$summary.hyperpar$mean)
tab[1,] = 1/tab[1,]
tab = rbind(c(StandardINLA$summary.fixed$mean,HyMiK$summary.fixed$mean,
              BolinWallin$summary.fixed$mean,HyPrick$summary.fixed$mean),
            tab,c(StandardINLA$cpu.used[4],HyMiK$cpu.used[4],BolinWallin$cpu.used[4],HyPrick$cpu.used[4]))
colnames(tab) = c("Standard","HyMiK","BW","HyPrick")
rownames(tab)=c("Mean","Overdispersion","Precision Temporal","Precision Spatial","Precision Interaction","CPU")

library(ggplot2)
library(patchwork)

plotDataE=data.table::data.table(StandardINLAE=StandardINLA$summary.random$interaction$mean,
                                 HyPrickE=HyPrick$summary.random$interaction$mean,
                                 HyMiKE=HyMiK$summary.random$interaction$mean[1:(nt*ns)],
                                 BWE=BolinWallin$summary.random$interaction$mean[1:(nt*ns)])

plotDataSd=data.table::data.table(StandardINLAE=StandardINLA$summary.random$interaction$sd,
                                  HyPrickE=HyPrick$summary.random$interaction$sd,
                                  HyMiKE=HyMiK$summary.random$interaction$sd[1:(nt*ns)],
                                  BWE=BolinWallin$summary.random$interaction$sd[1:(nt*ns)])

p1 = ggplot(data = plotDataE) +
  geom_point(aes(y = HyPrickE, x = StandardINLAE, colour = "HyPrick"), size = 1) +
  geom_point(aes(y = HyMiKE,   x = StandardINLAE, colour = "HyMiK"),   size = 1) +
  geom_point(aes(y = BWE,      x = StandardINLAE, colour = "BW"),      size = 1) +
  scale_colour_manual(
    name = "Method",
    values = c("HyPrick" = "black", "HyMiK" = "red", "BW" = "blue")
  ) +
  geom_abline(intercept = 0, slope = 1, size = 0.5,color="green") +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  coord_fixed() +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm")) +
  labs(x = "Standard", y = "Alternatives")

p2 = ggplot(data = plotDataSd) +
  geom_point(aes(y = HyPrickE, x = StandardINLAE, colour = "HyPrick"), size = 1) +
  geom_point(aes(y = HyMiKE,   x = StandardINLAE, colour = "HyMiK"),   size = 1) +
  geom_point(aes(y = BWE,      x = StandardINLAE, colour = "BW"),      size = 1) +
  scale_colour_manual(
    name = "Method",
    values = c("HyPrick" = "black", "HyMiK" = "red", "BW" = "blue")
  ) +
  geom_abline(intercept = 0, slope = 1, size = 0.5,color="green") +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  coord_fixed() +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm")) +
  labs(x = "Standard", y = "Alternatives")

(p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave("Sim_Interaction_E_sd.pdf")