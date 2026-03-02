library(Matrix)
library(sf)
library(spdep)
library(ggplot2)
library(xtable)
library(bigDM)
library(Rcpp)
library(dplyr)

rm(list=ls())
dir = "functions/"
#Update:
source(paste0(dir,"SpaceTimeProjConstr.R"))
source(paste0(dir,"GMRF_RW.R"))
#rcpp called.
Rcpp::sourceCpp(paste0(dir,"cp.cpp"))
nthreads = "3:2"


data("Carto_SpainMUN")

data(Data_LungCancer)
Data_LungCancer$ID
carto.nb <- poly2nb(Carto_SpainMUN)
Carto_SpainMUN_connected=Carto_SpainMUN[-2454,]
carto.nb2 <- poly2nb(Carto_SpainMUN_connected)
W=nb2mat(carto.nb2,style = "B")
W=-W
diag(W)=-rowSums(W)


nt=25
combined_data=inner_join(Data_LungCancer,Carto_SpainMUN_connected,by=c("ID"="ID"))
combined_data$year.index=combined_data$year-min(combined_data$year)+1
combined_data$region.index=as.numeric(as.factor(combined_data$ID))
combined_data=combined_data[combined_data$year.index<=nt,]
combined_data$interaction=(combined_data$year.index-1)*length(carto.nb2)+combined_data$region.index

#nt=max(combined_data$year.index)
ns =nrow(W)
Q_RW2=GMRF_RW(n=nt,order=2)
Q_st=kronecker(Q_RW2,W)

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
#TQT=(t(TMat$T)%*%((Q_st+Diagonal(ns*nt)*eps)[U,U]%*%(TMat$T)))
#TQT=drop0(TQT)
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

if(0)
{
  fil = paste0("Largedata_Standard_nt=",nt,".RDS")
if(file.exists(fil))
  StandardINLA = readRDS(fil)
if(!file.exists(fil))
{
  StandardINLA=inla(obs.x~offset(log(exp.x))+
                 f(year.index,model="generic0",diagonal=0,Cmatrix=GMRF_RW(n=nt,order=2)+Diagonal(nt)*1e-05,constr=T)+
                 f(region.index,constr=T,Cmatrix=W+Diagonal(nrow(W))*1e-05,diagonal =0,model="generic0")+
                 f(interaction,model="generic0",diagonal=0,Cmatrix=Q_st+Diagonal(nrow(Q_st))*1e-05,constr=F,
                   extraconstr=list(A=PC$A,e=rep(0,nrow(PC$A)))),
               verbose=T,family="poisson",control.fixed = list(prec.intercept =0.01),
               num.threads=nthreads,data=combined_data)
  saveRDS(StandardINLA,file=fil)
}}

fil = paste0("Largedata_HyPrick_nt=",nt,".RDS")
if(file.exists(fil))
  HyPrick = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick=inla(obs.x~offset(log(exp.x))+
                 f(year.index,model="generic0",diagonal=0,Cmatrix=GMRF_RW(n=nt,order=2)+Diagonal(nt)*1e-05,constr=T)+
                 f(region.index,constr=T,Cmatrix=W+Diagonal(nrow(W))*1e-05,diagonal =0,model="generic0")+
                 f(interaction,model="generic0",diagonal=0,Cmatrix=Q_st+10^6*(Diagonal(nrow(Q_st))-PC$P)+Diagonal(nrow(Q_st))*1e-05,constr=F,extraconstr=list(A=A_old,e=rep(0,nrow(A_old)))),
                   verbose=T,family="poisson",control.fixed = list(prec.intercept =0.01),
                 num.threads=nthreads,data=combined_data)
  saveRDS(HyPrick,file=fil)
}

fil = paste0("Largedata_Hymik_nt=",nt,".RDS")
if(file.exists(fil))
  HyMiK = readRDS(fil)
if(!file.exists(fil))
{
  HyMiK=inla(obs.x~offset(log(exp.x))+
                 f(year.index,model="generic0",diagonal=0,Cmatrix=GMRF_RW(n=nt,order=2)+Diagonal(nt)*1e-05,constr=T)+
                 f(region.index,constr=T,Cmatrix=W+Diagonal(nrow(W))*1e-05,diagonal =0,model="generic0")+
                 f(interaction,model="z",precision = 1e5, diagonal=0,Z=PC$P,Cmatrix=Q_st+Diagonal(nrow(Q_st))*1e-05,constr=F,
                   extraconstr=list(A=cbind(A_newHymiK,A_newHymiK*0),e=rep(0,nrow(A_newHymiK)))),verbose=T,family="poisson",
           control.fixed = list(prec.intercept =0.01),num.threads=nthreads,data=combined_data)
  saveRDS(HyMiK,file=fil)
}

fil = paste0("Largedata_Bolin_nt=",nt,".RDS")
if(file.exists(fil))
  BOLIN_WALLIN = readRDS(fil)
if(!file.exists(fil))
{
  BOLIN_WALLIN=inla(obs.x~offset(log(exp.x))+f(year.index,model="generic0",diagonal=0,Cmatrix=GMRF_RW(n=nt,order=2)+Diagonal(nt)*1e-05,constr=T)+
                        f(region.index,constr=T,diagonal=0,Cmatrix=W+Diagonal(nrow(W))*1e-05,model="generic0")+
                        f(interaction,Z=TMat$T[,U],precision = 1e5,model="z",diagonal = 0,Cmatrix=TQT[U,U]+1e-05*Diagonal(length(U)),
                          constr=F,extraconstr=list(A=cbind((A_newB),(A_new_BOLIN_sub)*0),e=rep(0,-1+ncol(A_BB)))),verbose=T,
                        family="poisson",control.fixed = list(prec.intercept =0.01),num.threads=nthreads,data=combined_data)
  saveRDS(BOLIN_WALLIN,file=fil)
}

tab =cbind(HyMiK$summary.hyperpar$mean,
           BOLIN_WALLIN$summary.hyperpar$mean,HyPrick$summary.hyperpar$mean)
tab[1,] = 1/tab[1,]
tab = rbind(c(HyMiK$summary.fixed$mean,
              BOLIN_WALLIN$summary.fixed$mean,HyPrick$summary.fixed$mean),
            tab,c(HyMiK$cpu.used[4],BOLIN_WALLIN$cpu.used[4],HyPrick$cpu.used[4]))
colnames(tab) = c("HyMiK","BW","HyPrick")
rownames(tab)=c("Mean","Precision Temporal","Precision Spatial","Precision Interaction","CPU")

library(ggplot2)
library(patchwork)

set.seed(34)
ind = sample(1:(ns*nt),1000,replace=FALSE)
plotDataE=data.table::data.table(HyPrickE=HyPrick$summary.random$interaction$mean,
                                 HyMiKE=HyMiK$summary.random$interaction$mean[1:(nt*ns)],
                                 BWE=BOLIN_WALLIN$summary.random$interaction$mean[1:(nt*ns)])

plotDataSd=data.table::data.table(HyPrickE=HyPrick$summary.random$interaction$sd,
                                  HyMiKE=HyMiK$summary.random$interaction$sd[1:(nt*ns)],
                                  BWE=BOLIN_WALLIN$summary.random$interaction$sd[1:(nt*ns)])

p1 = ggplot(data = plotDataE[ind,]) +
  geom_point(aes(x = HyMiKE, y = BWE, colour = "HyMiK"), size = 1) +
  geom_point(aes(x = HyMiKE,   y = HyPrickE, colour = "HyPrick"),   size = 1) +
  scale_colour_manual(
    name = "Method",
    values = c("HyMiK" = "black", "HyPrick" = "red")
  ) +
  geom_abline(intercept = 0, slope = 1, size = 0.5,color="green") +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  coord_fixed() +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm")) +
  labs(x = "HyMiK", y = "Alternatives")

p2 = ggplot(data = plotDataSd[ind,]) +
  geom_point(aes(x = HyMiKE, y = BWE, colour = "HyMiK"), size = 1) +
  geom_point(aes(x = HyMiKE,   y = HyPrickE, colour = "HyPrick"),   size = 1) +
  scale_colour_manual(
    name = "Method",
    values = c("HyMiK" = "black", "HyPrick" = "red")
  ) +
  geom_abline(intercept = 0, slope = 1, size = 0.5,color="green") +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  coord_fixed() +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm")) +
  labs(x = "HyMik", y = "Alternatives")

(p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave("Large_Interaction_E_sd.pdf")
