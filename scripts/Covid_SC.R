##INLA analysis of covid data
library(INLA)
#Can turn on pardiso if available
#INLA::inla.setOption("pardiso.license","~/sys/licences/pardiso.lic")
#INLA::inla.pardiso.check()
library(Matrix)
library(data.table)
library(dlnm)
library(ggplot2)
library(xtable)

rm(list=ls())
dir =""
dirFun = paste0(dir,"/functions/")
dirdata=paste0(dir,"/data/")
#Update:
source(paste0(dirFun,"SpaceTimeProjConstr.R"))
source(paste0(dirFun,"GMRF_RW.R"))
#rcpp called.
Rcpp::sourceCpp(paste0(dirFun,"cp.cpp"))
nthreads = "3:2"

#Read data and make indices for temporal and interaction terms
load(file=paste0(dirdata,"/coviddata.rda"))

df = coviddata$data
df$weekday = as.factor(weekdays(df$date))
df = df[df$date>as.Date("2020-10-01"),]
df$temp = as.numeric(df$date-df$date[1]+1)

#Important: how large nt?
nt = max(df$temp)

df$spat = as.numeric(as.factor(df$location_code))
ns = max(df$spat)


df = df[df$temp<(nt+1),]

df$interaction = (df$temp-1)*ns+df$spat  


#Make precision matrices
Q_ICAR = -coviddata$adj
for(i in 1:ns)
  Q_ICAR[i,i] = -sum(Q_ICAR[i,-i])
Q_RW2=GMRF_RW(n=nt,order=2)

eps = 1e-05
kap = 1e06


SCALE = FALSE
if(SCALE)
{
  #Scale model
  A_t=matrix(rep(1,nt),nrow=1)
  Q_RW2 = inla.scale.model(Q_RW2+diag(nt)*eps,list(A=A_t,e=0))
  Q_ICAR=inla.scale.model(Q_ICAR+diag(ns)*eps,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
}


Q_st=kronecker(Q_RW2,Q_ICAR)

PC = SpaceTimeProjConstr(ns = ns,nt = nt,dim="time",type="SC")
A21=((PC$A2[,1:(ns*nt)])%*%PC$P)


A_Bolin=t(kronecker(Diagonal(nt),matrix(rep(1,ns),nrow=1)))


#Alternative 1:
A_Bolin=as(A_Bolin,"dgCMatrix")

#Alternative 2:
A_Bolin=as(A_Bolin,"CsparseMatrix")

TMat=c_basis2(t(A_Bolin))
eps=1e-05

QT=(Q_st+Diagonal(ns*nt)*eps)%*%(TMat$T)
QT=drop0(QT,tol=0)
TQT=(t(TMat$T)%*%QT)
TQT=drop0(TQT,tol=0)

C=1:nrow(t(A_Bolin))
U=(1+nrow(t(A_Bolin))):(ncol(t(A_Bolin)))

A_BB=PC$A2[,seq(1,ns*nt)]
A_newHymiK=t(PC$P%*%t(A_BB))


A_newB=t(PC$P%*%t(PC$A2[,seq(1,ns*nt)]))
A_new_BOLIN=A_newB%*%t(TMat$T)
A_new_BOLIN_sub=A_new_BOLIN[,U]

A_old=t(A_BB)
A_BB_alternative=t((A_BB)%*%TMat$T[,U])

df$E = log(df$pop)

formula.standard = baseformula =cases~offset(E)+#weekday+
  f(temp,model="generic0",Cmatrix=Q_RW2+Diagonal(nt)*eps,constr=T,diagonal = 0)+
  f(spat,model="generic0",Cmatrix=Q_ICAR+Diagonal(ns)*eps,constr=T,diagonal = 0) +
  f(interaction,model="generic0",Cmatrix=Q_st+Diagonal(ns*nt)*eps,constr=F,
    extraconstr=list(A=PC$A,e=rep(0,nrow(PC$A))),diagonal = 0)

formula.HyPrick =cases~offset(E)+#weekday+
  f(temp,model="generic0",Cmatrix=Q_RW2+Diagonal(nt)*eps,constr=T,diagonal = 0)+
  f(spat,model="generic0",Cmatrix=Q_ICAR+Diagonal(ns)*eps,constr=T,diagonal = 0) +
  f(interaction,model="generic0",
    Cmatrix=Q_st+10^6*(Diagonal(nrow(Q_st))-PC$P)+Diagonal(nrow(Q_st))*1e-05,constr=F,extraconstr=list(A=t(A_old),e=rep(0,ncol(A_old))))

formula.HyMiK = baseformula =cases~offset(E)+#weekday+
  f(temp,model="generic0",Cmatrix=Q_RW2+Diagonal(nt)*eps,constr=T,diagonal = 0)+
  f(spat,model="generic0",Cmatrix=Q_ICAR+Diagonal(ns)*eps,constr=T,diagonal = 0) +
  f(interaction,model="z",Z=PC$P,Cmatrix=Q_st+Diagonal(ns*nt)*eps,constr=F,precision = 1e5,
    extraconstr=list(A=cbind(A_newHymiK,A_newHymiK*0),e=rep(0,nrow(A_newHymiK))),diagonal = 0)

fil = "covid.Standard.RDS"
if(file.exists(fil))
  StandardINLA = readRDS(fil)
if(!file.exists(fil))
{
  StandardINLA=inla(formula=formula.standard, family = "poisson",data =df,num.threads =nthreads,verbose=T,control.fixed=list(prec.intercept =0.01))
  saveRDS(StandardINLA,file=fil)
}



fil = "covid.HyMiK1.RDS"
if(file.exists(fil))
  HyMiK = readRDS(fil)
if(!file.exists(fil))
{
 HyMiK=inla(formula=formula.HyMiK, family = "poisson",data =df,num.threads =nthreads,verbose=T,control.fixed=list(prec.intercept =0.01))
 saveRDS(HyMiK,file=fil)
}

fil = "covid.HyPrick.RDS"
if(file.exists(fil))
  HyPrick = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick=inla(formula=formula.HyPrick, family = "poisson",data =df,num.threads =nthreads,verbose=T,control.fixed=list(prec.intercept =0.01))
  saveRDS(HyPrick,file=fil)
}

formula.BolinWallin=cases~offset(E)+
  f(temp,model="generic0",Cmatrix=Q_RW2+Diagonal(nt)*eps,constr=T,diagonal = 0)+
  f(spat,model="generic0",Cmatrix=Q_ICAR+Diagonal(ns)*eps,constr = T,diagonal = 0) +
  f(interaction,Z=TMat$T[,U],precision = 1e5,model="z",diagonal = 0,Cmatrix=TQT[U,U]+1e-05*Diagonal(length(U)),
    constr=F,extraconstr=list(A=cbind((A_newB),(A_new_BOLIN_sub)*0),e=rep(0,nrow(A_newB))))

fil = "covid.BolinWallin.RDS"
if(file.exists(fil))
  BolinWallin = readRDS(fil)
if(!file.exists(fil))
{
  BolinWallin = inla(formula.BolinWallin,family = "poisson",data =df,num.threads =nthreads,verbose=T,control.fixed=list(prec.intercept =0.01))
  saveRDS(BolinWallin,file=fil)
}

print(HyPrick$cpu.used)
print(HyMiK$cpu.used)
print(BolinWallin$cpu.used)

#Summary
tab =cbind(StandardINLA$summary.hyperpar$mean,HyMiK$summary.hyperpar$mean,
           BolinWallin$summary.hyperpar$mean,HyPrick$summary.hyperpar$mean)
tab[1,] = 1/tab[1,]
tab = rbind(c(StandardINLA$summary.fixed$mean,HyMiK$summary.fixed$mean,
              BolinWallin$summary.fixed$mean,HyPrick$summary.fixed$mean),
            tab,c(StandardINLA$cpu.used[4],HyMiK$cpu.used[4],BolinWallin$cpu.used[4],HyPrick$cpu.used[4]))
colnames(tab) = c("Standard","HyMiK","BW","HyPrick")
rownames(tab)=c("Mean","Precision Temporal","Precision Spatial","Precision Interaction","CPU")




#plotting:

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
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed() +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  labs(x = "Standard", y = "Alternatives")

p2 = ggplot(data = plotDataSd) +
  geom_point(aes(y = HyPrickE, x = StandardINLAE, colour = "HyPrick"), size = 1) +
  geom_point(aes(y = HyMiKE,   x = StandardINLAE, colour = "HyMiK"),   size = 1) +
  geom_point(aes(y = BWE,      x = StandardINLAE, colour = "BW"),      size = 1) +
  geom_abline(intercept = 0, slope = 1, size = 0.5,color="green") +
  scale_colour_manual(
    name = "Method",
    values = c("HyPrick" = "black", "HyMiK" = "red", "BW" = "blue")
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed() +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  labs(x = "Standard", y = "Alternatives")

(p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "right")


ggsave("Covid_Interaction_E_sd.pdf")