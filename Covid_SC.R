##INLA analysis of covid data
library(INLA)
#Can turn on pardiso if available
#INLA::inla.setOption("pardiso.license","~/sys/licences/pardiso.lic")
#INLA::inla.pardiso.check()
library(Matrix)
library(data.table)
library(dlnm)
library(ggplot2)
library(patchwork)
library(xtable)
rm(list=ls())

dir = "functions/"
source(paste0(dir,"GMRF_RW.R"))
source(paste0(dir,"inla_KnorrHeld4.R"))


#Read data and make indices for temporal and interaction terms
load("data/coviddata.rda")
df = coviddata$data
df$weekday = as.factor(weekdays(df$date))
df = df[df$date>as.Date("2020-10-01"),]
df$T1 = as.numeric(df$date-df$date[1]+1)
nt = max(df$T1)
df$county = as.numeric(as.factor(df$location_code))
ns = max(df$county)

##Reduce number of time points
#nt = 500
df = df[df$T1<(nt+1),]
df$E = log(df$pop)
df$S1T1 = (df$T1-1)*ns+df$county  

#Make precision matrices
Q_ICAR = -coviddata$adj
for(i in 1:ns)
  Q_ICAR[i,i] = -sum(Q_ICAR[i,-i])
Q_RW2=GMRF_RW(n=nt,order=2)

extracov = 'offset(E)'
indd = match(c("county","T1","S1T1","cases"),names(df))
sc = FALSE

fil = "res_new/covid.Standard.RDS"
if(file.exists(fil))
  StandardINLA = readRDS(fil)
if(!file.exists(fil))
{
  StandardINLA = inla_KnorrHeld4(df,Qtemp=Q_RW2,Qspat=Q_ICAR,indd=indd,extracov=extracov,family="poisson",scale=sc)
  #control.predictor=list(compute=TRUE
  saveRDS(StandardINLA,file=fil)
}

fil = "res_new/covid.HyMiK.RDS" 
if(file.exists(fil))
  HyMiK = readRDS(fil)
if(!file.exists(fil))
{
  HyMiK = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,indd=indd,extracov=extracov,family="poisson",method="hymik",scale=sc)
  saveRDS(HyMiK,file=fil)
}

fil = "res_new/covid.HyBW.RDS"
if(file.exists(fil))
  HyBW = readRDS(fil)
if(!file.exists(fil))
{
  #rcpp called.
  Rcpp::sourceCpp(paste0(dir,"cp.cpp"))
  HyBW = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,indd=indd,extracov=extracov,family="poisson",method="hybw",scale=sc)
  saveRDS(HyBW,file=fil)
}

fil = "res_new/covid.HyPrick.RDS"
if(file.exists(fil))
  HyPrick = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,indd=indd,extracov=extracov,family="poisson",method="hyprick",scale=sc)
  saveRDS(HyPrick,file=fil)
} 

if(0)
{ ## Not working!!
fil = "res_new/covid.BW.RDS"
if(file.exists(fil))
  BW = readRDS(fil)
if(!file.exists(fil))
{
  BW = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,indd=indd,extracov=extracov,family="poisson",method="bw")
  saveRDS(BW,file=fil)
}
}

if(0)
{
tab =cbind(StandardINLA$summary.hyperpar$mean,HyMiK$summary.hyperpar$mean,
           HyBW$summary.hyperpar$mean,HyPrick$summary.hyperpar$mean)
tab = rbind(c(StandardINLA$summary.fixed$mean,HyMiK$summary.fixed$mean,
              HyBW$summary.fixed$mean,HyPrick$summary.fixed$mean),
            tab,c(StandardINLA$cpu.used[4],HyMiK$cpu.used[4],HyBW$cpu.used[4],HyPrick$cpu.used[4]))
colnames(tab) = c("Standard","HyMiK","HyBW","HyPrick")
rownames(tab)=c("Mean","Precision Temporal","Precision Spatial","Precision Interaction","CPU")
print(tab)
print(xtable(tab,digits=3),type="latex",file="tab_covid.tex")

plotDataE=data.table::data.table(StandardINLAE=StandardINLA$summary.random$delta$mean,
                                 HyPrickE=HyPrick$summary.random$delta$mean,
                                 HyMiKE=HyMiK$summary.random$delta$mean[1:(nt*ns)],
                                 HyBW=HyBW$summary.random$delta$mean[1:(nt*ns)])

plotDataSd=data.table::data.table(StandardINLAE=StandardINLA$summary.random$delta$sd,
                                  HyPrickE=HyPrick$summary.random$delta$sd,
                                  HyMiKE=HyMiK$summary.random$delta$sd[1:(nt*ns)],
                                  HyBW=HyBW$summary.random$delta$sd[1:(nt*ns)])

p1 = ggplot(data = plotDataE) +
  geom_point(aes(y = HyPrickE, x = StandardINLAE, colour = "HyPrick"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyMiKE,   x = StandardINLAE, colour = "HyMiK"),   size = 1,alpha=0.2) +
  geom_point(aes(y = HyBW,      x = StandardINLAE, colour = "HyBW"),      size = 1,alpha=0.2) +
  scale_colour_manual(
    name = "Method",
    values = c("HyPrick" = "black", "HyMiK" = "red", "HyBW" = "blue")
  ) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.05,color="black") +
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
  geom_point(aes(y = HyPrickE, x = StandardINLAE, colour = "HyPrick"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyMiKE,   x = StandardINLAE, colour = "HyMiK"),   size = 1,alpha=0.2) +
  geom_point(aes(y = HyBW,      x = StandardINLAE, colour = "HyBW"),      size = 1,alpha=0.2) +
  scale_colour_manual(
    name = "Method",
    values = c("HyPrick" = "black", "HyMiK" = "red", "HyBW" = "blue")
  ) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.1,color="black") +
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

ggsave("res_new/covid_Interaction_E_sd.png",device="png",height=5,width=5)
}
