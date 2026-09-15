#library(INLAconstraints)
library(Matrix)
library(sf)
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
library(bigDM)
library(spdep)
rm(list=ls())

dir = "functions/"
source(paste0(dir,"SpaceTimeProjConstr.R"))
source(paste0(dir,"GMRF_RW.R"))
source(paste0(dir,"inla_KnorrHeld4.R"))
#rcpp called.
Rcpp::sourceCpp(paste0(dir,"cp.cpp"))

#Make data
data("Carto_SpainMUN")
data(Data_LungCancer)
carto.nb <- poly2nb(Carto_SpainMUN)
Carto_SpainMUN_connected=Carto_SpainMUN[-2454,]
carto.nb2 <- poly2nb(Carto_SpainMUN_connected)
Qspat=nb2mat(carto.nb2,style = "B")
Qspat=-Qspat
diag(Qspat)=-rowSums(Qspat)


#can be between 4 and 25
nt=25
library(dplyr)
combined_data=inner_join(Data_LungCancer,Carto_SpainMUN_connected,by=c("ID"="ID"))
combined_data$year.index=combined_data$year-min(combined_data$year)+1
combined_data$region.index=as.numeric(as.factor(combined_data$ID))
combined_data=combined_data[combined_data$year.index<=nt,]
combined_data$interaction=(combined_data$year.index-1)*length(carto.nb2)+combined_data$region.index

#nt=max(combined_data$year.index)
ns =nrow(Qspat)
Q_RW2=GMRF_RW(n=nt,order=2)

extracov = 'offset(log(exp.x))'

indd = match(c("region.index","year.index","interaction","obs.x"),names(combined_data))

sc = FALSE

fil = "res_new/Largedata.HyPrick1.RDS"
if(file.exists(fil))
  HyPrick1 = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick1 = inla_KnorrHeld4(combined_data,Q_RW2,Qspat,indd=indd,extracov=extracov,family="poisson",method="hyprick",kappa=1e+04,scale=sc)
  saveRDS(HyPrick1,file=fil)
} 

fil = "res_new/Largedata.HyPrick2.RDS"
if(file.exists(fil))
  HyPrick2 = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick2 = inla_KnorrHeld4(combined_data,Q_RW2,Qspat,indd=indd,extracov=extracov,family="poisson",method="hyprick",kappa=1e+05,scale=sc)
  saveRDS(HyPrick2,file=fil)
} 

fil = "res_new/Largedata.HyPrick.RDS"
if(file.exists(fil))
  HyPrick3 = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick3 = inla_KnorrHeld4(combined_data,Q_RW2,Qspat,indd=indd,extracov=extracov,family="poisson",method="hyprick",kappa=1e+06,scale=sc)
  saveRDS(HyPrick3,file=fil)
} 

fil = "res_new/Largedata.HyPrick4.RDS"
if(file.exists(fil))
  HyPrick4 = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick4 = inla_KnorrHeld4(combined_data,Q_RW2,Qspat,indd=indd,extracov=extracov,family="poisson",method="hyprick",kappa=1e+07,scale=sc)
  saveRDS(HyPrick4,file=fil)
} 

tab =cbind(HyPrick1$summary.hyperpar$mean,HyPrick2$summary.hyperpar$mean,
           HyPrick3$summary.hyperpar$mean,HyPrick4$summary.hyperpar$mean)
tab = rbind(c(HyPrick1$summary.fixed$mean,HyPrick2$summary.fixed$mean,
              HyPrick3$summary.fixed$mean,HyPrick4$summary.fixed$mean),
            tab,c(HyPrick1$cpu.used[4],HyPrick2$cpu.used[4],HyPrick3$cpu.used[4]),HyPrick4$cpu.used[4])
colnames(tab) = c("HyPrick1","HyPrick2","HyPrick3","HyPrick4")
rownames(tab)=c("Mean","Precision Temporal","Precision Spatial","Precision Interaction","CPU")

plotDataE=data.table::data.table(HyPrickE1=HyPrick1$summary.random$delta$mean,
                                 HyPrickE2=HyPrick2$summary.random$delta$mean,
                                 HyPrickE3=HyPrick3$summary.random$delta$mean,
                                 HyPrickE4=HyPrick4$summary.random$delta$mean)

plotDataSd=data.table::data.table(HyPrickE1=HyPrick1$summary.random$delta$sd,
                                  HyPrickE2=HyPrick2$summary.random$delta$sd,
                                  HyPrickE3=HyPrick3$summary.random$delta$sd)

p1 = ggplot(data = plotDataE) +
  geom_point(aes(y = HyPrickE1, x = HyPrickE3, colour = "HyPrick2"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyPrickE2,   x = HyPrickE3, colour = "HyPrick3"),   size = 1,alpha=0.2)  +
  geom_point(aes(y = HyPrickE4,   x = HyPrickE3, colour = "HyPrick4"),   size = 1,alpha=0.2)  +
  scale_colour_manual(
    name = "Method",
    values = c("HyPrick1" = "black", "HyPrick2" = "blue", "HyPrick4" = "red"),
    breaks = c("HyPrick1", "HyPrick2", "HyPrick4"),
    labels = c("HyPrick1" = "Kappa=1e04",
               "HyPrick2" = "Kappa=1e05",
               "HyPrick4" = "Kappa=1e07")
  ) +
  geom_abline(intercept = 0, slope = 1, size = 0.05,color="black") +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  coord_fixed() +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm")) +
  labs(x = "Kappa=1e06", y = "Alternatives")

p2 = ggplot(data = plotDataSd) +
  geom_point(aes(y = HyPrickE2, x = HyPrickE, colour = "HyPrick2"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyPrickE3,   x = HyPrickE, colour = "HyPrick3"),   size = 1,alpha=0.2) +
  scale_colour_manual(
    name = "Method",
    values = c("HyPrick1" = "black", "HyPrick2" = "blue", "HyPrick4" = "red"),
    breaks = c("HyPrick1", "HyPrick2", "HyPrick4"),
    labels = c("HyPrick1" = "Kappa=1e04",
               "HyPrick2" = "Kappa=1e05",
               "HyPrick4" = "Kappa=1e07")
  ) +
  geom_abline(intercept = 0, slope = 1, size = 0.1,color="black") +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)) +
  coord_fixed() +
  theme(plot.margin = grid::unit(c(0,0,0,0), "mm")) +
  labs(x = "Kappa=1e06", y = "Alternatives")

library(patchwork)
(p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave("res_new/LargeData_Interaction_E_sd_sens.pdf")

