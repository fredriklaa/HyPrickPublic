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

#Make data
data("Carto_SpainMUN")
data(Data_LungCancer)
carto.nb <- poly2nb(Carto_SpainMUN)
Carto_SpainMUN_connected=Carto_SpainMUN[-2454,]
#Carto_SpainMUN_connected=Carto_SpainMUN_connected[1:50,]
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

## Not working
#fil = "res_new/Largedata.Standard.RDS"
#if(file.exists(fil))
#  StandardINLA = readRDS(fil)
#if(!file.exists(fil))
#{
#  StandardINLA = inla_KnorrHeld4(combined_data,Q_RW2,Qspat,indd=indd,extracov=extracov,family="poisson",scale=sc)
#  saveRDS(StandardINLA,file=fil)
#}

fil = "res_new/Largedata.HyMiK.RDS" 
if(file.exists(fil))
  HyMiK = readRDS(fil)
if(!file.exists(fil))
{
  HyMiK = inla_KnorrHeld4(combined_data,Q_RW2,Qspat,indd=indd,extracov=extracov,family="poisson",method="hymik",scale=sc)
  saveRDS(HyMiK,file=fil)
}

fil = "res_new/Largedata.bolinWallin.RDS"
if(file.exists(fil))
  BolinWallin = readRDS(fil)
if(!file.exists(fil))
{
  #rcpp called.
  Rcpp::sourceCpp(paste0(dir,"cp.cpp"))
  BolinWallin = inla_KnorrHeld4(combined_data,Q_RW2,Qspat,indd=indd,extracov=extracov,family="poisson",method="hybw",scale=sc)
  saveRDS(BolinWallin,file=fil)
}

fil = "res_new/Largedata.HyPrick.RDS"
if(file.exists(fil))
  HyPrick = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick = inla_KnorrHeld4(combined_data,Q_RW2,Qspat,indd=indd,extracov=extracov,family="poisson",method="hyprick",scale=sc)
  saveRDS(HyPrick,file=fil)
} 

#{ ## Not working
#  fil = "res_new/simP.BW.RDS"
#  if(file.exists(fil))
#    BW = readRDS(fil)
#  if(!file.exists(fil))
#  {
#    BW = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson",method="bw",scale=sc)
#    saveRDS(BW,file=fil)
#  }


tab =cbind(HyMiK$summary.hyperpar$mean,
           BolinWallin$summary.hyperpar$mean,HyPrick$summary.hyperpar$mean)
tab = rbind(c(HyMiK$summary.fixed$mean,
              BolinWallin$summary.fixed$mean,HyPrick$summary.fixed$mean),
            tab,c(HyMiK$cpu.used[4],BolinWallin$cpu.used[4],HyPrick$cpu.used[4]))
colnames(tab) = c("HyMiK","HyBW","HyPrick")
rownames(tab)=c("Mean","Precision Temporal","Precision Spatial","Precision Interaction","CPU")
print(tab)

library(ggplot2)
library(patchwork)

plotDataE=data.table::data.table(HyPrickE=HyPrick$summary.random$delta$mean,
                                 HyMiKE=HyMiK$summary.random$delta$mean[1:(nt*ns)],
                                 HyBW=BolinWallin$summary.random$delta$mean[1:(nt*ns)])

plotDataE = plotDataE[seq(1, .N, by = 10)]
#plotDataE = plotDataE[order(plotDataE$HyPrickE),]
#ind = c(1:197650)*10
#plotDataE = plotDataE[ind,]

plotDataSd=data.table::data.table(HyPrickE=HyPrick$summary.random$delta$sd,
                                  HyMiKE=HyMiK$summary.random$delta$sd[1:(nt*ns)],
                                  HyBW=BolinWallin$summary.random$delta$sd[1:(nt*ns)])

#Picking every 10th obsevation to make graphical file smaller
plotDataSd = plotDataSd[seq(1, .N, by = 10)]


p1 = ggplot(data = plotDataE) +
  geom_point(aes(y = HyMiKE,   x = HyPrickE, colour = "HyMiK"),   size = 1,alpha=0.2) +
  geom_point(aes(y = HyBW,      x = HyPrickE, colour = "HyBW"),      size = 1,alpha=0.2) +
  scale_colour_manual(
    name = "Method",
    values = c("HyMiK" = "red", "HyBW" = "blue")
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
  labs(x = "HyPricK", y = "Alternatives")

p2 = ggplot(data = plotDataSd) +
  geom_point(aes(y = HyMiKE,   x = HyPrickE, colour = "HyMiK"),   size = 1,alpha=0.2) +
  geom_point(aes(y = HyBW,      x = HyPrickE, colour = "HyBW"),      size = 1,alpha=0.2) +
  scale_colour_manual(
    name = "Method",
    values = c("HyMiK" = "red", "HyBW" = "blue")
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
  labs(x = "HyPricK", y = "Alternatives")

(p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave("res_new/LargeData_Interaction_E_sd.pdf",height=5,width=5)




