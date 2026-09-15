library(Matrix)
library(INLA)
#Can include pardiso if available
#INLA::inla.setOption("pardiso.license","/nr/samba/user/storvik/sys/licences/pardiso.lic")
#INLA::inla.pardiso.check()
library(data.table)
library(ggplot2)
library(xtable)
library(Matrix)
dir = "functions/"
source(paste0(dir,"GMRF_RW.R"))
source(paste0(dir,"inla_KnorrHeld4.R"))

df = readRDS("data/SpatioTemporalData.RDS")
ns = max(df$main_spatial)
nt = max(df$main_temporal)

graph=system.file("demodata/germany.graph", package="INLA")

Q_ICAR=INLA::inla.graph2matrix(graph)
diag(Q_ICAR)=0
Q_ICAR=-Q_ICAR
diag(Q_ICAR)=-rowSums(Q_ICAR)

Q_RW2=GMRF_RW(n=nt,order=2)

sc = FALSE

fil = "res_new/simP.HyPrick1.RDS"
if(file.exists(fil))
  HyPrick1 = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick1 = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson",method="hyprick",kappa=1e+04,scale=sc)
  saveRDS(HyPrick1,file=fil)
} 

fil = "res_new/simP.HyPrick2.RDS"
if(file.exists(fil))
  HyPrick2 = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick2 = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson",method="hyprick",kappa=1e+05,scale=sc)
  saveRDS(HyPrick2,file=fil)
} 

fil = "res_new/simP.HyPrick.RDS"
if(file.exists(fil))
  HyPrick3 = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick3 = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson",method="hyprick",kappa=1e+06,scale=sc)
  saveRDS(HyPrick3,file=fil)
} 

fil = "res_new/simP.HyPrick4.RDS"
if(file.exists(fil))
  HyPrick4 = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick4 = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson",method="hyprick",kappa=1e+07,scale=sc)
  saveRDS(HyPrick4,file=fil)
} 

if(0)
{
library(xtable)
library(patchwork)
tab =cbind(HyPrick1$summary.hyperpar$mean,HyPrick2$summary.hyperpar$mean,
           HyPrick3$summary.hyperpar$mean,HyPrick4$summary.hyperpar$mean)
tab = rbind(c(HyPrick1$summary.fixed$mean,HyPrick2$summary.fixed$mean,
              HyPrick3$summary.fixed$mean,HyPrick4$summary.fixed$mean),
            tab,c(HyPrick1$cpu.used[4],HyPrick2$cpu.used[4],HyPrick3$cpu.used[4],HyPrick4$cpu.used[4]))
colnames(tab) = c("Kappa=1e04","Kappa=1e05","Kappa=1e06","Kappa=1e07")
rownames(tab)=c("Mean","Precision Temporal","Precision Spatial","Precision Interaction","CPU")
xtable(tab,digits=3)
print(xtable(tab,digits=3),type="latex",file="tab_simP_sens.tex")


plotDataE=data.table::data.table(HyPrickE1=HyPrick1$summary.random$delta$mean,
                                 HyPrickE2=HyPrick2$summary.random$delta$mean,
                                 HyPrickE3=HyPrick3$summary.random$delta$mean,
                                 HyPrickE4=HyPrick4$summary.random$delta$mean)

plotDataSd=data.table::data.table(HyPrickE1=HyPrick1$summary.random$delta$sd,
                                  HyPrickE2=HyPrick2$summary.random$delta$sd,
                                  HyPrickE3=HyPrick3$summary.random$delta$sd,
                                  HyPrickE4=HyPrick4$summary.random$delta$sd)

p1 = ggplot(data = plotDataE) +
  geom_point(aes(y = HyPrickE1, x = HyPrickE3, colour = "HyPrick1"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyPrickE2, x = HyPrickE3, colour = "HyPrick2"),   size = 1,alpha=0.2)  +
  geom_point(aes(y = HyPrickE4, x = HyPrickE3, colour = "HyPrick4"),   size = 1,alpha=0.2)  +
  scale_colour_manual(
    name = "Method",
    values = c("HyPrick1" = "black", "HyPrick2" = "blue", "HyPrick4" = "red"),
    breaks = c("HyPrick1", "HyPrick2", "HyPrick4"),
    labels = c("HyPrick1" = "Kappa=1e04",
               "HyPrick2" = "Kappa=1e05",
               "HyPrick4" = "Kappa=1e07")
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
  labs(x = "Kappa=1e06", y = "Alternatives")

p2 = ggplot(data = plotDataSd) +
  geom_point(aes(y = HyPrickE1, x = HyPrickE2, colour = "HyPrick1"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyPrickE2, x = HyPrickE2, colour = "HyPrick2"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyPrickE4, x = HyPrickE2, colour = "HyPrick4"),   size = 1,alpha=0.2) +
  scale_colour_manual(
    name = "Method",
    values = c("HyPrick1" = "black", "HyPrick2" = "blue", "HyPrick4" = "red"),
    breaks = c("HyPrick1", "HyPrick2", "HyPrick4"),
    labels = c("HyPrick1" = "Kappa=1e04",
               "HyPrick2" = "Kappa=1e05",
               "HyPrick4" = "Kappa=1e07")
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
  labs(x = "Kappa=1e06", y = "Alternatives")

(p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave("res_new/SimP_Interaction_E_sd_sens.png",device="png",height=5,width=5)
}