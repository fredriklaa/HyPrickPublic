library(Matrix)
library(INLA)
#Can include pardiso if available
#INLA::inla.setOption("pardiso.license","/nr/samba/user/storvik/sys/licences/pardiso.lic")
#INLA::inla.pardiso.check()
library(data.table)
library(ggplot2)
library(patchwork)
library(xtable)
library(Matrix)
dir = "functions/"
source(paste0(dir,"GMRF_RW.R"))
source(paste0(dir,"inla_KnorrHeld4.R"))


df = readRDS("data/SpatioTemporalDataSmall.RDS")
ns = max(df$main_spatial)
nt = max(df$main_temporal)

foo = readRDS("data/SpatioTemporalDataSmall_Q.RDS")
Q_ICAR = foo$Q_ICAR
Q_RW2 = foo$Q_RW2

fil = "res_new/SimPSmall.Standard.RDS"
if(file.exists(fil))
  StandardINLA = readRDS(fil)
if(!file.exists(fil))
{
  StandardINLA = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson")
  saveRDS(StandardINLA,file=fil)
}

fil = "res_new/SimPSmall.HyMiK.RDS" 
if(file.exists(fil))
  HyMiK = readRDS(fil)
if(!file.exists(fil))
{
  HyMiK = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson",method="hymik")
  saveRDS(HyMiK,file=fil)
}

fil = "res_new/SimPSmall.HyBW.RDS"
if(file.exists(fil))
  HyBW = readRDS(fil)
if(!file.exists(fil))
{
  #rcpp called.
  Rcpp::sourceCpp(paste0(dir,"cp.cpp"))
  HyBW = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson",method="hybw")
  saveRDS(HyBW,file=fil)
}

fil = "res_new/SimPSmall.HyPrick.RDS"
if(file.exists(fil))
  HyPrick = readRDS(fil)
if(!file.exists(fil))
{
  HyPrick = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson",method="hyprick")
  saveRDS(HyPrick,file=fil)
} 

fil = "res_new/SimPSmall.BW.RDS"
if(file.exists(fil))
  BW = readRDS(fil)
if(!file.exists(fil))
{
  BW = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson",method="bw")
  saveRDS(BW,file=fil)
}

fil = "res_new/SimPSmall.BW2.RDS"
if(file.exists(fil))
  BW2 = readRDS(fil)
if(!file.exists(fil))
{
  BW2 = inla_KnorrHeld4(df,Q_RW2,Q_ICAR,family="poisson",method="bw_A_formulation")
  saveRDS(BW2,file=fil)
}

if(0)
{
tab =cbind(StandardINLA$summary.hyperpar$mean,BW$summary.hyperpar$mean,BW2$summary.hyperpar$mean,HyMiK$summary.hyperpar$mean,
           HyBW$summary.hyperpar$mean,HyPrick$summary.hyperpar$mean)
tab = rbind(c(StandardINLA$summary.fixed$mean,BW$summary.fixed$mean,BW2$summary.fixed$mean,HyMiK$summary.fixed$mean,
              HyBW$summary.fixed$mean,HyPrick$summary.fixed$mean),
            tab,c(StandardINLA$cpu.used[4],BW$cpu.used[4],BW2$cpu.used[4],HyMiK$cpu.used[4],HyBW$cpu.used[4],HyPrick$cpu.used[4]))
colnames(tab) = c("Standard","BW","BW2","HyMiK","HyBW","HyPrick")
rownames(tab)=c("Mean","Precision Temporal","Precision Spatial","Precision Interaction","CPU")
print(tab)
print(xtable(tab,digits=3),type="latex",file="tab_simP_small.tex")


plotDataE=data.table::data.table(StandardINLAE=StandardINLA$summary.random$delta$mean,
                                 BWE=BW$summary.random$delta$mean[1:(nt*ns)],
                                 BWE2=BW2$summary.random$delta$mean[1:(nt*ns)],
                                 HyPrickE=HyPrick$summary.random$delta$mean,
                                 HyMiKE=HyMiK$summary.random$delta$mean[1:(nt*ns)],
                                 HyBW=HyBW$summary.random$delta$mean[1:(nt*ns)])

plotDataSd=data.table::data.table(StandardINLAE=StandardINLA$summary.random$delta$sd,
                                  BWE=BW$summary.random$delta$sd[1:(nt*ns)],
                                  BWE2=BW$summary.random$delta$sd[1:(nt*ns)],
                                  HyPrickE=HyPrick$summary.random$delta$sd,
                                  HyMiKE=HyMiK$summary.random$delta$sd[1:(nt*ns)],
                                  HyBW=HyBW$summary.random$delta$sd[1:(nt*ns)])

p1 = ggplot(data = plotDataE) +
  geom_point(aes(y = BWE, x = StandardINLAE, colour = "BW"), size = 1,alpha=0.2) +
  geom_point(aes(y = BWE2, x = StandardINLAE, colour = "BW2"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyPrickE, x = StandardINLAE, colour = "HyPrick"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyMiKE,   x = StandardINLAE, colour = "HyMiK"),   size = 1,alpha=0.2) +
  geom_point(aes(y = HyBW,      x = StandardINLAE, colour = "HyBW"),      size = 1,alpha=0.2) +
  scale_colour_manual(
    name = "Method",
    values = c("BW" = "pink", "BW2" = "orange", "HyPrick" = "black", "HyMiK" = "red", "HyBW" = "blue")
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
  geom_point(aes(y = BWE, x = StandardINLAE, colour = "BW"), size = 1,alpha=0.2) +
  geom_point(aes(y = BWE2, x = StandardINLAE, colour = "BW2"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyPrickE, x = StandardINLAE, colour = "HyPrick"), size = 1,alpha=0.2) +
  geom_point(aes(y = HyMiKE,   x = StandardINLAE, colour = "HyMiK"),   size = 1,alpha=0.2) +
  geom_point(aes(y = HyBW,      x = StandardINLAE, colour = "HyBW"),      size = 1,alpha=0.2) +
  scale_colour_manual(
    name = "Method",
    values = c("BW" = "pink", "BW2" = "orange", "HyPrick" = "black", "HyMiK" = "red", "HyBW" = "blue")
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

ggsave("res_new/SimPSmall_Interaction_E_sd.png",device="png",height=5,width=5)
}