### figure3.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2025 (11:43) 
## Version: 
## Last-Updated: mar 23 2025 (11:13) 
##           By: Brice Ozenne
##     Update #: 49
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * libraries
library(LMMstar)
library(ggplot2)
library(ggpubr)
source("FCT.R")
source("cor.testIID.R")

## * settings

## ** parametrisation
## param.lmm <- readRDS("results/strategy3-lmm-param.rds")
mu.asl <- c("amygdala" = 41.75461171, "antcin" = 63.97868936, "cau" = 51.03709735, "DLPFC" = 61.81180998, "entorhinc" = 37.93286726, "hippocampus" = 48.20552864, "ins" = 61.01470806, "medinffg" = 64.47818102, "occ" = 49.06232099, "orbfrc" = 54.13080092, "parc" = 57.03041177, "postcin" = 56.59856003, "put" = 53.58549712, "senmotc" = 57.3037015, "supfg" = 56.45369853, "suptempg" = 59.07503085, "th" = 47.54133327, "VLPFC" = 64.91376706) ## param.lmm$mu.asl
mu.pet <- c("amygdala" = 0.83208117, "antcin" = 1.56758472, "cau" = 0.24240839, "DLPFC" = 1.20843445, "entorhinc" = 0.54770079, "hippocampus" = 0.55927805, "ins" = 1.57503343, "medinffg" = 1.27606528, "occ" = 1.15537147, "orbfrc" = 1.44086414, "parc" = 1.1357671, "postcin" = 1.36909301, "put" = 0.35256114, "senmotc" = 0.84761643, "supfg" = 1.13908271, "suptempg" = 1.47829904, "th" = 0.42062905, "VLPFC" = 1.3450188) ## param.lmm$mu.pet
sigma.asl <- 9.270355  ## param.lmm$sigma.asl
sigma.pet <- 0.2546338  ## param.lmm$sigma.pet

rho.asl <- 0.4815182 ## param.lmm$rho.asl
rho.pet <- 0.6382445 ## param.lmm$rho.pet
rho.petasl <- 0.1381907 ## param.lmm$rho.petasl
rhoLag.petasl <- 0.1049642  ## param.lmm$rhoLag.petasl

rho.marginal <- 0.25 ## substantially differs from rho.petasl
rho.conditional <- 0.5 ## substantially differs from  (rho.petasl-rhoLag.petasl)/sqrt((1-rho.asl)*(1-rho.pet))


## ** seed
seed <- 12

## ** display
round.dotplot <- 2
xspace.dotplot <- 0.3
size.dotplot <- c(6,4)

round.factor <- 1+9*(round.dotplot %% 1)/10
round.dotplot <- floor(round.dotplot)

region.name <- readRDS("results/regionName.rds")
region.keep <- region.name[names(mu.pet)][c(3,6,12,18)]
id.keep <- c(1:3,24)

color.region <- setNames(rep(palette.colors()[1],length(region.name)), region.name)
color.region[region.keep] <- palette.colors()[c(2,6,4,8)]
shape.region <- setNames(rep(1,length(region.name)), region.name)
shape.region[region.keep] <- c(8,13,16,17)
shape.region[setdiff(names(shape.region),region.keep)] <- setdiff(1:18,shape.region[region.keep])

## * Data
## ** simulate
dtL.scenarioA <- simData(seed = seed, n.obs = 24,
                         mu.PET = rep(mean(mu.pet), length(mu.pet)),
                         mu.fMRI = rep(mean(mu.asl), length(mu.asl)),
                         rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = rho.marginal, rho.conditional = rho.conditional,
                         sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)
dtL.scenarioB <- simData(seed = seed, n.obs = 24,
                         mu.PET = mu.pet,
                         mu.fMRI = mu.asl,
                         rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = 0, rho.conditional = 0,
                         sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)
dtL.scenarioC <- simData(seed = seed, n.obs = 24,
                         mu.PET = mu.pet,
                         mu.fMRI = mu.asl,
                         rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = rho.marginal, rho.conditional = rho.conditional,
                         sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)


## M.res <- rbind(A=runCor(PET + fMRI ~ region|id, data = dtL.scenarioA)[type %in% c("marginal","unclear")]$estimate,
##                B=runCor(PET + fMRI ~ region|id, data = dtL.scenarioB)[type %in% c("marginal","unclear")]$estimate,
##                C=runCor(PET + fMRI ~ region|id, data = dtL.scenarioC)[type %in% c("marginal","unclear")]$estimate)[,c(1,2,4)]
## print(M.res)
## A 0.5552999  0.279276644  0.27951793
## B 0.7007177 -0.009044682 -0.01830424
## C 0.7032253  0.279276644  0.27951793

## ** normalize
dtL.scenarioA[, c("PET.norm","fMRI.norm") := .(scale(PET),scale(fMRI)), by = "region"]
dtL.scenarioB[, c("PET.norm","fMRI.norm") := .(scale(PET),scale(fMRI)), by = "region"]
dtL.scenarioC[, c("PET.norm","fMRI.norm") := .(scale(PET),scale(fMRI)), by = "region"]

dtL.scenarioA[, region := region.name[region]]
dtL.scenarioB[, region := region.name[region]]
dtL.scenarioC[, region := region.name[region]]

## * Figure

## ** sub-figure 1
dtL1.scenarioA <- rbind(dtL.scenarioA[id %in% id.keep,.(region,id,fMRI,PET)],
                        dtL.scenarioA[,.(id = "Average", fMRI = mean(fMRI), PET = mean(PET)), by = "region"])
dtL1.scenarioB <- rbind(dtL.scenarioB[id %in% c(1:3,24),.(region,id,fMRI,PET)],
                        dtL.scenarioB[,.(id = "Average", fMRI = mean(fMRI), PET = mean(PET)), by = "region"])
dtL1.scenarioC <- rbind(dtL.scenarioC[id %in% c(1:3,24),.(region,id,fMRI,PET)],
                        dtL.scenarioC[,.(id = "Average", fMRI = mean(fMRI), PET = mean(PET)), by = "region"])
dtL1.scenario <- rbind(cbind(scenario = "A", dtL1.scenarioA),
                       cbind(scenario = "B", dtL1.scenarioB),
                       cbind(scenario = "C", dtL1.scenarioC))
dtL1.scenario$id2 <- factor(dtL1.scenario$id,
                            levels = c(id.keep,"Average"),
                            labels = c(paste0("Patient ",id.keep[1:2]),"...",paste0("Patient ",id.keep[4]),"average"))
dtL1.scenario$scenario2 <- paste("Scenario", dtL1.scenario$scenario, sep = " ")

figure3.1 <- ggplot(dtL1.scenario, aes(x = PET, y = fMRI))
figure3.1 <- figure3.1 + geom_point(aes(shape=region,color=region,size=region)) + facet_grid(scenario2~id2)
figure3.1 <- figure3.1 + scale_shape_manual(values=shape.region,breaks=names(shape.region))
figure3.1 <- figure3.1 + scale_size_manual(values=2+region.name %in% region.keep, breaks=region.name)
figure3.1 <- figure3.1 + scale_color_manual(values=color.region,breaks=names(shape.region))
figure3.1 <- figure3.1 + labs(x = '[11C]Cimbi-36 BPND', y = 'Baseline CBF (ml/100g/min)',
                              shape = "Region", color = "Region", size = "Region")
figure3.1 <- figure3.1 + theme(text = element_text(size=12), 
                               axis.line = element_line(linewidth = 1.25),
                               axis.ticks = element_line(linewidth = 2),
                               axis.ticks.length=unit(.25, "cm"),
                               legend.key.size = unit(1,"line"))
figure3.1

## ** sub-figure 2
dtL2.scenarioA <- rbind(dtL.scenarioA[id %in% c(1:3,24),.(region,id,fMRI.norm,PET.norm)],
                        dtL.scenarioA[,.(id = "Average", fMRI.norm = mean(fMRI.norm), PET.norm = mean(PET.norm)), by = "region"])
dtL2.scenarioB <- rbind(dtL.scenarioB[id %in% c(1:3,24),.(region,id,fMRI.norm,PET.norm)],
                        dtL.scenarioB[,.(id = "Average", fMRI.norm = mean(fMRI.norm), PET.norm = mean(PET.norm)), by = "region"])
dtL2.scenarioC <- rbind(dtL.scenarioC[id %in% c(1:3,24),.(region,id,fMRI.norm,PET.norm)],
                        dtL.scenarioC[,.(id = "Average", fMRI.norm = mean(fMRI.norm), PET.norm = mean(PET.norm)), by = "region"])

dtL2.scenario <- rbind(cbind(scenario = "A", dtL2.scenarioA),
                       cbind(scenario = "B", dtL2.scenarioB),
                       cbind(scenario = "C", dtL2.scenarioC))
dtL2.scenario$id2 <- factor(dtL2.scenario$id,
                            levels = c(1,2,3,24,"Average"),
                            labels = c("Patient 1","Patient 2","...","Patient 24","Average"))
dtL2.scenario$scenario2 <- paste("Scenario", dtL2.scenario$scenario, sep = " ")

figure3.2 <- ggplot(dtL2.scenario, aes(x = PET.norm, y = fMRI.norm))
figure3.2 <- figure3.2 + geom_point(aes(shape=region,color=region,size=region)) + facet_grid(scenario2~id2)
figure3.2 <- figure3.2 + scale_shape_manual(values=shape.region,breaks=names(shape.region))
figure3.2 <- figure3.2 + scale_size_manual(values=2+region.name %in% region.keep, breaks=region.name)
figure3.2 <- figure3.2 + scale_color_manual(values=color.region,breaks=names(shape.region))
figure3.2 <- figure3.2 + labs(x = 'Normalized [11C]Cimbi-36 BPND', y = 'Normalized baseline CBF',
                              shape = "Region", color = "Region", size = "Region")
figure3.2 <- figure3.2 + theme(text = element_text(size=12), 
                               axis.line = element_line(linewidth = 1.25),
                               axis.ticks = element_line(linewidth = 2),
                               axis.ticks.length=unit(.25, "cm"),
                               legend.key.size = unit(1,"line"))
figure3.2
    
## ** sub-figure 3
dtLC.scenarioA <- rbind(dtL.scenarioA[region %in% region.keep])
dtLC.scenarioB <- rbind(dtL.scenarioB[region %in% region.keep])
dtLC.scenarioC <- rbind(dtL.scenarioC[region %in% region.keep])
dtLC.scenario <- rbind(cbind(scenario = "A", dtLC.scenarioA),
                       cbind(scenario = "B", dtLC.scenarioB),
                       cbind(scenario = "C", dtLC.scenarioC))

figure3.3a <- ggplot(dtLC.scenario, aes(x = PET, y = fMRI))
figure3.3a <- figure3.3a + geom_point(size = 2)
figure3.3a <- figure3.3a + geom_point(size = 5, data = dtLC.scenario[id %in% id.keep])
figure3.3a <- figure3.3a + geom_text(data = dtLC.scenario[id %in% id.keep], aes(label = id), color = "white", size = 3)
figure3.3a <- figure3.3a + geom_smooth(method = "lm", se = FALSE, aes(color = region), linewidth = 1.2)
figure3.3a <- figure3.3a + facet_grid(scenario~region)
figure3.3a <- figure3.3a + theme(strip.text.y = element_blank())
figure3.3a <- figure3.3a + scale_color_manual(values = color.region[region.keep], breaks = region.keep)
figure3.3a <- figure3.3a + guides(color = "none")
figure3.3a <- figure3.3a + labs(x = '[11C]Cimbi-36 BPND', y = 'Baseline CBF (ml/100g/min)')
figure3.3a
 
dtLCor.scenario <- rbind(dtL.scenarioA[,.(scenario = "A", estimand = "Correlation", cor=cor(PET,fMRI)),by="region"],
                         dtL.scenarioB[,.(scenario = "B", estimand = "Correlation", cor=cor(PET,fMRI)),by="region"],
                         dtL.scenarioC[,.(scenario = "C", estimand = "Correlation", cor=cor(PET,fMRI)),by="region"])
dtLCor.scenario[, scenario2 := paste("Scenario", scenario, sep = " ")]
dtLCor.scenario[, group := round(round.factor*cor,round.dotplot)]
dtLCor.scenario[, myx := seq(0, by = xspace.dotplot, length.out = .N) - (.N-1)*xspace.dotplot/2, by = c("group","scenario")]

figure3.3b <- ggplot(dtLCor.scenario, aes(y = cor))
figure3.3b <- figure3.3b + geom_boxplot() + geom_point(aes(x=myx, color = region, size = region))
figure3.3b <- figure3.3b + geom_point(data = dtLCor.scenario[region %in% region.keep],
                                      aes(x=myx, color = region, size = region))
figure3.3b <- figure3.3b + xlab("")
figure3.3b <- figure3.3b + facet_grid(scenario2~estimand)
figure3.3b <- figure3.3b + scale_color_manual(values = color.region, breaks = names(color.region))
figure3.3b <- figure3.3b + scale_size_manual(values= c(2,3)[region.name %in% region.keep+1], breaks = names(color.region))
figure3.3b <- figure3.3b + guides(shape = "none", size = "none", color = "none")
figure3.3b <- figure3.3b + theme(axis.title.y=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.x=element_text(color="white"))
figure3.3b

figure3.3 <- ggarrange(figure3.3a, figure3.3b + ylab(""), nrow = 1, widths = c(1.5,0.5))
figure3.3

## * export
pdf("./figures/figure3.1.pdf", width = 9, height = 6)
figure3.1
dev.off()

png("./figures/figure3.1.png", width = 9, height = 6, units = "in", res = 500)
figure3.1
dev.off()

pdf("./figures/figure3.2.pdf", width = 9, height = 6)
figure3.2
dev.off()

png("./figures/figure3.2.png", width = 9, height = 6, units = "in", res = 500)
figure3.2
dev.off()

pdf("./figures/figure3.3.pdf", width = 9, height = 6)
figure3.3
dev.off()

png("./figures/figure3.3.png", width = 9, height = 6, units = "in", res = 500)
figure3.3
dev.off()
##----------------------------------------------------------------------
### figure3.R ends here
