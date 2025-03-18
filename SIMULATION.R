### SIMULATION.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  7 2025 (11:48) 
## Version: 
## Last-Updated: mar 18 2025 (13:39) 
##           By: Brice Ozenne
##     Update #: 56
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code: cd /projects/biostat01/people/hpl802/article-mapping-fMRI-PET/


library(LMMstar)
library(pbapply)
library(data.table)
source("FCT.R")
source("cor.testIID.R")

## * Setting
n.sim <- 1000 ## 1e3
n.cpus <- 100 ## 

param.lmm <- readRDS("results/strategy3-lmm-param.rds")
param.lmmH <- readRDS("results/strategy3-lmmH-param.rds")

muRegion.asl <- c("amygdala" = 41.75461171, "antcin" = 63.97868936, "cau" = 51.03709735, "DLPFC" = 61.81180998, "entorhinc" = 37.93286726, "hippocampus" = 48.20552864, "ins" = 61.01470806, "medinffg" = 64.47818102, "occ" = 49.06232099, "orbfrc" = 54.13080092, "parc" = 57.03041177, "postcin" = 56.59856003, "put" = 53.58549712, "senmotc" = 57.3037015, "supfg" = 56.45369853, "suptempg" = 59.07503085, "th" = 47.54133327, "VLPFC" = 64.91376706) ## param.lmm$mu.asl
muRegion.pet <- c("amygdala" = 0.83208117, "antcin" = 1.56758472, "cau" = 0.24240839, "DLPFC" = 1.20843445, "entorhinc" = 0.54770079, "hippocampus" = 0.55927805, "ins" = 1.57503343, "medinffg" = 1.27606528, "occ" = 1.15537147, "orbfrc" = 1.44086414, "parc" = 1.1357671, "postcin" = 1.36909301, "put" = 0.35256114, "senmotc" = 0.84761643, "supfg" = 1.13908271, "suptempg" = 1.47829904, "th" = 0.42062905, "VLPFC" = 1.3450188) ## param.lmm$mu.pet
mu.asl <- mean(muRegion.asl)  
mu.pet <- mean(muRegion.pet) 
sigma.asl <- 9.270355  ## param.lmm$sigma.asl
sigma.pet <- 0.2546338  ## param.lmm$sigma.pet

rho.asl <- 0.4815182 ## param.lmm$rho.asl
rho.pet <- 0.6382445 ## param.lmm$rho.pet
rho.petasl <- 0.1381907 ## param.lmm$rho.petasl
rhoLag.petasl <- 0.1049642  ## param.lmm$rhoLag.petasl

rho.marginal <- rho.petasl
rho.conditional <- (rho.petasl-rhoLag.petasl)/sqrt((1-rho.asl)*(1-rho.pet))

## * Run simulation 

## ** Scenario A
lsA.sim <- pblapply(1:n.sim, function(iSim){ ## iSim <- 1
    iData26 <- simData(seed = iSim, n.obs = 26,
                       mu.PET = rep(mu.pet, length(muRegion.pet)),
                       mu.fMRI = rep(mu.asl, length(muRegion.asl)),
                       rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = 0.25, rho.conditional = 0.5,
                       sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)
    iRes26 <- runCor(PET + fMRI ~ region|id, data = iData26)[type!="latent"]

    iData1000 <- simData(seed = iSim, n.obs = 1000,
                         mu.PET = rep(mu.pet, length(muRegion.pet)),
                         mu.fMRI = rep(mu.asl, length(muRegion.asl)),
                         rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = 0.25, rho.conditional = 0.5,
                         sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)
    iRes1000 <- runCor(PET + fMRI ~ region|id, data = iData1000)[type!="latent"]

    iOut <- rbind(cbind(seed = iSim, scenario = "A", n = 26, iRes26),
                  cbind(seed = iSim, scenario = "A", n = 1000, iRes1000))
    iOut$truth <- c(unclear = 0.5, marginal = 0.25, conditional = 0.5)[iOut$type]
    return(iOut)
}, cl = n.cpus)

## ** Scenario B
lsB.sim <- pblapply(1:n.sim, function(iSim){
    iData26 <- simData(seed = iSim, n.obs = 26,
                       mu.PET = muRegion.pet,
                       mu.fMRI = muRegion.asl,
                       rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = 0, rho.conditional = 0,
                       sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)
    iRes26 <- runCor(PET + fMRI ~ region|id, data = iData26)[type!="latent"]

    iData1000 <- simData(seed = iSim, n.obs = 1000,
                         mu.PET = muRegion.pet,
                         mu.fMRI = muRegion.asl,
                         rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = 0, rho.conditional = 0,
                         sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)
    iRes1000 <- runCor(PET + fMRI ~ region|id, data = iData1000)[type!="latent"]

    iOut <- rbind(cbind(seed = iSim, scenario = "B", n = 26, iRes26),
                  cbind(seed = iSim, scenario = "B", n = 1000, iRes1000))
    iOut$truth <- c(unclear = 0, marginal = 0, conditional = 0)[iOut$type]
    return(iOut)
}, cl = n.cpus)

## ** Scenario C
lsC.sim <- pblapply(1:n.sim, function(iSim){
    iData26 <- simData(seed = iSim, n.obs = 26,
                       mu.PET = muRegion.pet,
                       mu.fMRI = muRegion.asl,
                       rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = 0.25, rho.conditional = 0.5,
                       sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)
    iRes26 <- runCor(PET + fMRI ~ region|id, data = iData26)[type!="latent"]

    iData1000 <- simData(seed = iSim, n.obs = 1000,
                         mu.PET = muRegion.pet,
                         mu.fMRI = muRegion.asl,
                         rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = 0.25, rho.conditional = 0.5,
                         sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)
    iRes1000 <- runCor(PET + fMRI ~ region|id, data = iData1000)[type!="latent"]

    iOut <- rbind(cbind(seed = iSim, scenario = "C", n = 26, iRes26),
                  cbind(seed = iSim, scenario = "C", n = 1000, iRes1000))
    iOut$truth <- c(unclear = 0.5, marginal = 0.25, conditional = 0.5)[iOut$type]
    return(iOut)
}, cl = n.cpus)

## ** Scenario D (extra scenario considering model misspecification - heteroschedasticity)
lsD.sim <- pblapply(1:n.sim, function(iSim){
    iData26 <- simData(seed = iSim, n.obs = 26,
                       mu.PET = param.lmmH$mu.pet,
                       mu.fMRI = param.lmmH$mu.asl,
                       rho.PET = param.lmmH$rho.pet, rho.fMRI = param.lmmH$rho.asl, rho.marginal = param.lmmH$rho.petasl,
                       rho.conditional = (param.lmmH$rho.petasl-param.lmmH$rhoLag.petasl)/sqrt((1-param.lmmH$rho.asl)*(1-param.lmmH$rho.pet)),
                       sigma.PET = param.lmmH$sigma.pet, sigma.fMRI = param.lmmH$sigma.asl)
    iRes26 <- runCor(PET + fMRI ~ region|id, data = iData26)[type!="latent"]

    iData1000 <- simData(seed = iSim, n.obs = 1000,
                         mu.PET = param.lmmH$mu.pet,
                         mu.fMRI = param.lmmH$mu.asl,
                         rho.PET = param.lmmH$rho.pet, rho.fMRI = param.lmmH$rho.asl, rho.marginal = param.lmmH$rho.petasl,
                         rho.conditional = (param.lmmH$rho.petasl-param.lmmH$rhoLag.petasl)/sqrt((1-param.lmmH$rho.asl)*(1-param.lmmH$rho.pet)),
                         sigma.PET = param.lmmH$sigma.pet, sigma.fMRI = param.lmmH$sigma.asl)
    iRes1000 <- runCor(PET + fMRI ~ region|id, data = iData1000)[type!="latent"]

    iOut <- rbind(cbind(seed = iSim, scenario = "D", n = 26, iRes26),
                  cbind(seed = iSim, scenario = "D", n = 1000, iRes1000))
    iOut$truth <- c(unclear = 0.5, marginal = 0.25, conditional = 0.5)[iOut$type]
    return(iOut)
}, cl = n.cpus)

## * export results 
dt.sim <- as.data.table(do.call(rbind,c(lsA.sim,lsB.sim,lsC.sim,lsD.sim)))[method %in% c("averageSignal","averageCor","averageIdNorm","lmm")]
saveRDS(dt.sim, file = "results/sim-scenarioABC_n26_1000.rds")

##----------------------------------------------------------------------
### SIMULATION.R ends here
