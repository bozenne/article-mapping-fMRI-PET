## * Header 
## path <- "c:/Users/hpl802/Documents/Github/article-mapping-fMRI-PET/"
## setwd(path)
## source("BATCH_simulationABC.R")
## sbatch -a 1-1 -J 'simulationABC' --output=/dev/null --error=/dev/null R CMD BATCH --vanilla BATCH_simulationABC.R /dev/null 

rm(list = ls())
gc()

## * number of jobs
iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
if(is.na(iter_sim)){iter_sim <- 1}
if(is.na(n.iter_sim)){n.iter_sim <- 10}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

## * prepare export
path <- "."
path.res <- file.path(path,"results","simulationABC")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"results"))==FALSE){
    dir.create(file.path(path,"results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulationABC")
if(dir.exists(path.output)==FALSE){
    if(dir.exists(file.path(path,"output"))==FALSE){
    dir.create(file.path(path,"output"))
    }
    dir.create(path.output)
}

## * libraries
library(LMMstar)
library(data.table)
source("FCT.R")
source("cor.testIID.R")

## * settings

## ** number of simulations
n.sim <- 100

## ** seeds
nsimAll <- n.sim * n.iter_sim
allSeeds <- 1:nsimAll
currentSeeds <- (iter_sim-1)*n.sim+(1:n.sim)

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

grid.scenario <-  expand.grid(n = c(24,1000), scenario = c("A","B","C"), stringsAsFactors = FALSE)
n.grid <- NROW(grid.scenario)

## * function to execute
res <- NULL
for(iSim in 1:n.sim){ ## iSim <- 1
    iSeed <- currentSeeds[iSim]
    cat("simulation ",iSim," (seed = ",iSeed,"): ", sep = "")

    for(iG in 1:n.grid){ ## iG <- 1
        
        iN.obs <- grid.scenario[iG,"n"]
        iScenario <- grid.scenario[iG,"scenario"]
        cat(iG,"(n=",iN.obs,",",iScenario,") ", sep = "")
        
        iMu.PET <- switch(grid.scenario[iG,"scenario"],
                          "A" = rep(mean(mu.pet), length(mu.pet)),
                          "B" = mu.pet,
                          "C" = mu.pet)
        iMu.fMRI <- switch(grid.scenario[iG,"scenario"],
                           "A" = rep(mean(mu.asl), length(mu.asl)),
                           "B" = mu.asl,
                           "C" = mu.asl)
        iRho.marginal <- ifelse(iScenario=="B", 0, rho.marginal)
        iRho.conditional <- ifelse(iScenario=="B", 0, rho.conditional)

        iData <- simData(seed = iSeed, n.obs = iN.obs,
                         mu.PET = iMu.PET,
                         mu.fMRI = iMu.fMRI,
                         rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = iRho.marginal, rho.conditional = iRho.conditional,
                         sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)

        iRes <- cbind(seed = iSeed, n = iN.obs, scenario = iScenario, 
                      runCor(PET + fMRI ~ region|id, data = iData)[type!="latent"])
        iRes$truth <- c(unclear = iRho.conditional, marginal = iRho.marginal, conditional = iRho.conditional)[iRes$type]

        res <- rbind(res, iRes)
    }
    cat("\n")
    saveRDS(res, file = file.path(path.res,paste0("simul_",iter_sim,"(tempo).rds")))
}

## * export
saveRDS(res, file = file.path(path.res,paste0("simul_",iter_sim,".rds")))

## * R version
print(sessionInfo())

	
