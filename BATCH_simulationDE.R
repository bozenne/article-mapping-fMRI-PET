## * Header 
## path <- "c:/Users/hpl802/Documents/Github/article-mapping-fMRI-PET/"
## setwd(path)
## source("BATCH_simulationDE.R")
## sbatch -a 1-1 -J 'simulationDE' --output=/dev/null --error=/dev/null R CMD BATCH --vanilla BATCH_simulationDE.R /dev/null 

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
path.res <- file.path(path,"results","simulationDE")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"results"))==FALSE){
    dir.create(file.path(path,"results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","simulationDE")
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
## param.lmmH <- readRDS("results/strategy3-lmmH-param.rds")
mu.asl <- c("amygdala" = 41.75461171, "antcin" = 63.97868936, "cau" = 51.03709735, "DLPFC" = 61.81180998, "entorhinc" = 37.93286726, "hippocampus" = 48.20552864, "ins" = 61.01470806, "medinffg" = 64.47818102, "occ" = 49.06232099, "orbfrc" = 54.13080092, "parc" = 57.03041177, "postcin" = 56.59856003, "put" = 53.58549712, "senmotc" = 57.3037015, "supfg" = 56.45369853, "suptempg" = 59.07503085, "th" = 47.54133327, "VLPFC" = 64.91376706) ## param.lmmH$mu.asl
mu.pet <- c("amygdala" = 0.83845228, "antcin" = 1.56758472, "cau" = 0.2503588, "DLPFC" = 1.20843445, "entorhinc" = 0.54554437, "hippocampus" = 0.5574931, "ins" = 1.57503343, "medinffg" = 1.27606528, "occ" = 1.15442359, "orbfrc" = 1.44086414, "parc" = 1.1357671, "postcin" = 1.36909301, "put" = 0.34928437, "senmotc" = 0.84761643, "supfg" = 1.13908271, "suptempg" = 1.47829904, "th" = 0.42230255, "VLPFC" = 1.3450188) ## param.lmmH$mu.pet
sigma.asl <- c("amygdala" = 8.99965093, "antcin" = 8.09495239, "cau" = 13.20924235, "DLPFC" = 7.69483976, "entorhinc" = 10.67055643, "hippocampus" = 6.35501704, "ins" = 9.1784132, "medinffg" = 7.4033057, "occ" = 9.10531456, "orbfrc" = 10.5240258, "parc" = 6.82659697, "postcin" = 12.80013308, "put" = 9.62445417, "senmotc" = 6.86834497, "supfg" = 7.5001178, "suptempg" = 7.27306258, "th" = 16.99779638, "VLPFC" = 9.12007523)  ## param.lmmH$sigma.asl
sigma.pet <- c("amygdala" = 0.31454483, "antcin" = 0.29756436, "cau" = 0.14952683, "DLPFC" = 0.2476599, "entorhinc" = 0.27166285, "hippocampus" = 0.18168736, "ins" = 0.28450984, "medinffg" = 0.2379356, "occ" = 0.20805809, "orbfrc" = 0.2928512, "parc" = 0.22924271, "postcin" = 0.39920728, "put" = 0.12105493, "senmotc" = 0.1881363, "supfg" = 0.24928567, "suptempg" = 0.28677391, "th" = 0.10523051, "VLPFC" = 0.27782552)  ## param.lmmH$sigma.pet

rho.asl <- 0.55434116 ## param.lmmH$rho.asl
rho.pet <- 0.6843263  ## param.lmmH$rho.pet
rho.petasl <- 0.1353023  ## param.lmmH$rho.petasl
rhoLag.petasl <- 0.1137632   ## param.lmmH$rhoLag.petasl

rho.marginal <- rho.petasl
rho.conditional <- (rho.petasl-rhoLag.petasl)/sqrt((1-rho.asl)*(1-rho.pet))

grid.scenario <-  expand.grid(n = c(24,1000), scenario = c("D","E"), stringsAsFactors = FALSE)
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
                          "D" = rep(mean(mu.pet), length(mu.pet)),
                          "E" = mu.pet)
        iMu.fMRI <- switch(grid.scenario[iG,"scenario"],
                           "D" = rep(mean(mu.asl), length(mu.asl)),
                           "E" = mu.asl)
        iRho.marginal <- ifelse(iScenario=="D", 0, rho.marginal)
        iRho.conditional <- ifelse(iScenario=="D", 0, rho.conditional)

        iData <- simData(seed = iSeed, n.obs = iN.obs,
                         mu.PET = iMu.PET,
                         mu.fMRI = iMu.fMRI,
                         rho.PET = rho.pet, rho.fMRI = rho.asl, rho.marginal = iRho.marginal, rho.conditional = iRho.conditional,
                         sigma.PET = sigma.pet, sigma.fMRI = sigma.asl)

        iRes <- cbind(seed = iSeed, n = iN.obs, scenario = iScenario, 
                      runCor(PET + fMRI ~ region|id, data = iData, lmmH = FALSE)[type!="latent"])

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

	
