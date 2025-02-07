### SIMULATION-rep.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  7 2025 (11:48) 
## Version: 
## Last-Updated: feb  7 2025 (19:17) 
##           By: Brice Ozenne
##     Update #: 25
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


library(LMMstar)
library(pbapply)
library(data.table)
source("FCT.R")
source("cor.testIID.R")

format.pval2 <- function(pv, ...){
    ## WARNING: no difference between 0<x<0.001 and -0.001<x<0 (both are coded <0.001)
    sign <- sign(pv)
    out <- format.pval(pv = abs(pv), ...)
    out[sign<0 & !grepl("<",out)] <- paste0("-",out[sign<0 & !grepl("<",out)])
    return(out)
}

## * Setting
n.region <- 18
n.sim <- 1000 ## 1e3
n.cpus <- 110 ## 

## * Run simulation 

## ** Scenario A
lsA.sim <- pblapply(1:n.sim, function(iSim){ ## iSim <- 1
    iData26 <- simData(seed = iSim, n.obs = 26,
                       mu.PET = rep(4.8,n.region),
                       mu.fMRI = 1+rep(5,n.region),
                       rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5,
                       sigma.PET = 3, sigma.fMRI = 5)
    iRes26 <- runCor(PET + fMRI ~ region|id, data = iData26)[type!="latent"]

    iData1000 <- simData(seed = iSim, n.obs = 1000,
                       mu.PET = rep(4.8,n.region),
                       mu.fMRI = 1+rep(5,n.region),
                       rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5,
                       sigma.PET = 3, sigma.fMRI = 5)
    iRes1000 <- runCor(PET + fMRI ~ region|id, data = iData1000)[type!="latent"]

    iOut <- rbind(cbind(seed = iSim, scenario = "A", n = 26, iRes26),
                  cbind(seed = iSim, scenario = "A", n = 1000, iRes1000))
    iOut$truth <- c(unclear = 0.5, marginal = 0.25, conditional = 0.5)[iOut$type]
    return(iOut)
}, cl = n.cpus)

## ** Scenario B
lsB.sim <- pblapply(1:n.sim, function(iSim){
    iData26 <- simData(seed = iSim, n.obs = 26,
                       mu.PET = 1:n.region/2,
                       mu.fMRI = 1+5+(1:n.region)/5,
                       rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0, rho.conditional = 0)
    iRes26 <- runCor(PET + fMRI ~ region|id, data = iData26)[type!="latent"]

    iData1000 <- simData(seed = iSim, n.obs = 1000,
                         mu.PET = 1:n.region/2,
                         mu.fMRI = 1+5+(1:n.region)/5,
                         rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0, rho.conditional = 0)
    iRes1000 <- runCor(PET + fMRI ~ region|id, data = iData1000)[type!="latent"]

    iOut <- rbind(cbind(seed = iSim, scenario = "B", n = 26, iRes26),
                  cbind(seed = iSim, scenario = "B", n = 1000, iRes1000))
    iOut$truth <- c(unclear = 0, marginal = 0, conditional = 0)[iOut$type]
    return(iOut)
}, cl = n.cpus)

## ** Scenario C
lsC.sim <- pblapply(1:n.sim, function(iSim){
    iData26 <- simData(seed = iSim, n.obs = 26,
                       mu.PET = 1:n.region/2,
                       mu.fMRI = 1+5+(1:n.region)/5,
                       rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5)
    iRes26 <- runCor(PET + fMRI ~ region|id, data = iData26)[type!="latent"]

    iData1000 <- simData(seed = iSim, n.obs = 1000,
                         mu.PET = 1:n.region/2,
                         mu.fMRI = 1+5+(1:n.region)/5,
                         rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5)
    iRes1000 <- runCor(PET + fMRI ~ region|id, data = iData1000)[type!="latent"]

    iOut <- rbind(cbind(seed = iSim, scenario = "C", n = 26, iRes26),
                  cbind(seed = iSim, scenario = "C", n = 1000, iRes1000))
    iOut$truth <- c(unclear = 0.5, marginal = 0.25, conditional = 0.5)[iOut$type]
    return(iOut)
}, cl = n.cpus)


## * export results 
dt.sim <- as.data.table(do.call(rbind,c(lsA.sim,lsB.sim,lsC.sim)))[method %in% c("averageSignal","averageCor","averageIdNorm","lmm")]
saveRDS(dt.sim, file = "results/sim-scenarioABC_n26_1000.rds")

## * process results
dt.sim <-  readRDS(file = "results/sim-scenarioABC_n26_1000.rds")
## dt.sim$truth <- 0 ## scenario B
## dt.sim[scenario == "A", truth := c(unclear = 0.5, marginal = 0.25, conditional = 0.5)[.SD$type]]
## dt.sim[scenario == "C", truth := c(unclear = 0.5, marginal = 0.25, conditional = 0.5)[.SD$type]]
dt.sim[, strategy := as.character(factor(type, levels = c("unclear","marginal","conditional"), labels = c("1","2","3")))]
dt.sim[method == "lmm", strategy := paste0(strategy, " (LMM)")]
dt.sim[, strategy := factor(strategy, levels = c("1","2","2 (LMM)", "3", "3 (LMM)"),
                            labels = c("1      ","2      ","2 (LMM)", "3      ", "3 (LMM)"))]

dtSL.sim <- dt.sim[, .(.N, bias = mean(truth-estimate), coverage = mean(between(truth, lower, upper)), rejection.rate = mean(p.value<=0.05) ), by = c("scenario","n","strategy")]
dtSW.sim <- dcast(dtSL.sim, scenario+strategy~n, value.var = c("bias","coverage","rejection.rate"))
## options(width = 140)
table.sim <- cbind(dtSW.sim[,.(scenario,strategy)],
                   Bias = paste0(format.pval2(dtSW.sim$bias_26, digits = 1, eps = 1e-3), "/",format.pval2(dtSW.sim$bias_1000, digits = 1, eps = 1e-3)),
                   Coverage = paste0(formatC(dtSW.sim$coverage_26, digits = 3, format = "f"), "/",formatC(dtSW.sim$coverage_1000, digits = 3, format = "f")),
                   "Type 1 error" = paste0(formatC(dtSW.sim$rejection.rate_26, digits = 3, format = "f"), "/",formatC(dtSW.sim$rejection.rate_1000, digits = 3, format = "f")),
                   Power = paste0(formatC(dtSW.sim$rejection.rate_26, digits = 3, format = "f"), "/",formatC(dtSW.sim$rejection.rate_1000, digits = 3, format = "f"))
)
table.sim[scenario != "B", "Type 1 error" := ""]
table.sim[scenario == "B", "Power" := ""]
names(table.sim)[3:6] <- paste(names(table.sim)[3:6], "(n=26/n=1000)")
## WARNING: no difference between 0<x<0.001 and -0.001<x<0 (both are coded <0.001)

table.sim
##     scenario strategy Bias (n=26/n=1000) Coverage (n=26/n=1000) Type 1 error (n=26/n=1000) Power (n=26/n=1000)
##       <char>   <fctr>             <char>                 <char>                     <char>              <char>
##  1:        A  1              0.005/0.006            0.948/0.949                                    0.594/0.590
##  2:        A  2             0.005/<0.001            0.957/0.954                                    0.408/1.000
##  3:        A  2 (LMM)      <0.001/<0.001            0.955/0.955                                    0.440/1.000
##  4:        A  3              0.025/0.011            0.891/0.542                                    1.000/1.000
##  5:        A  3 (LMM)      <0.001/<0.001            0.958/0.955                                    1.000/1.000
##  6:        B  1            -0.869/-0.996            0.000/0.000                1.000/1.000                    
##  7:        B  2             0.001/<0.001            0.954/0.954                0.046/0.046                    
##  8:        B  2 (LMM)       0.001/<0.001            0.956/0.953                0.044/0.047                    
##  9:        B  3            -0.002/<0.001            0.946/0.951                0.054/0.049                    
## 10:        B  3 (LMM)      -0.002/<0.001            0.953/0.953                0.047/0.047                    
## 11:        C  1            -0.376/-0.496            0.043/0.000                                    1.000/1.000
## 12:        C  2             0.004/<0.001            0.954/0.956                                    0.413/1.000
## 13:        C  2 (LMM)      -0.002/<0.001            0.955/0.957                                    0.446/1.000
## 14:        C  3              0.025/0.011            0.894/0.547                                    1.000/1.000
## 15:        C  3 (LMM)      <0.001/<0.001            0.958/0.950                                    1.000/1.000


## relative bias
100*dtSW.sim$bias_26[1:5]/c(0.5,0.25,0.25,0.5,0.5)


0.005 / 0.25
##----------------------------------------------------------------------
### SIMULATION-rep.R ends here
