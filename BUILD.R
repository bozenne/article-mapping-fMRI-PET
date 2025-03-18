### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 18 2025 (16:26) 
## Version: 
## Last-Updated: mar 18 2025 (18:33) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(data.table)
source("FCT.R")

## * Collect and process results
dt.simABC <- loadRes("results/simulationABC")
dt.simDE <- loadRes("results/simulationDE")

dtS.sim <- rbind(dt.simABC,dt.simDE)[, .(.N, n.NA = sum(is.na(estimate)),
                                         truth = unique(truth),
                                         estimate = mean(estimate, na.rm = TRUE),
                                         std = sd(estimate, na.rm = TRUE),
                                         coverage = mean(between(truth, lower, upper), na.rm = TRUE),
                                         rejection.rate = mean(p.value<=0.05, na.rm = TRUE) ),
                                     by = c("scenario","n","method","type")]

## * Export
saveRDS(dtS.sim, file = "results/simulation-summary.rds")

##----------------------------------------------------------------------
### BUILD.R ends here
