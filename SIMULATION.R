### SIMULATION.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2025 (11:43) 
## Version: 
## Last-Updated: jan 23 2025 (11:39) 
##           By: Brice Ozenne
##     Update #: 8
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(LMMstar)
source("FCT.R")

## * Setting
n.region <- 18


## * Analysis

## ** Scenario A
dtL.scenarioA <- simData(seed = 1, 
                         mu.PET = rep(4.8,n.region),
                         mu.fMRI = 1+rep(5,n.region),
                         rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5,
                         sigma.PET = 3, sigma.fMRI = 5)

runCor(PET + fMRI ~ region|id, data = dtL.scenarioA)[type=="conditional"]
##           method        type   estimate      lower     upper      p.value
##           <char>      <char>      <num>      <num>     <num>        <num>
## 1: averageSignal conditional 0.35523747 -0.1338279 0.7051574 1.480007e-01
## 2:    averageCor conditional 0.06724627 -0.1299778 0.2593515 5.051268e-01
## 3:     averageId conditional 0.46098596  0.4026315 0.5156086 0.000000e+00
## 4: averageIdNorm conditional 0.46194101  0.4039113 0.5162680 0.000000e+00
## 5:          lmm0 conditional 0.45521698  0.3727159 0.5305780 3.729945e-09
## 6:           lmm conditional 0.45818886  0.3742144 0.5347101 4.781979e-09
plotCor(dtL.scenarioA, type.norm = 2)


## ** Scenario B
dtL.scenarioB <- simData(seed = 1,
                         mu.PET = 1:n.region/2,
                         mu.fMRI = 1+5+(1:n.region)/5,
                         rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0, rho.conditional = 0)

runCor(PET + fMRI ~ region|id, data = dtL.scenarioB)[type=="conditional"]
##           method        type    estimate      lower       upper      p.value
##           <char>      <char>       <num>      <num>       <num>        <num>
## 1: averageSignal conditional  0.87715246  0.6948242 0.953534624 1.761851e-06
## 2:    averageCor conditional -0.19241825 -0.3771513 0.007039887 5.854188e-02
## 3:     averageId conditional  0.27814586  0.2019486 0.350996205 4.508172e-12
## 4: averageIdNorm conditional -0.04499694 -0.1173474 0.027828671 2.258184e-01
## 5:          lmm0 conditional  0.27395755  0.1830694 0.360202276 5.260763e-07
## 6:           lmm conditional -0.06950383 -0.1699130 0.032333633 1.685021e-01

plotCor(dtL.scenarioB, type.norm = 2)

## ** Scenario C
dtL.scenarioC <- simData(seed = 1,
                         mu.PET = 1:n.region/2,
                         mu.fMRI = 1+5+(1:n.region)/5,
                         rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5)

runCor(PET + fMRI ~ region|id, data = dtL.scenarioC)[type=="conditional"]
##           method        type   estimate      lower     upper      p.value
##           <char>      <char>      <num>      <num>     <num>        <num>
## 1: averageSignal conditional 0.88649462  0.7159336 0.9572006 9.647869e-07
## 2:    averageCor conditional 0.07516425 -0.1165997 0.2615253 4.430905e-01
## 3:     averageId conditional 0.35941359  0.2871448 0.4276139 0.000000e+00
## 4: averageIdNorm conditional 0.47297607  0.4168770 0.5254856 0.000000e+00
## 5:          lmm0 conditional 0.35962481  0.2733142 0.4401911 7.697434e-10
## 6:           lmm conditional 0.46933238  0.3863745 0.5447577 2.957791e-09

plotCor(dtL.scenarioC, type.norm = 2)

## * Figure simulation
figCor <- plotCor(list("Scenario A" = dtL.scenarioA,"Scenario B" = dtL.scenarioB, "Scenario C" = dtL.scenarioC),
                  type.norm = 2, xspace.dotplot = 0.15, round.dotplot = 1.9, size.dotplot = c(3,2))

pdf("./figures/fig-correlationEstimators.pdf", width = 12, height = 8)
figCor
dev.off()

png("./figures/fig-correlationEstimators.png", width = 12, height = 8, units = "in", res = 500)
figCor
dev.off()
##----------------------------------------------------------------------
### SIMULATION.R ends here
