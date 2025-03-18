### figure3.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2025 (11:43) 
## Version: 
## Last-Updated: mar 18 2025 (11:46) 
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
library(ggplot2)
library(ggpubr)
source("FCT.R")
source("cor.testIID.R")

## * Setting
n.region <- 18


## * Analysis

## ** Scenario A
dtL.scenarioA <- simData(seed = 1, n.obs = 26,
                         mu.PET = rep(4.8,n.region),
                         mu.fMRI = 1+rep(5,n.region),
                         rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5,
                         sigma.PET = 3, sigma.fMRI = 5)

runCor(PET + fMRI ~ region|id, data = dtL.scenarioA)[type!="latent"]
##           method        type   estimate      lower     upper   p.value
##           <char>      <char>      <num>      <num>     <num>     <num>
## 1: averageSignal     unclear 0.35523747 -0.1338279 0.7051574 0.1480007
## 2:    averageCor    marginal 0.06724627 -0.1762087 0.3029487 0.5906723
## 3:     averageId conditional 0.46098596  0.4026315 0.5156086 0.0000000
## 4: averageIdNorm conditional 0.46194101  0.4039113 0.5162680 0.0000000
## 5:          lmm0    marginal 0.06699296 -0.2090322 0.3331310 0.6377102
## 6:          lmm0 conditional 0.45521698  0.3782534 0.5259304 0.0000000
## 7:           lmm    marginal 0.06895623 -0.2065305 0.3343131 0.6270641
## 8:           lmm conditional 0.45818886  0.3798878 0.5299707 0.0000000
plotCor(dtL.scenarioA, type.norm = 2)

dtL1000.scenarioA <- simData(seed = 1, n.obs = 1000,
                             mu.PET = rep(4.8,n.region),
                             mu.fMRI = 1+rep(5,n.region),
                             rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5,
                             sigma.PET = 3, sigma.fMRI = 5)
runCor(PET + fMRI ~ region|id, data = dtL1000.scenarioA)[type!="latent"]
##          method        type  estimate      lower     upper    p.value
##           <char>      <char>     <num>      <num>     <num>      <num>
## 1: averageSignal     unclear 0.5226332 0.07376088 0.7954195 0.02606444
## 2:    averageCor    marginal 0.2341701 0.19123069 0.2762150 0.00000000
## 3:     averageId conditional 0.4929088 0.48146203 0.5041875 0.00000000
## 4: averageIdNorm conditional 0.4920914 0.48065444 0.5033609 0.00000000
## 5:          lmm0    marginal 0.2344417 0.19136205 0.2766199 0.00000000
## 6:          lmm0 conditional 0.5025698 0.49124787 0.5137220 0.00000000
## 7:           lmm    marginal 0.2344234 0.19134019 0.2766051 0.00000000
## 8:           lmm conditional 0.5025487 0.49122076 0.5137068 0.00000000

## ** Scenario B
dtL.scenarioB <- simData(seed = 1,
                         mu.PET = 1:n.region/2,
                         mu.fMRI = 1+5+(1:n.region)/5,
                         rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0, rho.conditional = 0)

runCor(PET + fMRI ~ region|id, data = dtL.scenarioB)[type!="latent"]
##           method        type    estimate       lower      upper      p.value
##           <char>      <char>       <num>       <num>      <num>        <num>
## 1: averageSignal     unclear  0.87715246  0.69482415 0.95353462 1.761851e-06
## 2:    averageCor    marginal -0.19241825 -0.41409138 0.05080116 1.200991e-01
## 3:     averageId conditional  0.27814586  0.20194861 0.35099620 4.508172e-12
## 4: averageIdNorm conditional -0.04499694 -0.11734742 0.02782867 2.258184e-01
## 5:          lmm0    marginal  0.10685976 -0.01371068 0.22436726 8.224102e-02
## 6:          lmm0 conditional  0.27395755  0.18572707 0.35780567 3.413854e-09
## 7:           lmm    marginal -0.19472621 -0.43633655 0.07307638 1.528796e-01
## 8:           lmm conditional -0.06950383 -0.16321342 0.02544865 1.512295e-01
plotCor(dtL.scenarioB, type.norm = 2)

dtL1000.scenarioB <- simData(seed = 1, n.obs = 1000,
                             mu.PET = 1:n.region/2,
                             mu.fMRI = 1+5+(1:n.region)/5,
                             rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0, rho.conditional = 0)
runCor(PET + fMRI ~ region|id, data = dtL1000.scenarioB)[type!="latent"]
##           method        type     estimate       lower      upper      p.value
##           <char>      <char>        <num>       <num>      <num>        <num>
## 1: averageSignal     unclear  0.994426413  0.98473915 0.99797070 4.601269e-17
## 2:    averageCor    marginal -0.015713127 -0.05958441 0.02821874 4.833414e-01
## 3:     averageId conditional  0.317328451  0.30357429 0.33095039 0.000000e+00
## 4: averageIdNorm conditional  0.004084752 -0.01092864 0.01909630 5.938663e-01
## 5:          lmm0    marginal  0.189618749  0.17095615 0.20814533 0.000000e+00
## 6:          lmm0 conditional  0.321129336  0.30758230 0.33454621 0.000000e+00
## 7:           lmm    marginal -0.015554549 -0.05944429 0.02839522 4.879461e-01
## 8:           lmm conditional  0.003255548 -0.01178367 0.01829329 6.713760e-01

## ** Scenario C
dtL.scenarioC <- simData(seed = 1,
                         mu.PET = 1:n.region/2,
                         mu.fMRI = 1+5+(1:n.region)/5,
                         rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5)

runCor(PET + fMRI ~ region|id, data = dtL.scenarioC)[type!="latent"]
##           method        type   estimate       lower     upper      p.value
##           <char>      <char>      <num>       <num>     <num>        <num>
## 1: averageSignal     unclear 0.88649462  0.71593356 0.9572006 9.647869e-07
## 2:    averageCor    marginal 0.07516425 -0.16298066 0.3050331 5.381438e-01
## 3:     averageId conditional 0.35941359  0.28714480 0.4276139 0.000000e+00
## 4: averageIdNorm conditional 0.47297607  0.41687698 0.5254856 0.000000e+00
## 5:          lmm0    marginal 0.20151191  0.08427719 0.3132440 8.326941e-04
## 6:          lmm0 conditional 0.35962481  0.27588645 0.4379453 2.442491e-15
## 7:           lmm    marginal 0.07689656 -0.19908544 0.3415776 5.880967e-01
## 8:           lmm conditional 0.46933238  0.39196976 0.5401039 0.000000e+00
plotCor(dtL.scenarioC, type.norm = 2)

dtL1000.scenarioC <- simData(seed = 1, n.obs = 1000,
                             mu.PET = 1:n.region/2,
                             mu.fMRI = 1+5+(1:n.region)/5,
                             rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5)
runCor(PET + fMRI ~ region|id, data = dtL1000.scenarioC)[type!="latent"]
##           method        type  estimate     lower     upper      p.value
##           <char>      <char>     <num>     <num>     <num>        <num>
## 1: averageSignal     unclear 0.9947922 0.9857361 0.9981041 2.676362e-17
## 2:    averageCor    marginal 0.2328200 0.1898082 0.2749400 0.000000e+00
## 3:     averageId conditional 0.3906901 0.3775123 0.4037098 0.000000e+00
## 4: averageIdNorm conditional 0.4931067 0.4818268 0.5042232 0.000000e+00
## 5:          lmm0    marginal 0.2779337 0.2601777 0.2955018 0.000000e+00
## 6:          lmm0 conditional 0.3996423 0.3869348 0.4121980 0.000000e+00
## 7:           lmm    marginal 0.2330568 0.1899371 0.2752791 0.000000e+00
## 8:           lmm conditional 0.5025507 0.4912227 0.5137087 0.000000e+00

## * Figure simulation
figure3 <- plotCor(list("Scenario A" = dtL.scenarioA,"Scenario B" = dtL.scenarioB, "Scenario C" = dtL.scenarioC),
                  type.norm = 2, xspace.dotplot = 0.15, round.dotplot = 1.9, size.dotplot = c(3,2))

pdf("./figures/figure3.pdf", width = 12, height = 8)
figure3
dev.off()

png("./figures/figure3.png", width = 12, height = 8, units = "in", res = 500)
figure3
dev.off()
##----------------------------------------------------------------------
### figure3.R ends here
