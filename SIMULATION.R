### SIMULATION.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2025 (11:43) 
## Version: 
## Last-Updated: feb  4 2025 (17:59) 
##           By: Brice Ozenne
##     Update #: 17
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
##           method        type   estimate      lower     upper      p.value
##           <char>      <char>      <num>      <num>     <num>        <num>
## 1: averageSignal     unclear 0.35523747 -0.1338279 0.7051574 1.480007e-01
## 2:    averageCor    marginal 0.06724627 -0.1762087 0.3029487 5.906723e-01
## 3:     averageId conditional 0.46098596  0.4026315 0.5156086 0.000000e+00
## 4: averageIdNorm conditional 0.46194101  0.4039113 0.5162680 0.000000e+00
## 5:          lmm0    marginal 0.06699296 -0.3278222 0.4418869 6.639662e-01
## 6:          lmm0 conditional 0.45521698  0.3727159 0.5305780 3.729945e-09
## 7:           lmm    marginal 0.06895623 -0.3232252 0.4409231 6.539359e-01
## 8:           lmm conditional 0.45818886  0.3742144 0.5347101 4.781979e-09
plotCor(dtL.scenarioA, type.norm = 2)

dtL1000.scenarioA <- simData(seed = 1, n.obs = 1000,
                             mu.PET = rep(4.8,n.region),
                             mu.fMRI = 1+rep(5,n.region),
                             rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5,
                             sigma.PET = 3, sigma.fMRI = 5)
runCor(PET + fMRI ~ region|id, data = dtL1000.scenarioA)[type!="latent"]
##           method        type  estimate      lower     upper    p.value
##           <char>      <char>     <num>      <num>     <num>      <num>
## 1: averageSignal     unclear 0.5226332 0.07376088 0.7954195 0.02606444
## 2:    averageCor    marginal 0.2341701 0.19123069 0.2762150 0.00000000
## 3:     averageId conditional 0.4929088 0.48146203 0.5041875 0.00000000
## 4: averageIdNorm conditional 0.4920914 0.48065444 0.5033609 0.00000000
## 5:          lmm0    marginal 0.2344417 0.19091624 0.2770472 0.00000000
## 6:          lmm0 conditional 0.5025698 0.49122820 0.5137411 0.00000000
## 7:           lmm    marginal 0.2344234 0.19089460 0.2770321 0.00000000
## 8:           lmm conditional 0.5025487 0.49120107 0.5137259 0.00000000


## ** Scenario B
dtL.scenarioB <- simData(seed = 1,
                         mu.PET = 1:n.region/2,
                         mu.fMRI = 1+5+(1:n.region)/5,
                         rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0, rho.conditional = 0)

runCor(PET + fMRI ~ region|id, data = dtL.scenarioB)[type!="latent"]
##           method        type    estimate      lower      upper      p.value
##           <char>      <char>       <num>      <num>      <num>        <num>
## 1: averageSignal     unclear  0.87715246  0.6948242 0.95353462 1.761851e-06
## 2:    averageCor    marginal -0.19241825 -0.4140914 0.05080116 1.200991e-01
## 3:     averageId conditional  0.27814586  0.2019486 0.35099620 4.508172e-12
## 4: averageIdNorm conditional -0.04499694 -0.1173474 0.02782867 2.258184e-01
## 5:          lmm0    marginal  0.10685976 -0.0383722 0.24767062 1.253738e-01
## 6:          lmm0 conditional  0.27395755  0.1830694 0.36020228 5.260763e-07
## 7:           lmm    marginal -0.19472621 -0.4764036 0.12320251 1.926280e-01
## 8:           lmm conditional -0.06950383 -0.1699130 0.03233363 1.685021e-01
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
## 5:          lmm0    marginal  0.189618749  0.17087978 0.20822059 0.000000e+00
## 6:          lmm0 conditional  0.321129336  0.30757227 0.33455605 0.000000e+00
## 7:           lmm    marginal -0.015554549 -0.05961766 0.02856906 4.884705e-01
## 8:           lmm conditional  0.003255548 -0.01180898 0.01831860 6.715027e-01

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
## 5:          lmm0    marginal 0.20151191  0.06046077 0.3346729 1.189716e-02
## 6:          lmm0 conditional 0.35962481  0.27331419 0.4401911 7.697434e-10
## 7:           lmm    marginal 0.07689656 -0.32366329 0.4540823 6.198342e-01
## 8:           lmm conditional 0.46933238  0.38637452 0.5447577 2.957791e-09
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
## 5:          lmm0    marginal 0.2779337 0.2601057 0.2955722 0.000000e+00
## 6:          lmm0 conditional 0.3996423 0.3869251 0.4122075 0.000000e+00
## 7:           lmm    marginal 0.2330568 0.1894907 0.2757070 0.000000e+00
## 8:           lmm conditional 0.5025507 0.4912031 0.5137278 0.000000e+00

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
