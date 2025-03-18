### ANALYSIS.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2025 (14:06) 
## Version: 
## Last-Updated: mar 18 2025 (13:28) 
##           By: Brice Ozenne
##     Update #: 69
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * load code
library(data.table)
library(LMMstar)
source("cor.testIID.R")
source("FCT.R")

## * load data
df.joint <- read.csv("source/data.csv")
region.name <- c('DLPFC' = 'Dorsolateral PFC', 
                 'VLPFC' = 'Ventrolateral PFC', 
                 'amygdala' = 'Amygdala',
                 'antcin' = 'Anterior Cingulate Cortex', 
                 'cau' = 'Caudate',                      
                 'entorhinc' = 'Entorhinal Cortex',            
                 'hippocampus' = 'Hippocampus',                  
                 'ins' = 'Insula',                       
                 'medinffg' = 'Medial-Inferior Frontal Gyrus',
                 'occ' = 'Occipital Cortex',             
                 'orbfrc' = 'Orbitofrontal Cortex',         
                 'parc' = 'Parietal Cortex',              
                 'postcin' = 'Posterior Cingulate Cortex',   
                 'put' = 'Putamen',                      
                 'senmotc' = 'Sensory-Motor Cortex',         
                 'supfg' = 'Superior Frontal Gyrus',       
                 'suptempg' = 'Superior Temporal Gyrus',      
                 'th' = 'Thalamus')                     
all(df.joint$region %in% names(region.name))
df.joint$region.full <- region.name[df.joint$region]

## * Correlation analysis
## ** All at once
strategyAll.cor <- runCor(asl + pet ~ region|cimbi, data = df.joint)[type!="latent"]
strategyAll.cor
##            method        type   estimate       lower     upper      p.value
##           <char>      <char>      <num>       <num>     <num>        <num>
## 1: averageSignal     unclear 0.67160959  0.29826037 0.8667168 2.271529e-03
## 2:    averageCor    marginal 0.14038232 -0.14331502 0.4027605 3.321895e-01
## 3:     averageId conditional 0.49725395  0.42017948 0.5672053 0.000000e+00
## 4: averageIdNorm conditional 0.03189227 -0.07990364 0.1428958 5.765669e-01
## 5:          lmm0    marginal 0.42675938  0.30797055 0.5324466 8.415757e-11
## 6:          lmm0 conditional 0.51211078  0.42906508 0.5865716 0.000000e+00
## 7:           lmm    marginal 0.13819072 -0.10048466 0.3618267 2.558514e-01
## 8:           lmm conditional 0.07672043 -0.03011947 0.1818271 1.591055e-01

## ** Strategy 1
df.2cohorts <- aggregate(cbind(pet,asl)~region, data = df.joint,
                         FUN = "mean", na.action = na.pass, na.rm = TRUE)

strategy1.cor <- cor.test(df.2cohorts$pet,df.2cohorts$asl)
strategy1.cor
## 	Pearson's product-moment correlation

## data:  df.2cohorts$pet and df.2cohorts$asl
## t = 3.6259, df = 16, p-value = 0.002272
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.2982604 0.8667168
## sample estimates:
##       cor 
## 0.6716096 

## ** Strategy 2
ls.cor <- by(df.joint, INDICES = df.joint$region, FUN = function(data){cor.test(x = data$pet, y = data$asl, use.na = "pairwise")})
strategy2.cor0 <- cbind(do.call(rbind,lapply(ls.cor, "[[","estimate")),do.call(rbind,lapply(ls.cor, "[[","conf.int")), do.call(rbind,lapply(ls.cor, "[[","p.value")))
colnames(strategy2.cor0) <- c("estimate","lower","upper","p.value")
strategy2.cor <- rbind(strategy2.cor0,
                       "average" = as.data.frame(strategyAll.cor[strategyAll.cor$method=="averageCor",.(estimate,lower,upper,p.value)]))
## Model based se
##                estimate       lower     upper   p.value
## amygdala     0.32467739 -0.18480740 0.6966126 0.2035508
## antcin       0.33620879 -0.07773011 0.6512759 0.1082047
## cau         -0.11204262 -0.50957949 0.3249142 0.6195993
## DLPFC        0.21740643 -0.20386987 0.5707480 0.3074987
## entorhinc    0.29738063 -0.33342585 0.7442612 0.3478742
## hippocampus  0.03879789 -0.37947304 0.4439009 0.8604845
## ins          0.16695887 -0.25351205 0.5343657 0.4355316
## medinffg     0.18215496 -0.23878806 0.5454701 0.3942720
## occ          0.24009490 -0.19101120 0.5935544 0.2698181
## orbfrc       0.12029423 -0.29754104 0.4994543 0.5755485
## parc         0.15988550 -0.26029984 0.5291533 0.4554967
## postcin      0.15508873 -0.26487945 0.5256025 0.4693034
## put          0.12808294 -0.36033444 0.5613834 0.6125110
## senmotc      0.01200828 -0.39329359 0.4134024 0.9555891
## supfg        0.18580245 -0.23522485 0.5481167 0.3847075
## suptempg     0.06913242 -0.34385361 0.4597089 0.7482235
## th          -0.05540597 -0.45716715 0.3651353 0.8017429
## VLPFC        0.06035589 -0.35160139 0.4527298 0.7793602
## average      0.14038232 -0.14331502 0.4027605 0.3321895

## sanity check
mean(strategy2.cor0[,"estimate"])
## [1] 0.1403823
range(strategy2.cor[,"estimate"])
## [1] -0.1120426  0.3362088
min(strategy2.cor[,"p.value"])
## 0.1082047

## number of regions
NROW(strategy2.cor0)
## [1] 18

## ** Strategy 3
strategy3.corlmm <- partialCor(asl+pet~region, repetition = ~region|cimbi, data = df.joint, structure = "CS", df = FALSE)

## split output
strategy3.cor <- strategy3.corlmm
attr(strategy3.cor,"lmm") <- NULL
strategy3.lmm <- attr(strategy3.corlmm,"lmm")

strategy3.cor
##             estimate     se    df   lower upper p.value
## marginal      0.1382 0.1201 Inf -0.1005 0.362   0.256
## conditional   0.0767 0.0543 Inf -0.0301 0.182   0.159
## latent        0.1893 0.2108 Inf -0.2325 0.551   0.381

## partialCor(asl+pet~region, repetition = ~region|cimbi, data = df.joint, structure = "CS", df = TRUE)
##             estimate     se    df   lower upper p.value
## marginal      0.1382 0.1201  6.27 -0.1561 0.410   0.297
## conditional   0.0767 0.0543 21.91 -0.0364 0.188   0.173
## latent        0.1893 0.2108  5.36 -0.3446 0.631   0.418

model.tables(strategy3.lmm, effects = c("variance","correlation"))
##                 estimate          se        df       lower       upper      p.value
## sigma         9.27035474 0.716169306  4.895445  7.59091508 11.32135929           NA
## k.pet         0.02746754 0.003372901 11.107797  0.02096915  0.03597978 7.220446e-12
## rho(1,2,dt=0) 0.13819072 0.120065540  6.269444 -0.15605055  0.40990521 2.974110e-01
## rho1          0.48151823 0.080286274  5.702418  0.25983972  0.65498593 2.776387e-03
## rho(1,2,dt=1) 0.10496416 0.119125262  5.726731 -0.19047218  0.38296658 4.169382e-01
## rho2          0.63824446 0.072511481  8.241502  0.44182255  0.77627125 2.366809e-04

## * Extra analyses

## ** Strategy 3 (simplified)
dt.joint <- as.data.table(df.joint)
dt.joint[,asl.norm := scale(asl), by = "region"]
dt.joint[,pet.norm := scale(pet), by = "region"]

dtS.joint <- dt.joint[, .(cor = cor(asl.norm, pet.norm, use = "pairwise")), by = "cimbi"]
t.test(dtS.joint$cor)
## 	One Sample t-test

## data:  dtS.joint$cor
## t = 0.55879, df = 23, p-value = 0.5817
## alternative hypothesis: true mean is not equal to 0
## 95 percent confidence interval:
##  -0.08617488  0.14995943
## sample estimates:
##  mean of x 
## 0.03189227 

##  4: averageIdNorm conditional 0.03189227 -0.07990364 0.1428958 5.765669e-01

## ** Strategy 3 without assuming same variance
## WARNING: time consuming ~2 min
system.time(
    strategy3.corlmmH <- partialCor(asl+pet~region, repetition = ~region|cimbi, data = df.joint, structure = "HCS", df = FALSE)
)

## split output
strategy3.corH <- strategy3.corlmmH
attr(strategy3.corH,"lmm") <- NULL
strategy3.lmmH <- attr(strategy3.corlmmH,"lmm")

strategy3.corH
##             estimate     se  df  lower upper p.value
## marginal      0.1353 0.1314 Inf -0.125 0.379   0.309
## conditional   0.0574 0.0557 Inf -0.052 0.165   0.304
## latent        0.1847 0.2090 Inf -0.233 0.545   0.388

## * Gather summary statistics for parametrizing simulations
gridASL.region <- data.frame(CCvariableCC = "asl", region = sort(names(region.name)))
gridPET.region <- data.frame(CCvariableCC = "pet", region = sort(names(region.name)))
     
param.lmm <- list(mu.asl = setNames(predict(strategy3.lmm, newdata = gridASL.region), sort(names(region.name))),
                     mu.pet = setNames(predict(strategy3.lmm, newdata = gridPET.region), sort(names(region.name))),
                     sigma.asl = coef(strategy3.lmm, effects = "variance", transform.k = "sd")["sigma.asl"],
                     sigma.pet = coef(strategy3.lmm, effects = "variance", transform.k = "sd")["sigma.pet"],
                     rho.asl = coef(strategy3.lmm, effects = "correlation")["rho1"],
                     rho.pet = coef(strategy3.lmm, effects = "correlation")["rho2"],
                     rho.petasl = coef(strategy3.lmm, effects = "correlation")["rho(1,2,dt=0)"],
                     rhoLag.petasl = coef(strategy3.lmm, effects = "correlation")["rho(1,2,dt=1)"]
                     )

allSigma <- coef(strategy3.lmmH, effects = "variance", transform.k = "sd")
param.lmmH <- list(mu.asl = setNames(predict(strategy3.lmmH, newdata = gridASL.region), sort(names(region.name))),
                     mu.pet = setNames(predict(strategy3.lmmH, newdata = gridPET.region), sort(names(region.name))),
                     sigma.asl = setNames(allSigma[grepl("asl",names(allSigma))], gsub("sigma.asl:","",grep("asl",names(allSigma),value = TRUE),fixed = TRUE)),
                     sigma.pet = setNames(allSigma[grepl("pet",names(allSigma))], gsub("sigma.pet:","",grep("pet",names(allSigma),value = TRUE),fixed = TRUE)),
                     rho.asl = coef(strategy3.lmmH, effects = "correlation")["rho1"],
                     rho.pet = coef(strategy3.lmmH, effects = "correlation")["rho2"],
                     rho.petasl = coef(strategy3.lmmH, effects = "correlation")["rho(1,2,dt=0)"],
                     rhoLag.petasl = coef(strategy3.lmmH, effects = "correlation")["rho(1,2,dt=1)"]
                     )


## * Export
saveRDS(region.name, file = "results/regionName.rds")
saveRDS(df.2cohorts, file = "results/df2cohorts.rds")

saveRDS(strategy1.cor, file = "results/strategy1.rds")
saveRDS(strategy2.cor, file = "results/strategy2.rds")
saveRDS(strategy3.cor, file = "results/strategy3.rds")
saveRDS(strategy3.corH, file = "results/strategy3H.rds")

saveRDS(param.lmm, file = "results/strategy3-lmm-param.rds")
saveRDS(param.lmmH, file = "results/strategy3-lmmH-param.rds")


saveRDS(strategy3.lmm, file = "source/strategy3-lmm.rds") ## not anonymized and heavy
saveRDS(strategy3.lmmH, file = "source/strategy3-lmmH.rds") ## not anonymized and heavy


 
##----------------------------------------------------------------------
### analysis.R ends here
