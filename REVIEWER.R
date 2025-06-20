### REVIEWER.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 20 2025 (11:01) 
## Version: 
## Last-Updated: jun 20 2025 (17:28) 
##           By: Brice Ozenne
##     Update #: 22
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * load code
library(ggpubr)
library(ggplot2)
library(pbapply)
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


## * Single cohort
strategyAll.cor <- runCor(asl + pet ~ region|cimbi, data = df.joint)[type!="latent"]
strategyAll.cor
##           method        type   estimate       lower     upper      p.value
##           <char>      <char>      <num>       <num>     <num>        <num>
## 1: averageSignal     unclear 0.67160959  0.29826037 0.8667168 2.271529e-03
## 2:    averageCor    marginal 0.14038232 -0.14331502 0.4027605 3.321895e-01
## 3:     averageId conditional 0.49725395  0.42017948 0.5672053 0.000000e+00
## 4: averageIdNorm conditional 0.03189227 -0.07990364 0.1428958 5.765669e-01
## 5:          lmm0    marginal 0.42675938  0.30797055 0.5324466 8.415757e-11
## 6:          lmm0 conditional 0.51211078  0.42906508 0.5865716 0.000000e+00
## 7:           lmm    marginal 0.13819072 -0.10048466 0.3618267 2.558514e-01
## 8:           lmm conditional 0.07672043 -0.03011947 0.1818271 1.591055e-01

system.time(
    strategy3.corlmmH <- partialCor(asl+pet~region, repetition = ~region|cimbi, data = df.joint, structure = "HCS", df = FALSE)
)

## * Initialize cluster
cl <- parallel::makeCluster(6)
parallel::clusterExport(cl, varlist = c("df.joint"))
tempo <- parallel::clusterCall(cl = cl, function(x){
    library(data.table)
    library(LMMstar)
    source("cor.testIID.R")
    source("FCT.R")
})

## * Split cohort

## ** illustration
set.seed(5)
sort(sample.int(24, size = 6, replace = FALSE))

df.joint1 <- df.joint[df.joint$cimbi %in% unique(df.joint$cimbi)[c(2,7,9,11,15,23)],]
df.joint2 <- df.joint[df.joint$cimbi %in% unique(df.joint$cimbi)[c(1,3:6,8)],]

dfmean.cohort12 <- merge(aggregate(pet~region, data = df.joint1, FUN = "mean", na.action = na.pass, na.rm = TRUE),
                         aggregate(asl~region, data = df.joint2, FUN = "mean", na.action = na.pass, na.rm = TRUE),
                         by="region")
cor.test(dfmean.cohort12$pet, dfmean.cohort12$asl)
## 	Pearson's product-moment correlation

## data:  dfmean.cohort12$pet and dfmean.cohort12$asl
## t = 2.3007, df = 16, p-value = 0.0352
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.04133887 0.78315698
## sample estimates:
##       cor 
## 0.4985863 

dfmean.cohort21 <- merge(aggregate(pet~region, data = df.joint2, FUN = "mean", na.action = na.pass, na.rm = TRUE),
                         aggregate(asl~region, data = df.joint1, FUN = "mean", na.action = na.pass, na.rm = TRUE),
                         by="region")
cor.test(dfmean.cohort21$pet, dfmean.cohort21$asl)
## 	Pearson's product-moment correlation

## data:  dfmean.cohort21$pet and dfmean.cohort21$asl
## t = 3.409, df = 16, p-value = 0.003591
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.2607323 0.8562222
## sample estimates:
##       cor 
## 0.6486441 

dfmean.cohort1 <- aggregate(cbind(pet,asl)~region, data = df.joint1,
                            FUN = "mean", na.action = na.pass, na.rm = TRUE)
cor.test(dfmean.cohort1$pet, dfmean.cohort1$asl)
## 	Pearson's product-moment correlation

## data:  dfmean.cohort1$pet and dfmean.cohort1$asl
## t = 3.8552, df = 16, p-value = 0.0014
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.3359229 0.8767661
## sample estimates:
##       cor 
## 0.6939579 

dfmean.cohort2 <- aggregate(cbind(pet,asl)~region, data = df.joint2,
                            FUN = "mean", na.action = na.pass, na.rm = TRUE)
cor.test(dfmean.cohort2$pet, dfmean.cohort2$asl)
## 	Pearson's product-moment correlation

## data:  dfmean.cohort2$pet and dfmean.cohort2$asl
## t = 1.7801, df = 16, p-value = 0.09405
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.07441109  0.73410479
## sample estimates:
##       cor 
## 0.4065838 

## ** simulation
n.split <- 1000

df.split <- do.call(rbind, pblapply(1:n.split, function(iSeed){ ## iSeed <- 1
    iOut <- NULL
    for(iSize in 4:12){ ## iSize <- 6
        iData <- randomSplitHalf(df.joint, cluster = "cimbi", size = c(iSize,iSize), seed = iSeed)
        iDataMean.same <- aggregate(cbind(pet,asl)~region, data = iData[iData$cohort == 1, ],
                                    FUN = "mean", na.action = na.pass, na.rm = TRUE)
        iCor.same <- cor.test(x = iDataMean.same$pet, y = iDataMean.same$asl)
        iDataMean.diff <- aggregate(cbind(cohort.pet,cohort.asl)~region, data = iData,
                                    FUN = "mean", na.action = na.pass, na.rm = TRUE)
        iCor.diff <- cor.test(x = iDataMean.diff$cohort.pet, y = iDataMean.diff$cohort.asl)
        iOut <- rbind(iOut,
                      data.frame(size = iSize, cohort = "same", estimate = iCor.same$estimate, lower = iCor.same$conf.int[1], upper = iCor.same$conf.int[2], p.value = iCor.same$p.value),
                      data.frame(size = iSize, cohort = "diff", estimate = iCor.diff$estimate, lower = iCor.diff$conf.int[1], upper = iCor.diff$conf.int[2], p.value = iCor.diff$p.value),
                      data.frame(size = iSize, cohort = "contrast", estimate = iCor.same$estimate - iCor.diff$estimate, lower = NA, upper = NA, p.value = NA)
                      )        
    }
    return(cbind(seed = iSeed, iOut))
}, cl = cl)) ## elapsed=20s  

dfS.split <- summarize(estimate ~ cohort + size, data = df.split)
print(dfS.split, digits =  2)
##      cohort size observed missing     mean    sd   min       q1   median    q3  max
## 1  contrast    4     1000       0  0.13583 0.131 -0.34  0.04925  0.15479 0.227 0.46
## 2      diff          1000       0  0.48905 0.054  0.30  0.46455  0.50234 0.502 0.59
## 3      same          1000       0  0.62488 0.106  0.23  0.56121  0.63500 0.707 0.85
## 4  contrast    5     1000       0  0.14908 0.122 -0.30  0.07596  0.16916 0.235 0.41
## 5      diff          1000       0  0.48591 0.059  0.23  0.45857  0.48664 0.519 0.64
## 6      same          1000       0  0.63500 0.094  0.31  0.57413  0.64563 0.702 0.83
## 7  contrast    6     1000       0  0.15128 0.128 -0.32  0.07236  0.16634 0.247 0.44
## 8      diff          1000       0  0.49036 0.068  0.28  0.44453  0.49809 0.534 0.65
## 9      same          1000       0  0.64165 0.087  0.31  0.58320  0.64778 0.704 0.84
## 10 contrast    7     1000       0  0.13779 0.122 -0.25  0.05532  0.14376 0.227 0.43
## 11     diff          1000       0  0.50709 0.065  0.27  0.46828  0.50693 0.551 0.66
## 12     same          1000       0  0.64488 0.079  0.39  0.59485  0.64821 0.704 0.83
## 13 contrast    8     1000       0  0.12941 0.118 -0.30  0.05576  0.13618 0.216 0.41
## 14     diff          1000       0  0.52062 0.060  0.33  0.48063  0.52238 0.561 0.69
## 15     same          1000       0  0.65003 0.072  0.40  0.60392  0.65486 0.704 0.82
## 16 contrast    9     1000       0  0.08855 0.129 -0.32  0.00087  0.09614 0.179 0.41
## 17     diff          1000       0  0.56342 0.069  0.35  0.51603  0.56290 0.611 0.74
## 18     same          1000       0  0.65197 0.067  0.41  0.60564  0.65755 0.701 0.81
## 19 contrast   10     1000       0  0.03786 0.120 -0.28 -0.05090  0.04317 0.122 0.37
## 20     diff          1000       0  0.61739 0.065  0.43  0.57561  0.61869 0.665 0.77
## 21     same          1000       0  0.65525 0.061  0.48  0.61320  0.65876 0.698 0.82
## 22 contrast   11     1000       0  0.02012 0.109 -0.27 -0.05535  0.02296 0.098 0.32
## 23     diff          1000       0  0.63813 0.055  0.47  0.60035  0.64075 0.677 0.78
## 24     same          1000       0  0.65825 0.057  0.49  0.62097  0.66166 0.699 0.80
## 25 contrast   12     1000       0 -0.00052 0.102 -0.26 -0.07096 -0.00072 0.071 0.28
## 26     diff          1000       0  0.66136 0.052  0.50  0.62685  0.66321 0.698 0.78
## 27     same          1000       0  0.66084 0.052  0.52  0.62561  0.66326 0.697 0.78

df.split$cohort.label <- factor(df.split$cohort, levels = c("diff","same","contrast"), labels = c("Different subjects","Same subjects","Same - Different"))
gg.split <- ggplot(df.split[df.split$cohort %in% c("diff","same"),], aes(x = as.factor(size), y = estimate, fill = cohort.label))
gg.split <- gg.split + geom_boxplot()
gg.split <- gg.split + labs(x = "Number of available values per modality", y = "Pearson \u03c1", fill = "Availability of CBF vs. PET values")
gg.split

gg.split2 <- ggplot(df.split[df.split$cohort %in% c("contrast"),], aes(x = as.factor(size), y = estimate))
gg.split2 <- gg.split2 + geom_boxplot()
gg.split2 <- gg.split2 + labs(x = "Number of available values per modality", y = "Difference in Pearson \u03c1")
gg.split2

gg.split2abs <- ggplot(df.split[df.split$cohort %in% c("contrast"),], aes(x = as.factor(size), y = abs(estimate)))
gg.split2abs <- gg.split2abs + geom_boxplot()
gg.split2abs <- gg.split2abs + labs(x = "Number of available values per modality", y = "Absolute difference in Pearson \u03c1")
gg.split2abs

gg.splitAll <- ggarrange(gg.split + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2, byrow = TRUE)), gg.split2, gg.split2abs, nrow = 1)


## ** export
ggsave(gg.splitAll, filename = file.path("figures","figureReviewer.png"),
       width = 12, height = 5)


intersect(df.split[df.split$size == 6 & df.split$cohort == "diff" & df.split$estimate < 0.45,"seed"], df.split[df.split$size == 6 & df.split$cohort == "same" & df.split$estimate > 0.6,"seed"])

##----------------------------------------------------------------------
### REVIEWER.R ends here
