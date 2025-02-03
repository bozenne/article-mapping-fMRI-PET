### ANALYSIS.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2025 (14:06) 
## Version: 
## Last-Updated: jan 29 2025 (18:39) 
##           By: Brice Ozenne
##     Update #: 9
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * load data
df.joint <- read.csv("source/data.csv")



## * Correlation analysis
## ** Strategy 1
df.2cohorts <- aggregate(cbind(pet,asl)~region, data = df.joint,
                         FUN = "mean", na.action = na.pass, na.rm = TRUE)

cor.test(df.2cohorts$pet,df.2cohorts$asl)
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
cor.test(df.joint[df.joint$region=="amygdala","asl"],df.joint[df.joint$region=="amygdala","pet"], use.na = "pairwise")
## 	Pearson's product-moment correlation

## data:  df.joint[df.joint$region == "amygdala", "asl"] and df.joint[df.joint$region == "amygdala", "pet"]
## t = 1.3295, df = 15, p-value = 0.2036
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.1848074  0.6966126
## sample estimates:
##       cor 
## 0.3246774 
cor.testIID(df.joint[df.joint$region=="amygdala","asl"],df.joint[df.joint$region=="amygdala","pet"])

cor.test(x,y)

e.pCor <- partialCor(asl+pet~1, data = df.joint[df.joint$region=="suptempg",])
e.lmmCor <- lmm(CCvalueCC ~ CCvariableCC, repetition = ~CCvariableCC|CCindexCC, structure = "UN", data = attr(e.pCor,"lmm")$data.original, method = "ML")
model.tables(e.lmmCor, effects = "all", transform.rho = FALSE)
##                     estimate          se       df       lower        upper   p.value
## (Intercept)      59.07503085 1.580717075 24.00495  55.8126268  62.33743492 0.0000000
## CCvariableCCpet -57.59673181 1.577691188 24.00478 -60.8528921 -54.34057152 0.0000000
## sigma             7.74390052 1.117735763 23.99051   5.7488577  10.43128891        NA
## k.pet             0.03823814 0.007786652 47.78049   0.0253898   0.05758828 0.0000000
## rho(asl,pet)      0.06913242 0.203148577 12.05997  -0.3585842   0.47285877 0.7402783


cor.testIID(df.joint[df.joint$region=="suptempg","asl"], df.joint[df.joint$region=="suptempg","pet"], transform = FALSE)$se

resCor <- runCor(asl + pet ~ region|cimbi, data = df.joint)
resCor
##            method        type   estimate       lower     upper      p.value
##            <char>      <char>      <num>       <num>     <num>        <num>
##  1: averageSignal conditional 0.29249012 -0.37279449 0.7591641 3.827605e-01
##  2:    averageCor conditional 0.14038232 -0.10331623 0.3681811 2.582675e-01
##  3:     averageId conditional 0.49725395  0.42017948 0.5672053 0.000000e+00
##  4: averageIdNorm conditional 0.03189227 -0.07990364 0.1428958 5.765669e-01
##  5:          lmm0    marginal 0.42675938  0.28279207 0.5519221 3.025870e-04
##  6:          lmm0 conditional 0.51211078  0.42623045 0.5888425 1.052047e-12
##  7:          lmm0      latent 0.09838959 -0.57941824 0.6957443 7.301897e-01
##  8:           lmm    marginal 0.13819072 -0.15605055 0.4099052 2.974110e-01
##  9:           lmm conditional 0.07672043 -0.03635723 0.1878584 1.731352e-01
## 10:           lmm      latent 0.18933946 -0.34464758 0.6307502 4.183106e-01




##----------------------------------------------------------------------
### ANALYSIS.R ends here
