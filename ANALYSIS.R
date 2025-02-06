### ANALYSIS.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2025 (14:06) 
## Version: 
## Last-Updated: feb  5 2025 (17:48) 
##           By: Brice Ozenne
##     Update #: 42
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
library(ggplot2)
library(ggrepel)
library(LMMstar)
library(officer)
library(flextable)
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
strategyAll.cor <- runCor(asl + pet ~ region|cimbi, data = df.joint)
strategyAll.cor
##            method        type   estimate       lower     upper      p.value
##            <char>      <char>      <num>       <num>     <num>        <num>
##  1: averageSignal     unclear 0.67160959  0.29826037 0.8667168 2.271529e-03
##  2:    averageCor    marginal 0.14038232 -0.14331502 0.4027605 3.321895e-01
##  3:     averageId conditional 0.49725395  0.42017948 0.5672053 0.000000e+00
##  4: averageIdNorm conditional 0.03189227 -0.07990364 0.1428958 5.765669e-01
##  5:          lmm0    marginal 0.42675938  0.28279207 0.5519221 3.025870e-04
##  6:          lmm0 conditional 0.51211078  0.42623045 0.5888425 1.052047e-12
##  7:          lmm0      latent 0.09838959 -0.57941824 0.6957443 7.301897e-01
##  8:           lmm    marginal 0.13819072 -0.15605055 0.4099052 2.974110e-01
##  9:           lmm conditional 0.07672043 -0.03635723 0.1878584 1.731352e-01
## 10:           lmm      latent 0.18933946 -0.34464758 0.6307502 4.183106e-01

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

## ** Strategy 3
strategy3.cor <- partialCor(asl+pet~region, repetition = ~region|cimbi, data = df.joint, structure = "CS")
strategy3.cor
##             estimate     se    df   lower upper p.value
## marginal      0.1382 0.1201  6.27 -0.1561 0.410   0.297
## conditional   0.0767 0.0543 21.91 -0.0364 0.188   0.173
## latent        0.1893 0.2108  5.36 -0.3446 0.631   0.418

model.tables(attr(strategy3.cor,"lmm"), effects = c("variance","correlation"))
##                 estimate          se        df       lower       upper      p.value
## sigma         9.27035474 0.716169306  4.895445  7.59091508 11.32135929           NA
## k.pet         0.02746754 0.003372901 11.107797  0.02096915  0.03597978 7.220446e-12
## rho(1,2,dt=0) 0.13819072 0.120065540  6.269444 -0.15605055  0.40990521 2.974110e-01
## rho1          0.48151823 0.080286274  5.702418  0.25983972  0.65498593 2.776387e-03
## rho(1,2,dt=1) 0.10496416 0.119125262  5.726731 -0.19047218  0.38296658 4.169382e-01
## rho2          0.63824446 0.072511481  8.241502  0.44182255  0.77627125 2.366809e-04

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

## * Figures
## ** strategy 1
nsig <- 2
sub.text.corr <- paste0("Pearson's rho: ", signif(strategy1.cor$estimate,nsig),
                        ', 95% CI: ', paste(signif(strategy1.cor$conf.int,nsig), collapse = '; '),
                        ', p = ', signif(strategy1.cor$p.value,nsig))

gg.strategy1 <- ggplot(cbind(df.2cohorts, region.full = region.name[df.2cohorts$region]), aes(x = pet, y = asl)) 
gg.strategy1 <- gg.strategy1 + geom_smooth(method = 'lm', alpha = .8, color = 'black') 
gg.strategy1 <- gg.strategy1 + geom_point(aes(color = region)) 
gg.strategy1 <- gg.strategy1 + theme_classic() 
gg.strategy1 <- gg.strategy1 + theme(axis.text = element_text(size=20),
                                     axis.title = element_text(size=24,face="bold"), 
                                     legend.position = 'none',
                                     plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
                                     plot.subtitle = element_text(hjust = 0.5, size = 20))
gg.strategy1 <- gg.strategy1 + labs(x = '[11C]Cimbi-36 BPND', y = 'Baseline CBF (ml/100g/min)', title = '5-HT2AR vs. CBF', subtitle = sub.text.corr)
gg.strategy1 <- gg.strategy1 + scale_x_continuous(breaks = seq(0.2,1.6,0.2))
gg.strategy1 <- gg.strategy1 + scale_y_continuous(limits = c(35,70), breaks = seq(35,70,5))
gg.strategy1 <- gg.strategy1 + geom_label_repel(aes(label = region.full, 
                                                    fill = region), size = 5)
gg.strategy1
    
if(FALSE){
    ggsave(gg.strategy1, filename = file.path("figures","gg-analysis-strategy1.pdf"), units = "in", width = 12.5, height = 7.5)
    ggsave(gg.strategy1, filename = file.path("figures","gg-analysis-strategy1.png"), units = "in", width = 12.5, height = 7.5, dpi=300)
}

## ** strategy 2
dfGG.strategy2 <- cbind(region = rownames(strategy2.cor),
                        region.full = c("average" = "Average", region.name)[rownames(strategy2.cor)],
                        strategy2.cor)
dfGG.strategy2$region.full <- factor(dfGG.strategy2$region.full, levels = unique(dfGG.strategy2$region.full))
dfGG.strategy2$style <- dfGG.strategy2$region=="average"


gg.strategy2 <- ggplot(dfGG.strategy2, aes(x = region.full, y =  estimate, ymin = lower, ymax = upper) )
gg.strategy2 <- gg.strategy2 + geom_point(aes(shape = style, size = style)) + geom_linerange(aes(linewidth = style))
gg.strategy2 <- gg.strategy2 + geom_hline(yintercept = 0, linetype = "dashed")
gg.strategy2 <- gg.strategy2 + coord_cartesian(ylim = c(-1,1))
gg.strategy2 <- gg.strategy2 + labs(x="",y="Pearson's correlation coefficient") + guides(linewidth = "none", shape = "none", size = "none")
gg.strategy2 <- gg.strategy2 + scale_size_manual(breaks = c(FALSE,TRUE), values = c(2,3))
gg.strategy2 <- gg.strategy2 + scale_shape_manual(breaks = c(FALSE,TRUE), values = c(19,15))
gg.strategy2 <- gg.strategy2 + scale_linewidth_manual(breaks = c(FALSE,TRUE), values = c(1,1.5))
gg.strategy2 <- gg.strategy2 + theme(text = element_text(size=15), ,
                                     axis.text.x=element_text(angle=90, size = 8),                                     
                                     axis.line = element_line(linewidth = 1.25),
                                     axis.ticks = element_line(linewidth = 2),
                                     axis.ticks.length=unit(.25, "cm"))
gg.strategy2 <- gg.strategy2 + scale_x_discrete(breaks = levels(dfGG.strategy2$region.full),
                                                labels = c(levels(dfGG.strategy2$region.full)[-length(dfGG.strategy2$region.full)],
                                                           expression(bold("average"))))
gg.strategy2

if(FALSE){
    ggsave(gg.strategy2, filename = file.path("figures","gg-analysis-strategy2.pdf"), width = 8, height = 5)
    ggsave(gg.strategy2, filename = file.path("figures","gg-analysis-strategy2.png"), width = 8, height = 5)
}

## * Tables
## ** strategy 2
n.digit <- 3

df.table1 <- cbind(" "=c("average"="Average",region.name)[rownames(strategy2.cor)],
                   "Estimated correlation" = formatC(round(strategy2.cor$estimate,n.digit), format = "f", digits = n.digit),
                   "Confidence interval" = paste0("[",formatC(round(strategy2.cor$lower,n.digit), format = "f", digits = n.digit),
                                                  "; ",formatC(round(strategy2.cor$upper,n.digit), format = "f", digits = n.digit),"]"),
                   "P-value" = formatC(round(strategy2.cor$p.value,n.digit), format = "f", digits = n.digit))
## formatC keeps trailing 0

table1 <- read_docx()
table1 <- body_add_flextable(table1,set_table_properties(autofit(flextable(as.data.frame(df.table1))), width = 0.9, layout = "autofit"))
print(table1, target = "tables/table1.docx")

##----------------------------------------------------------------------
### analysis.R ends here
