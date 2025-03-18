### table2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 18 2025 (11:03) 
## Version: 
## Last-Updated: mar 18 2025 (16:39) 
##           By: Brice Ozenne
##     Update #: 13
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(data.table)
library(officer)
library(flextable)
format.pval2 <- function(pv, ...){
    ## WARNING: no difference between 0<x<0.001 and -0.001<x<0 (both are coded <0.001)
    sign <- sign(pv)
    out <- format.pval(pv = abs(pv), ...)
    out[sign<0 & !grepl("<",out)] <- paste0("-",out[sign<0 & !grepl("<",out)])
    return(out)
}

## * load data
dtS.sim <-  readRDS(file = "results/simulation-summary.rds")

dtL.table2 <- dtS.sim[scenario %in% c("A","B","C") & method %in% c("averageSignal","averageCor","lmm") & type %in% c("unclear","marginal")]
dtL.table2[, strategy := factor(method, levels = c("averageSignal","averageCor","lmm"), labels = 1:3)]

## * prepare table

## ** process simulation results
dtW.table2 <- dcast(dtL.table2, scenario+strategy~n, value.var = c("bias","std","coverage","rejection.rate"))
## options(width = 140)
dt.table2 <- cbind(dtW.table2[,.(scenario,strategy)],
                   Bias = paste0(format.pval2(dtW.table2$bias_26, digits = 1, eps = 1e-3), "/",format.pval2(dtW.table2$bias_1000, digits = 1, eps = 1e-3)),
                   Std = paste0(formatC(dtW.table2$std_26, digits = 3, format = "f"), "/",formatC(dtW.table2$std_1000, digits = 3, format = "f")),
                   Coverage = paste0(formatC(dtW.table2$coverage_26, digits = 3, format = "f"), "/",formatC(dtW.table2$coverage_1000, digits = 3, format = "f")),
                   "Type 1 error" = paste0(formatC(dtW.table2$rejection.rate_26, digits = 3, format = "f"), "/",formatC(dtW.table2$rejection.rate_1000, digits = 3, format = "f")),
                   Power = paste0(formatC(dtW.table2$rejection.rate_26, digits = 3, format = "f"), "/",formatC(dtW.table2$rejection.rate_1000, digits = 3, format = "f"))
)
dt.table2[scenario != "B", "Type 1 error" := ""]
dt.table2[scenario == "B", "Power" := ""]
names(dt.table2)[3:7] <- paste(names(dt.table2)[3:7], "(n=26/n=1000)")
## WARNING: no difference between 0<x<0.001 and -0.001<x<0 (both are coded <0.001)

##    scenario strategy Bias (n=26/n=1000) Std (n=26/n=1000) Coverage (n=26/n=1000) Type 1 error (n=26/n=1000) Power (n=26/n=1000)
##      <char>   <fctr>             <char>            <char>                 <char>                     <char>              <char>
## 1:        A        1      -0.008/-0.009       0.187/0.187            0.952/0.953                                    0.597/0.593
## 2:        A        2      -0.005/<0.001       0.116/0.018            0.939/0.949                                    0.522/1.000
## 3:        A        3       0.002/<0.001       0.115/0.018            0.940/0.951                                    0.559/1.000
## 4:        B        1        0.664/0.674       0.033/0.005            0.000/0.000                1.000/1.000                    
## 5:        B        2      -0.001/<0.001       0.116/0.018            0.943/0.952                0.057/0.048                    
## 6:        B        3      -0.001/<0.001       0.115/0.018            0.945/0.952                0.055/0.048                    
## 7:        C        1        0.169/0.174       0.028/0.005            1.000/1.000                                    1.000/1.000
## 8:        C        2      -0.005/<0.001       0.116/0.018            0.939/0.949                                    0.522/1.000
## 9:        C        3       0.002/<0.001       0.115/0.018            0.940/0.951                                    0.559/1.000

## NOTE 1: estimand 1 corresponds to strategy 1
##         estimand 2 (marginal) corresponds to strategy 2 and 3 
##         estimand 3 (conditional) is not described in the article
##
## NOTE 2: scenario D is an extra simulation not described in the article with non-homogeneous variance across regions.
##         the LMM ignores this issue and is fitted as if all the regions have the same variance within modality.
##         Using partialCor(asl+pet~region, repetition = ~region|cimbi, data = df.joint, structure = "HCS", df = FALSE) does not make this restriction
##         but is not investigated in the simulation study (partly because it is much slower, partly because it is beyond the scope of the paper).


## ** format table for article

dt.table2 <- table.sim[scenario %in% c("A","B","C") & trimws(estimand) %in% c("1","2","2 (LMM)")]
dt.table2$estimand <- droplevels(dt.table2$estimand) ## drop estimand 3 levels
dt.table2$estimand <- factor(dt.table2$estimand, labels = trimws(levels(dt.table2$estimand))) ## remove white space
dt.table2$estimand <- factor(dt.table2$estimand, levels = c("1","2","2 (LMM)"), labels = 1:3) ## rename
names(dt.table2)[names(dt.table2)=="estimand"] <- "strategy"

table2 <- read_docx()
table2 <- body_add_flextable(table2,set_table_properties(autofit(flextable(as.data.frame(dt.table2))), width = 0.9, layout = "autofit"))

## * Result as text
truth.marginal <- dt.sim[type == "marginal", .(truth = unique(truth)), by = "scenario"]


merge(dtSW.sim, truth.marginal, by = "scenario")[scenario %in% c("A","C") & trimws(estimand) %in% c("1","2","2 (LMM)"),
                                                 .("bias (%)" = 100*bias_26/truth),
                                                 by = c("estimand","scenario")]
##    estimand scenario   bias (%)
##      <fctr>   <char>      <num>
## 1:  1              A   2.132179
## 2:  2              A   1.609760
## 3:  2 (LMM)        A  -1.274466
## 4:  1              C -67.181214
## 5:  2              C   1.609760
## 6:  2 (LMM)        C  -1.274466


## * export results 
print(table2, target = "tables/table2.docx")

##----------------------------------------------------------------------
### table2.R ends here
