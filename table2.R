### table2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 18 2025 (11:03) 
## Version: 
## Last-Updated: mar 21 2025 (17:03) 
##           By: Brice Ozenne
##     Update #: 27
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
format.pval2 <- function(pv, digits, eps){
    ## WARNING: no difference between 0<x<0.001 and -0.001<x<0 (both are coded <0.001)
    sign <- sign(pv)
    small <- abs(pv)<eps
    out <- trimws(formatC(pv, digits = digits, format = "f"))
    out[is.infinite(pv)] <- NA
    out[sign>0 & small] <- paste0("<",eps)
    out[sign<0 & small] <- paste0("<-",eps)
    return(out)
}

## * load data
dtS.sim <-  readRDS(file = "results/simulation-summary.rds")

dtL.table2 <- dtS.sim[method %in% c("averageSignal","averageCor","lmm") & type %in% c("unclear","marginal")]
dtL.table2[, strategy := factor(method, levels = c("averageSignal","averageCor","lmm"), labels = 1:3)]
## dtL.table2[, scenario := paste0(scenario, " (truth=",formatC(truth, format = "f", digits = 2),")")]
## dtL.table2[, scenario := factor(scenario, unique(scenario))]

## * prepare table

## ** process simulation results
dtW.table2 <- dcast(dtL.table2,
                    formula = scenario+truth+strategy~n,
                    value.var = c("estimate","std","coverage","widthCI","rejection.rate"))

dt.table2 <- cbind(dtW.table2[,.(scenario, truth = round(truth,3), variance = ifelse(scenario %in% c("A","B","C"),"modality", "region & modality"), strategy)],
                   Estimate = paste0(format.pval2(dtW.table2$estimate_24, digits = 3, eps = 1e-3), "/",
                                     format.pval2(dtW.table2$estimate_1000, digits = 3, eps = 1e-3)),
                   "Bias (%)" = paste0(format.pval2(100*(dtW.table2$estimate_24-dtW.table2$truth)/dtW.table2$truth, digits = 2, eps = 1e-3), "/",
                                       format.pval2(100*(dtW.table2$estimate_1000-dtW.table2$truth)/dtW.table2$truth, digits = 2, eps = 1e-3)),
                   Std = paste0(formatC(dtW.table2$std_24, digits = 3, format = "f"), "/",
                                formatC(dtW.table2$std_1000, digits = 3, format = "f")),
                   Coverage = paste0(formatC(dtW.table2$coverage_24, digits = 3, format = "f"), "/",
                                     formatC(dtW.table2$coverage_1000, digits = 3, format = "f")),
                   Width = paste0(formatC(dtW.table2$widthCI_24, digits = 3, format = "f"), "/",
                                  formatC(dtW.table2$widthCI_1000, digits = 3, format = "f")),
                   "Type 1 error" = paste0(formatC(dtW.table2$rejection.rate_24, digits = 3, format = "f"), "/",
                                           formatC(dtW.table2$rejection.rate_1000, digits = 3, format = "f")),
                   Power = paste0(formatC(dtW.table2$rejection.rate_24, digits = 3, format = "f"), "/",
                                  formatC(dtW.table2$rejection.rate_1000, digits = 3, format = "f"))
                   )
setkeyv(dt.table2, c("scenario","strategy"))
dt.table2[truth==0, "Bias (%)" := "."]
dt.table2[truth==0, "Power" := "."]
dt.table2[truth!=0, "Type 1 error" := "."]
## names(dt.table2)[3:7] <- paste(names(dt.table2)[3:7], "(n=24/n=1000)")
## options(width = 150)
dt.table2

##     scenario truth          variance strategy       Estimate        Bias (%)         Std    Coverage       Width Type 1 error       Power
##       <char> <num>            <char>   <fctr>         <char>          <char>      <char>      <char>      <char>       <char>      <char>
##  1:        A 0.500          modality        1    0.491/0.491     -1.74/-1.83 0.188/0.187 0.951/0.953 0.713/0.714            . 0.596/0.593
##  2:        A 0.250          modality        2    0.244/0.250     -2.27/-0.01 0.121/0.018 0.942/0.949 0.463/0.072            . 0.490/1.000
##  3:        A 0.250          modality        3    0.253/0.250       1.00/0.06 0.120/0.018 0.941/0.951 0.450/0.072            . 0.533/1.000
##  4:        B 0.000          modality        1    0.663/0.674               . 0.034/0.005 0.000/0.000 0.577/0.566  1.000/1.000           .
##  5:        B 0.000          modality        2  -0.001/<0.001               . 0.121/0.018 0.943/0.952 0.464/0.073  0.057/0.048           .
##  6:        B 0.000          modality        3  -0.001/<0.001               . 0.120/0.018 0.944/0.952 0.453/0.073  0.056/0.048           .
##  7:        C 0.500          modality        1    0.668/0.674     33.70/34.80 0.029/0.005 1.000/1.000 0.571/0.566            . 1.000/1.000
##  8:        C 0.250          modality        2    0.244/0.250     -2.27/-0.01 0.121/0.018 0.942/0.949 0.463/0.072            . 0.490/1.000
##  9:        C 0.250          modality        3    0.253/0.250       1.00/0.06 0.120/0.018 0.941/0.951 0.450/0.072            . 0.533/1.000
## 10:        D 0.000 region & modality        1   0.004/<0.001               . 0.243/0.245 0.948/0.947 0.889/0.888  0.051/0.052           .
## 11:        D 0.000 region & modality        2 <-0.001/<0.001               . 0.132/0.020 0.943/0.950 0.506/0.080  0.057/0.050           .
## 12:        D 0.000 region & modality        3 <-0.001/<0.001               . 0.120/0.019 0.948/0.953 0.461/0.074  0.052/0.048           .
## 13:        E 0.057 region & modality        1    0.663/0.674 1054.08/1072.97 0.038/0.006 0.000/0.000 0.577/0.566            . 1.000/1.000
## 14:        E 0.135 region & modality        2    0.132/0.135      -2.59/0.08 0.130/0.020 0.941/0.950 0.495/0.078            . 0.177/1.000
## 15:        E 0.135 region & modality        3    0.120/0.123    -11.40/-9.12 0.118/0.018 0.945/0.902 0.453/0.073            . 0.175/1.000

## NOTE 2: scenario D-E are extra simulations not described in the article with non-homogeneous variance across regions.
##         the LMM ignores this issue and is fitted as if all the regions have the same variance within modality.
##         Using partialCor(asl+pet~region, repetition = ~region|cimbi, data = df.joint, structure = "HCS", df = FALSE) does not make this restriction
##         but is not investigated in the simulation study (partly because it is much slower, partly because it is beyond the scope of the paper).

keep.col <- c("scenario","truth","strategy","Estimate","Bias (%)","Coverage","Type 1 error","Power")
table2 <- read_docx()
table2 <- body_add_flextable(table2,set_table_properties(autofit(flextable(as.data.frame(dt.table2[,.SD,.SDcols = keep.col]))), width = 0.9, layout = "autofit"))

## * export results 
print(table2, target = "tables/table2.docx")

##----------------------------------------------------------------------
### table2.R ends here
