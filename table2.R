### table2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 18 2025 (11:03) 
## Version: 
## Last-Updated: mar 18 2025 (18:32) 
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

dtL.table2 <- dtS.sim[method %in% c("averageSignal","averageCor","lmm") & type %in% c("unclear","marginal")]
dtL.table2[, strategy := factor(method, levels = c("averageSignal","averageCor","lmm"), labels = 1:3)]
## dtL.table2[, scenario := paste0(scenario, " (truth=",formatC(truth, format = "f", digits = 2),")")]
## dtL.table2[, scenario := factor(scenario, unique(scenario))]

## * prepare table

## ** process simulation results
dtW.table2 <- dcast(dtL.table2,
                    formula = scenario+truth+strategy~n,
                    value.var = c("estimate","std","coverage","rejection.rate"))

dt.table2 <- cbind(dtW.table2[,.(scenario, truth = round(truth,3), variance = ifelse(scenario %in% c("A","B","C"),"modality", "region & modality"), strategy)],
                   Estimate = paste0(format.pval2(dtW.table2$estimate_26, digits = 1, eps = 1e-3), "/",format.pval2(dtW.table2$estimate_1000, digits = 3, eps = 1e-3)),
                   Std = paste0(formatC(dtW.table2$std_26, digits = 3, format = "f"), "/",formatC(dtW.table2$std_1000, digits = 3, format = "f")),
                   Coverage = paste0(formatC(dtW.table2$coverage_26, digits = 3, format = "f"), "/",formatC(dtW.table2$coverage_1000, digits = 3, format = "f")),
                   "Type 1 error" = paste0(formatC(dtW.table2$rejection.rate_26, digits = 3, format = "f"), "/",formatC(dtW.table2$rejection.rate_1000, digits = 3, format = "f")),
                   Power = paste0(formatC(dtW.table2$rejection.rate_26, digits = 3, format = "f"), "/",formatC(dtW.table2$rejection.rate_1000, digits = 3, format = "f"))
                   )
setkeyv(dt.table2, c("scenario","strategy"))
dt.table2[truth==0, "Power" := ""]
dt.table2[truth!=0, "Type 1 error" := ""]
## names(dt.table2)[3:7] <- paste(names(dt.table2)[3:7], "(n=26/n=1000)")
## options(width = 150)
dt.table2
## WARNING: no difference between 0<x<0.001 and -0.001<x<0 (both are coded <0.001)

##     scenario truth          variance strategy      Estimate         Std    Coverage Type 1 error       Power
##       <char> <num>            <char>   <fctr>        <char>      <char>      <char>       <char>      <char>
##  1:        A 0.500          modality        1   0.492/0.491 0.187/0.187 0.952/0.953              0.597/0.593
##  2:        A 0.250          modality        2   0.245/0.250 0.116/0.018 0.939/0.949              0.522/1.000
##  3:        A 0.250          modality        3   0.252/0.250 0.115/0.018 0.940/0.951              0.559/1.000
##  4:        B 0.000          modality        1   0.664/0.674 0.033/0.005 0.000/0.000  1.000/1.000            
##  5:        B 0.000          modality        2 -0.001/<0.001 0.116/0.018 0.943/0.952  0.057/0.048            
##  6:        B 0.000          modality        3 -0.001/<0.001 0.115/0.018 0.945/0.952  0.055/0.048            
##  7:        C 0.500          modality        1   0.669/0.674 0.028/0.005 1.000/1.000              1.000/1.000
##  8:        C 0.250          modality        2   0.245/0.250 0.116/0.018 0.939/0.949              0.522/1.000
##  9:        C 0.250          modality        3   0.252/0.250 0.115/0.018 0.940/0.951              0.559/1.000
## 10:        D 0.000 region & modality        1  0.005/<0.001 0.243/0.245 0.948/0.947  0.051/0.052            
## 11:        D 0.000 region & modality        2 -0.001/<0.001 0.127/0.020 0.942/0.950  0.058/0.050            
## 12:        D 0.000 region & modality        3 -0.001/<0.001 0.115/0.019 0.948/0.953  0.052/0.048            
## 13:        E 0.057 region & modality        1   0.664/0.674 0.037/0.006 0.000/0.000              1.000/1.000
## 14:        E 0.135 region & modality        2   0.132/0.135 0.124/0.020 0.943/0.950              0.190/1.000
## 15:        E 0.135 region & modality        3   0.120/0.123 0.113/0.018 0.945/0.902              0.183/1.000

## NOTE 2: scenario D-E are extra simulations not described in the article with non-homogeneous variance across regions.
##         the LMM ignores this issue and is fitted as if all the regions have the same variance within modality.
##         Using partialCor(asl+pet~region, repetition = ~region|cimbi, data = df.joint, structure = "HCS", df = FALSE) does not make this restriction
##         but is not investigated in the simulation study (partly because it is much slower, partly because it is beyond the scope of the paper).

keep.col <- c("scenario","truth","strategy","Estimate","Std","Coverage","Type 1 error","Power")
table2 <- read_docx()
table2 <- body_add_flextable(table2,set_table_properties(autofit(flextable(as.data.frame(dt.table2[scenario %in% c("A","B","C"),.SD,.SDcols = keep.col]))), width = 0.9, layout = "autofit"))

## ** format table for article

## * Result as text
rBias.table2 <- dtL.table2[scenario %in% c("A","C","E"),.("bias (%)" = 100*(estimate-truth)/truth), by = c("strategy","scenario","n")]
dcast(rBias.table2, formula = strategy+n~scenario)
##    strategy     n           A           C           E
##      <fctr> <num>       <num>       <num>       <num>
## 1:        1    26 -1.60518506 33.79706186 1055.770239
## 2:        1  1000 -1.83492753 34.79760080 1072.969076
## 3:        2    26 -2.09814713 -2.09814713   -2.556856
## 4:        2  1000 -0.01419166 -0.01419166    0.082628
## 5:        3    26  0.90051817  0.90051817  -11.450650
## 6:        3  1000  0.06139116  0.06139116   -9.124178


## * export results 
print(table2, target = "tables/table2.docx")

##----------------------------------------------------------------------
### table2.R ends here
