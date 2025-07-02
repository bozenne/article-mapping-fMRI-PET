### table2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 18 2025 (11:03) 
## Version: 
## Last-Updated: jul  2 2025 (15:44) 
##           By: Brice Ozenne
##     Update #: 53
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * function and packages
library(data.table)
library(officer)
library(flextable)
set_flextable_defaults(
    font.size = 5, theme_fun = theme_vanilla)

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

table2word <- function(data, n){

    if(length(n)==1){
        df.table <- data.frame(data$scenario,
                               character(factor(data$strategy, levels = c("1 (single cohort)","1 (two cohorts)","2","3"), labels = c("1.1","1.2","2","3"))),
                               ifelse(data$scenario=="A","No","Yes"),
                               round(data$truth,3),
                               format.pval2(data[[paste0("estimate_",n)]], digits = 3, eps = 1e-3),
                               format.pval2(100*(data[[paste0("estimate_",n)]]-data$truth)/data$truth, digits = 2, eps = 1e-3),
                               formatC(data[[paste0("coverage_",n)]], digits = 3, format = "f"), 
                               formatC(data[[paste0("rejection.rate_",n)]], digits = 3, format = "f"),
                               formatC(data[[paste0("rejection.rate_",n)]], digits = 3, format = "f"),
                               check.names = FALSE)
        names(df.table) <- c("Scenario","Strategy","Region confound","True correlation",
                             paste(c("Estimate correlation","Bias (%)","Coverage","Type-1 error","Statistical power"),n, sep="_N="),
                             )
        index. <- list(7,9,10)
    
    }else if(length(n)==2){

        df.table <- data.frame(data$scenario,
                               as.character(factor(data$strategy, levels = c("1 (single cohort)","1 (two cohorts)","2","3"), labels = c("1.1","1.2","2","3"))),
                               ifelse(data$scenario=="A","No","Yes"),
                               round(data$truth,3),
                               format.pval2(data[[paste0("estimate_",n[1])]], digits = 3, eps = 1e-3),
                               format.pval2(data[[paste0("estimate_",n[2])]], digits = 3, eps = 1e-3),
                               format.pval2(100*(data[[paste0("estimate_",n[1])]]-data$truth)/data$truth, digits = 2, eps = 1e-3),
                               format.pval2(100*(data[[paste0("estimate_",n[2])]]-data$truth)/data$truth, digits = 2, eps = 1e-3),
                               formatC(data[[paste0("coverage_",n[1])]], digits = 3, format = "f"), 
                               formatC(data[[paste0("coverage_",n[2])]], digits = 3, format = "f"), 
                               formatC(data[[paste0("rejection.rate_",n[1])]], digits = 3, format = "f"),
                               formatC(data[[paste0("rejection.rate_",n[2])]], digits = 3, format = "f"),
                               formatC(data[[paste0("rejection.rate_",n[1])]], digits = 3, format = "f"),
                               formatC(data[[paste0("rejection.rate_",n[2])]], digits = 3, format = "f"),
                               check.names = FALSE)
        names(df.table) <- c("Scenario","Strategy","Region confound","True correlation",
                             unlist(lapply(c("Estimate correlation","Bias (%)","Coverage","Type-1 error","Statistical power"), paste, n, sep="_N="))
                             )

        index. <- list(7:8,11:12,13:14)

    }else if(length(n)>2){
        stop("Argument \'n\' must have length at most 2. \n")
    }

    
    df.table <- df.table[order(df.table$Scenario,df.table$Strategy),]
    df.table$Scenario <- ifelse(duplicated(df.table$Scenario),"",df.table$Scenario)
    df.table[df.table[["True correlation"]]==0, index.[[1]]] <-  "--" ## Bias
    df.table[df.table[["True correlation"]]!=0, index.[[2]]] <- "--" ## Type 1
    df.table[df.table[["True correlation"]]==0, index.[[3]]] <- "--" ## Power

    flextable.table <- flextable(df.table)
    flextable.table <- theme_booktabs(flextable.table)
    
    flextable.table <- vline(flextable.table, j = 4, border = fp_border(color = "black", width = 1.5))
    flextable.table <- vline(flextable.table, j = c(6,8,10,12,14), border = fp_border(color = "black", width = 1))
    ## flextable.table <- hline(flextable.table, i = c(1:4)*4, border = brdr)
    flextable.table <- align(flextable.table, align = "center", part = "all")
    flextable.table <- set_header_df(flextable.table, mapping = data.frame(
                                                          col_keys = c("Region confound","True correlation","Estimated correlation"),
                                                          labels = c("Region\nconfound","True\nCorrelation","Estimated\ncorrelation"),
                                                          stringsAsFactors = FALSE
                                                      ))
    index.grey <- unlist(lapply(c(0,2,4), function(iIndex){iIndex*4+1:4}))
    index.white <- unlist(lapply(c(1,3), function(iIndex){iIndex*4+1:4}))
    flextable.table <- bg(flextable.table, i = index.grey, bg = gray.colors(1, start = 0.85, end = 0.85))
    flextable.table <- bg(flextable.table, i = index.white, bg = gray.colors(1, start = 0.99, end = 0.99))

    out <- body_add_flextable(read_docx(),
                              set_table_properties(autofit(separate_header(flextable.table, split = "_")),
                                                   width = 0.9, layout = "autofit"))
    return(out)
}

## * load data
dtS.sim <-  readRDS(file = "results/simulation-summary.rds")
## unique(dtS.sim$method)
## [1] "averageSignal"  "averageCor"     "averageId"      "averageIdNorm"  "lmm0"           "lmm"            "averageSignal2"
dtL.table2 <- dtS.sim[method %in% c("averageSignal","averageSignal2","averageCor","lmm") & type %in% c("unclear","marginal")]
dtL.table2[, strategy := factor(method, levels = c("averageSignal","averageSignal2","averageCor","lmm"), labels = c("1 (single cohort)","1 (two cohorts)","2","3"))]
## dtL.table2[, scenario := paste0(scenario, " (truth=",formatC(truth, format = "f", digits = 2),")")]
## dtL.table2[, scenario := factor(scenario, unique(scenario))]

## * prepare table

dtW.table2 <- dcast(dtL.table2,
                    formula = scenario+truth+strategy~n,
                    value.var = c("estimate","std","coverage","widthCI","rejection.rate"))

## ** simulation results (large-ish n)

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

##     scenario truth          variance          strategy       Estimate        Bias (%)         Std    Coverage       Width Type 1 error       Power
##       <char> <num>            <char>            <fctr>         <char>          <char>      <char>      <char>      <char>       <char>      <char>
##  1:        A 0.500          modality 1 (single cohort)    0.491/0.491     -1.74/-1.83 0.188/0.187 0.951/0.953 0.713/0.714            . 0.596/0.593
##  2:        A 0.500          modality   1 (two cohorts)   0.003/<0.001   -99.42/-99.95 0.242/0.244 0.432/0.433 0.889/0.888            . 0.047/0.051
##  3:        A 0.250          modality                 2    0.244/0.250     -2.27/-0.01 0.121/0.018 0.942/0.949 0.463/0.072            . 0.490/1.000
##  4:        A 0.250          modality                 3    0.253/0.250       1.00/0.06 0.120/0.018 0.941/0.951 0.450/0.072            . 0.533/1.000
##  5:        B 0.000          modality 1 (single cohort)    0.663/0.674               . 0.034/0.005 0.000/0.000 0.577/0.566  1.000/1.000           .
##  6:        B 0.000          modality   1 (two cohorts)    0.663/0.674               . 0.034/0.005 0.000/0.000 0.578/0.566  1.000/1.000           .
##  7:        B 0.000          modality                 2  -0.001/<0.001               . 0.121/0.018 0.943/0.952 0.464/0.073  0.057/0.048           .
##  8:        B 0.000          modality                 3  -0.001/<0.001               . 0.120/0.018 0.944/0.952 0.453/0.073  0.056/0.048           .
##  9:        C 0.500          modality 1 (single cohort)    0.668/0.674     33.70/34.80 0.029/0.005 1.000/1.000 0.571/0.566            . 1.000/1.000
## 10:        C 0.500          modality   1 (two cohorts)    0.663/0.674     32.55/34.76 0.034/0.005 1.000/1.000 0.578/0.566            . 1.000/1.000
## 11:        C 0.250          modality                 2    0.244/0.250     -2.27/-0.01 0.121/0.018 0.942/0.949 0.463/0.072            . 0.490/1.000
## 12:        C 0.250          modality                 3    0.253/0.250       1.00/0.06 0.120/0.018 0.941/0.951 0.450/0.072            . 0.533/1.000
## 13:        D 0.000 region & modality 1 (single cohort)   0.004/<0.001               . 0.243/0.245 0.948/0.947 0.889/0.888  0.051/0.052           .
## 14:        D 0.000 region & modality   1 (two cohorts) -0.002/<-0.001               . 0.245/0.244 0.951/0.947 0.888/0.889  0.048/0.052           .
## 15:        D 0.000 region & modality                 2 <-0.001/<0.001               . 0.132/0.020 0.943/0.950 0.506/0.080  0.057/0.050           .
## 16:        D 0.000 region & modality                 3 <-0.001/<0.001               . 0.120/0.019 0.948/0.953 0.461/0.074  0.052/0.048           .
## 17:        E 0.057 region & modality 1 (single cohort)    0.663/0.674 1054.08/1072.97 0.038/0.006 0.000/0.000 0.577/0.566            . 1.000/1.000
## 18:        E 0.057 region & modality   1 (two cohorts)    0.662/0.674 1053.00/1072.88 0.039/0.006 0.000/0.000 0.578/0.566            . 1.000/1.000
## 19:        E 0.135 region & modality                 2    0.132/0.135      -2.59/0.08 0.130/0.020 0.941/0.950 0.495/0.078            . 0.177/1.000
## 20:        E 0.135 region & modality                 3    0.120/0.123    -11.40/-9.12 0.118/0.018 0.945/0.902 0.453/0.073            . 0.175/1.000
##     scenario truth          variance          strategy       Estimate        Bias (%)         Std    Coverage       Width Type 1 error       Power

## NOTE 2: scenario D-E are extra simulations not described in the article with non-homogeneous variance across regions.
##         the LMM ignores this issue and is fitted as if all the regions have the same variance within modality.
##         Using partialCor(asl+pet~region, repetition = ~region|cimbi, data = df.joint, structure = "HCS", df = FALSE) does not make this restriction
##         but is not investigated in the simulation study (partly because it is much slower, partly because it is beyond the scope of the paper).

## ** process simulation results (small n)
dt.table2bis <- cbind(dtW.table2[,.(scenario, truth = round(truth,3), variance = ifelse(scenario %in% c("A","B","C"),"modality", "region & modality"), strategy)],
                      Estimate = paste0(format.pval2(dtW.table2$estimate_6, digits = 3, eps = 1e-3), "/",
                                        format.pval2(dtW.table2$estimate_12, digits = 3, eps = 1e-3)),
                      "Bias (%)" = paste0(format.pval2(100*(dtW.table2$estimate_6-dtW.table2$truth)/dtW.table2$truth, digits = 2, eps = 1e-3), "/",
                                          format.pval2(100*(dtW.table2$estimate_12-dtW.table2$truth)/dtW.table2$truth, digits = 2, eps = 1e-3)),
                      Std = paste0(formatC(dtW.table2$std_6, digits = 3, format = "f"), "/",
                                   formatC(dtW.table2$std_12, digits = 3, format = "f")),
                      Coverage = paste0(formatC(dtW.table2$coverage_6, digits = 3, format = "f"), "/",
                                        formatC(dtW.table2$coverage_12, digits = 3, format = "f")),
                      Width = paste0(formatC(dtW.table2$widthCI_6, digits = 3, format = "f"), "/",
                                     formatC(dtW.table2$widthCI_12, digits = 3, format = "f")),
                      "Type 1 error" = paste0(formatC(dtW.table2$rejection.rate_6, digits = 3, format = "f"), "/",
                                              formatC(dtW.table2$rejection.rate_12, digits = 3, format = "f")),
                      Power = paste0(formatC(dtW.table2$rejection.rate_6, digits = 3, format = "f"), "/",
                                     formatC(dtW.table2$rejection.rate_12, digits = 3, format = "f"))
                      )
setkeyv(dt.table2bis, c("scenario","strategy"))
dt.table2bis[, scenario := ifelse(duplicated(scenario),"",scenario)]
dt.table2bis[truth==0, "Bias (%)" := "."]
dt.table2bis[truth==0, "Power" := "."]
dt.table2bis[truth!=0, "Type 1 error" := "."]
## names(dt.table2)[3:7] <- paste(names(dt.table2)[3:7], "(n=24/n=1000)")
## options(width = 150)
dt.table2bis

##     scenario truth          variance          strategy       Estimate        Bias (%)         Std    Coverage       Width Type 1 error       Power
##       <char> <num>            <char>            <fctr>         <char>          <char>      <char>      <char>      <char>       <char>      <char>
##  1:        A 0.500          modality 1 (single cohort)    0.486/0.489     -2.74/-2.11 0.189/0.188 0.950/0.951 0.717/0.715            . 0.584/0.589
##  2:        A 0.500          modality   1 (two cohorts)  <0.001/<0.001   -99.85/-99.95 0.244/0.241 0.433/0.431 0.889/0.890            . 0.051/0.049
##  3:        A 0.250          modality                 2    0.227/0.240     -9.24/-3.92 0.250/0.173 0.931/0.930 0.943/0.651            . 0.190/0.299
##  4:        A 0.250          modality                 3    0.264/0.257       5.45/2.90 0.242/0.170 0.900/0.927 0.801/0.612            . 0.274/0.354
##  5:        B 0.000          modality 1 (single cohort)    0.631/0.652               . 0.066/0.047 0.009/0.000 0.610/0.589  0.990/1.000           .
##  6:        B 0.000          modality   1 (two cohorts)    0.630/0.652               . 0.066/0.047 0.011/0.000 0.611/0.589  0.988/1.000           .
##  7:        B 0.000          modality                 2 <-0.001/<0.001               . 0.248/0.173 0.946/0.936 0.944/0.651  0.054/0.064           .
##  8:        B 0.000          modality                 3 <-0.001/<0.001               . 0.243/0.171 0.924/0.936 0.805/0.615  0.076/0.064           .
##  9:        C 0.500          modality 1 (single cohort)    0.653/0.663     30.52/32.63 0.057/0.041 0.994/1.000 0.587/0.576            . 0.997/1.000
## 10:        C 0.500          modality   1 (two cohorts)    0.630/0.652     25.98/30.33 0.066/0.047 0.994/0.998 0.611/0.589            . 0.989/1.000
## 11:        C 0.250          modality                 2    0.227/0.240     -9.24/-3.92 0.250/0.173 0.931/0.930 0.943/0.651            . 0.190/0.299
## 12:        C 0.250          modality                 3    0.264/0.257       5.45/2.90 0.242/0.170 0.900/0.927 0.801/0.612            . 0.274/0.354
## 13:        D 0.000 region & modality 1 (single cohort)  -0.003/<0.001               . 0.243/0.246 0.950/0.948 0.889/0.888  0.049/0.051           .
## 14:        D 0.000 region & modality   1 (two cohorts) <-0.001/-0.001               . 0.247/0.245 0.943/0.943 0.887/0.888  0.056/0.056           .
## 15:        D 0.000 region & modality                 2  <0.001/<0.001               . 0.272/0.190 0.944/0.934 1.013/0.707  0.056/0.066           .
## 16:        D 0.000 region & modality                 3  <0.001/<0.001               . 0.245/0.172 0.927/0.936 0.823/0.627  0.073/0.064           .
## 17:        E 0.057 region & modality 1 (single cohort)    0.632/0.652 1001.04/1036.24 0.076/0.053 0.064/0.009 0.607/0.588            . 0.976/0.998
## 18:        E 0.057 region & modality   1 (two cohorts)    0.628/0.651  993.53/1032.94 0.076/0.054 0.073/0.009 0.612/0.590            . 0.971/0.999
## 19:        E 0.135 region & modality                 2    0.124/0.130     -8.72/-3.91 0.268/0.186 0.941/0.934 0.998/0.693            . 0.077/0.125
## 20:        E 0.135 region & modality                 3    0.115/0.119   -15.23/-12.27 0.240/0.169 0.925/0.938 0.813/0.617            . 0.102/0.122
##     scenario truth          variance          strategy       Estimate        Bias (%)         Std    Coverage       Width Type 1 error       Power

## NOTE 2: scenario D-E are extra simulations not described in the article with non-homogeneous variance across regions.
##         the LMM ignores this issue and is fitted as if all the regions have the same variance within modality.
##         Using partialCor(asl+pet~region, repetition = ~region|cimbi, data = df.joint, structure = "HCS", df = FALSE) does not make this restriction
##         but is not investigated in the simulation study (partly because it is much slower, partly because it is beyond the scope of the paper).


## * export results 
flextable2 <- table2word(dtW.table2, n = c(24,1000))
print(flextable2, target = "tables/table2.docx")
flextable2SM <- table2word(dtW.table2, n = c(6,12))
print(flextable2SM, target = "tables/table2SM.docx")

##----------------------------------------------------------------------
### table2.R ends here
