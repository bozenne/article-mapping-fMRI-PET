### cor.testIID.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  4 2025 (13:04) 
## Version: 
## Last-Updated: feb  7 2025 (10:57) 
##           By: Brice Ozenne
##     Update #: 35
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * cor.testIID (documentation)
##' @title Pearson's correlation test
##' @description Test for an association between two paried samples using Pearson's correlation.
##' Compared to \code{base::cor.test}, it also returns the standard error and influence information
##' @param x,y [numeric vector] paired measurments.
##' @param conf.level [numeric, 0-1] confidence level for the confidence interval.
##' @param se [character] type of standard error: model-based (\code{"model"}) or robust (\code{"robust"}).
##' @param transform [logical] should an atanh transformation be used for evaluating the confidence interval (followed by backtransformation) and the p-value.
##' @param fast [logical] approach used to evaluate the influence function.
##'
##' @details It seems that \code{base::cor.test} is computing the p-value and the confidence intervals based on different, possibly incompatible, uncertainty formulas: \itemize{
##' \item p-value: based student approximation (n-2 degrees of freedom) on the original scale.
##' \item confidence intervals: computed using a atanh transform, with normal approximation  and variance (1/(n-3)), then back-transformed.
##' }
##' 
##' 
##' @examples
##' library(lava)
##' set.seed(350)
##' df <- sim(lvm(Y~X), n = 10)
##'
##' #### default
##' cor.test(df$X,df$Y)
##' ## 0.6319313 ([0.003828097 0.902478809], p-value = 0.04998)
##'
##' #### equivalent to default (same CI, nearly identical p-value)
##' cor.testIID(df$X,df$Y)
##' ## 0.6319313 ([0.003828097 0.902478809], p-value = 0.04883)
##' 
##' #### robust alternative (same p-value, nearly identical CI)
##' cor.testIID(df$X,df$Y, transform = FALSE)
##' ## 0.681272 ([5.727936e-05 1.263805e+00] p-value = 0.04998)

## * cor.testIID (code)
cor.testIID <- function(x, y, conf.level = 0.95, se = "model", transform = TRUE, fast = TRUE){


    se <- match.arg(se, c("model","robust"))
    
    ## ** take care of missing data
    if(any(is.na(x)) || any(is.na(y))){
        out <- cor.testIID(x[!is.na(x) & !is.na(y)],y[!is.na(x) & !is.na(y)])
        iid.save <- out$iid
        out$iid <- matrix(0, nrow = length(x), ncol = 2, dimnames = list(NULL, c("original","atanh")))
        out$iid[!is.na(x) & !is.na(y),] <- iid.save
        return(out)
    }

    ## ** no missing data
    n.obs <- length(x)

    if(fast){
        ## *** standardize
        x.norm <- as.double(scale(x))
        y.norm <- as.double(scale(y))

        e.cor <- sum(x.norm*y.norm)/(n.obs-1)
        ## equation 3 in Robust Estimation and Outlier Detection with Correlation Coefficients (Devlin et al. 1975, https://doi.org/10.2307/2335508)
        IF.rho <- n.obs/(n.obs-1) * ((x.norm*y.norm) - e.cor*(x.norm^2+y.norm^2)/2)

    }else{
        ## *** original
        ## precompute
        mu.x <- mean(x)
        mu.x2 <- mean(x^2)
        mu.y <- mean(y)
        mu.y2 <- mean(y^2)
        mu.xy <- mean(x*y)
        x.center <- (x - mu.x)
        y.center <- (y - mu.y)

        ## point estimate
        e.num <- mu.xy-mu.x*mu.y
        e.denom1 <- (mu.x2-mu.x^2)
        e.denom2 <- (mu.y2-mu.y^2)
        e.cor <- e.num/sqrt(e.denom1*e.denom2)

        ## influence function
        ## https://cran.r-project.org/web/packages/lava/vignettes/influencefunction.html#example-pearson-correlation
        IF.num <- x.center*y.center - e.num
        IF.denom1 <- x.center^2 - e.denom1
        IF.denom2 <- y.center^2 - e.denom2
        IF.denom <- IF.denom1*e.denom2+IF.denom2*e.denom1
        IF.rho <- IF.num/sqrt(e.denom1*e.denom2) - 0.5 * e.num * IF.denom / (e.denom1*e.denom2)^(3/2)
    }
    
    ## ** export
    out <- cor.test(x, y, conf.level = conf.level)
    lower.level <- (1-conf.level)/2
    upper.level <- 1-lower.level

    out$estimate[] <- e.cor
    out$iid <- cbind(original = IF.rho, atanh = IF.rho/(1-e.cor^2))
    if(se == "model"){
        
        sigmaATANH.model <- (atanh(out$conf.int)[2] - atanh(out$estimate))/qnorm(upper.level)
        ## sigmaATANH.model - sqrt(1/(n.obs - 3))
        sigmaATANH.robust <- sqrt(sum(out$iid[,2]^2))/n.obs
        out$iid[,2] <- sigmaATANH.model/sigmaATANH.robust *out$iid[,2]
        sigmaTANH.model <- sigmaATANH.model*sqrt(1-out$estimate^2)*sqrt((n.obs-3)/(n.obs-2))
        ## correction by sqrt((n.obs-2)/(n.obs-3)) to retrieve the same p.value as cor.test
        ## sqrt( (1 - out$estimate^2) / (n.obs - 2) )
        sigmaTANH.robust <- sqrt(sum(out$iid[,1]^2))/n.obs
        out$iid[,1] <- sigmaTANH.model/sigmaTANH.robust *out$iid[,1]
    }
    out$se <- sqrt(colSums(out$iid^2))/n.obs

    if(transform){
        out$parameter[] <- Inf
        out$statistic[] <- atanh(out$estimate)/out$se[2]
        out$conf.int[] <- tanh(atanh(out$estimate) + qnorm(c(lower.level,upper.level)) * out$se[2])
    }else{
        out$parameter[] <- n.obs-2
        out$statistic[] <- out$estimate/out$se[1]
        out$conf.int[] <- out$estimate + qt(c(lower.level,upper.level), df = out$parameter[]) * out$se[1]
    }
    out$p.value[] <- 2*(1-pt(abs(out$statistic), df = out$parameter))
    return(out)
}

##----------------------------------------------------------------------
### cor.testIID.R ends here
