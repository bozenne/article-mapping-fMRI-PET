### FCT.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2025 (11:40) 
## Version: 
## Last-Updated: mar 21 2025 (10:11) 
##           By: Brice Ozenne
##     Update #: 126
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * simData (documentation)
##' @title Simulate regional fMRI/PET data
##' @description
##'
##' @param seed [interger, >0] Initial state of the random number generator. Used for reproducibility.
##' @param mu.PET [numeric vector] mean PET value in each region
##' @param mu.fMRI [numeric vector] mean fMRI value in each region. Should have the same length as \code{mu.PET}.
##' @param n.obs  [integer, >0] number of independent replicates.
##' @param rho.PET [numeric, between -1 and 1] correlation between any pair of PET regions. 
##' @param rho.fMRI [numeric, between -1 and 1] correlation between any pair of fMRI regions. 
##' @param rho.marginal [numeric, between -1 and 1] correlation between fMRI and PET in the same region belonging to two different individuals. 
##' @param rho.conditional [numeric, between -1 and 1] correlation between fMRI and PET in the same region belonging to the same individual. 
##' @param sigma.PET [numeric vector, >0] variability of the PET value in each region. Should have the same length as \code{mu.PET}.
##' @param sigma.fMRI [numeric vector, >0] variability of the fMRI value in each region. Should have the same length as \code{mu.PET}.
##'
##' @return A data.table in the long format.
##' 
##' @examples
##' df <- simData(seed = 1, n.obs = 50,
##'              mu.PET = 1:10,
##'              mu.fMRI = 1:10,
##'              rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5,
##'              sigma.PET = 3, sigma.fMRI = 5)
##' head(df)
##'
##' dfS <- aggregate(cbind(PET, fMRI) ~ region, data = df, FUN = "mean")
##' dfS

## * simData (code)
simData <- function(seed, mu.PET, mu.fMRI, n.obs, 
                    rho.PET, rho.fMRI, rho.marginal, rho.conditional,
                    sigma.PET = 1, sigma.fMRI = 5){

    require(data.table)
    require(mvtnorm)

    ## *** normalize user-input
    set.seed(seed)
    n.region <- length(mu.PET)
    if(length(mu.fMRI)!=n.region){
        stop("Mismatch between the length of argument \'mu.PET\' and \'mu.fMRI\'. \n")
    }
    if(length(sigma.fMRI)==1){sigma.fMRI <- rep(sigma.fMRI, n.region)}
    if(length(sigma.PET)==1){sigma.PET <- rep(sigma.PET, n.region)}
    
    ## *** prepare 
    mu <- c(mu.PET, mu.fMRI)
    sigma <- c(sigma.PET, sigma.fMRI)
    R.PET <- rho.PET + diag(1-rho.PET,n.region,n.region)
    R.fMRI <- rho.fMRI + diag(1-rho.fMRI,n.region,n.region)
    R.PETfMRI <- diag(rho.marginal,n.region,n.region)    
    R.PETfMRI[R.PETfMRI==0] <- rho.marginal - rho.conditional*sqrt((1-rho.PET)*(1-rho.fMRI))
    ## rho_conditional = (rho_marginal-rho)/sqrt(...)
    ## so rho = rho_marginal - rho_conditional * sqrt(...)

    R.all <- rbind(cbind(R.PET,R.PETfMRI),
                   cbind(R.PETfMRI,R.fMRI))
    Sigma <- tcrossprod(sigma)*R.all

    ## *** simulate
    dtW <- data.table(1:n.obs, rmvnorm(n.obs, mean = mu, sigma = Sigma))
    names(dtW) <- c("id",paste0("PET",1:n.region),paste0("fMRI",1:n.region))
    dtW$id <- as.factor(dtW$id)

    ## *** reshape
    out <- melt(dtW, id.vars = c("id"), measure.vars = patterns("PET","fMRI"),
                value.name = c("PET","fMRI"), variable.name = "region")

    ## *** export
    attr(out,"mu") <- mu
    attr(out,"Sigma") <- Sigma
    return(out)

}

## * runCor (documentation)
##' @title Analyse regional fMRI/PET data
##' @description Evaluate the correlation between the two modalities using different methods.
##'
##' @param formula [formula] 
##' @param data [data.frame] Dataset containing a column id, region, PET, and fMRI
##' @param df.lmm [logical] Should a Satterthwaite approximation be used to evaluate the degrees of freedom?
##' Can be quite conservative in small samples for marginal/latent correlation.
##' @param lmmH [logical] Should a linear mixed model with heteroscedastic residual variance be fitted?
##' Can be time consuming.
##'
##' @details Available methods:
##' \itemize{
##' \item{region}{Pearson's correlation between the average (over subjects) regional values.}
##' \item{id}{Average (over regions) of the region-specific Pearson's correlation between the two modalities.}
##' \item{idNorm}{Same as before but the signal in each modality and each region is centered and scaled.}
##' \item{lmm0}{Linear mixed model modeling modality dependent mean, variance, and correlation.}
##' \item{lmm}{Linear mixed model modeling a regional and modality dependent mean. The variance and correlation are modeled only dependent on the modality}
##' }
##'
##' @return A data.table with the method used to compute the correlation, type of correlation, and estimated correlation value.
##'
##' @examples
##' df <- simData(seed = 1, n.obs = 50,
##'              mu.PET = 1:10,
##'              mu.fMRI = 1:10,
##'              rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5,
##'              sigma.PET = 3, sigma.fMRI = 5)
##' runCor(PET + fMRI ~ region|id, data = df)

## * runCor (code)
runCor <- function(formula, data, df.lmm = FALSE, lmmH = FALSE){

    require(LMMstar)
    require(data.table)

    ## ** read formula
    if(!inherits(formula,"formula")){
        stop("The argument \'formula\' should contain a formula object. \n",
             "Something like PET + fMRI ~ region|id. \n")
    }
    ls.formula <- LMMstar:::formula2var(formula)
    if(!inherits(data,"data.frame")){
        stop("The argument \'data\' should contain a data.frame object. \n")
    }
    if(any(ls.formula$vars$all %in% names(data) == FALSE)){
        stop("Mismatch between argument \'formula\' and argument \'data\' \n",
             "Cannot find column \"",paste(setdiff(var.f,names(data)), collapse = "\", \""),"\" in argument \'data'. \n")
    }
    outcome <- ls.formula$vars$response
    if(length(outcome)!=2 || any(duplicated(outcome))){
        stop("There should be 2 distinct variables on the left hand side of argument \'formula\'. \n")
    }
    if(any(paste0(outcome,".norm") %in% names(data))){
        stop("The outcome variables with the extension .norm should not be in the data as these names will be used internally. \n")
    }
    region <- ls.formula$vars$time
    Uregion <- unique(data[[region]])
    if(length(region)!=1){
        stop("There should be exactly 1 regressor variables on the left hand side of argument \'formula\', the other should be a grouping variable. \n",
             "Something like PET + fMRI ~ region |id. \n")
    }
    cluster <- ls.formula$vars$cluster
    if(length(region)!=1){
        stop("There should be exactly 1 grouping variables on the left hand side of argument \'formula\'. \n",
             "Something like PET + fMRI ~ region |id. \n")
    }

    ## ** normalize data
    if(is.data.table(data)){
        data2 <- data.table::copy(data)
    }else{
        data2 <- data.table::as.data.table(data)
    }
    data2[, paste0(outcome[1],".norm") := scale(.SD[[outcome[1]]]), by = region]
    data2[, paste0(outcome[2],".norm") := scale(.SD[[outcome[2]]]), by = region]

    n.NNA <- data2[,.(NNA = sum(!is.na(.SD[[outcome[1]]]))+sum(!is.na(.SD[[outcome[2]]]))), by = cluster]
    if(any(n.NNA$NNA==0)){ ## remove clusters with only missing values
        data2 <- data2[data2[[cluster]] %in% n.NNA[NNA!=0,cluster]]
    }
    
    ## ** compute correlation
    ## correlation between the mean regional values
    cor.meanRegion <- cor.test(data2[,mean(.SD[[outcome[1]]], na.rm=TRUE), by = region][[2]],
                               data2[,mean(.SD[[outcome[2]]], na.rm=TRUE), by = region][[2]])

    ## average regional correlation
    ls.regionalCor <- by(data2,data2[[region]], function(iDF){
        cor.testIID(iDF[[outcome[1]]],iDF[[outcome[2]]])
    })
    mean.regionalCor <- mean(sapply(ls.regionalCor,"[[","estimate"))
    iidATANH.regionalCor <- rowMeans(do.call(cbind,lapply(ls.regionalCor,function(iLS){iLS$iid[,"atanh"]})))
    seATANH.regionalCor <- sqrt(sum(iidATANH.regionalCor^2))/length(iidATANH.regionalCor)
    pvalue.regionalCor <- 2*(1-pnorm(abs(atanh(mean.regionalCor)/seATANH.regionalCor)))
    ci.regionalCor <- tanh(atanh(mean.regionalCor) + qnorm(c(0.025,0.975)) * seATANH.regionalCor)
    
    ## average individual correlation
    cor.individual <- data2[,.(cor = cor(.SD[[outcome[1]]],.SD[[outcome[2]]], use = "pairwise.complete.obs")),by=cluster]$cor
    meanCor.individual <- mean(cor.individual)
    seCor.individual <- sd(cor.individual)/sqrt(length(cor.individual))
    seCorATANH.individual <- seCor.individual/(1-meanCor.individual^2)
    pvalue.individual <- 2*(1-pnorm(abs(atanh(meanCor.individual)/seCorATANH.individual)))
    ci.individual <- tanh(atanh(meanCor.individual) + qnorm(c(0.025,0.975)) * seCorATANH.individual)

    ## average individual correlation on normalized signal
    cor.individualNorm <- data2[,.(cor = cor(.SD[[paste0(outcome[1],".norm")]],.SD[[paste0(outcome[2],".norm")]], use = "pairwise.complete.obs")),by=cluster]$cor
    meanCor.individualNorm <- mean(cor.individualNorm)
    seCor.individualNorm <- sd(cor.individualNorm)/sqrt(length(cor.individualNorm))
    seCorATANH.individualNorm <- seCor.individualNorm/(1-meanCor.individualNorm^2)
    pvalue.individualNorm <- 2*(1-pnorm(abs(atanh(meanCor.individualNorm)/seCorATANH.individualNorm)))
    ci.individualNorm <- tanh(atanh(meanCor.individualNorm) + qnorm(c(0.025,0.975)) * seCorATANH.individualNorm)

    ## mixed model approach
    f.lmm0 <- reformulate(termlabels = "1", response = paste(outcome,collapse="+"))
    f.lmm <- reformulate(termlabels = region, response = paste(outcome,collapse="+"))
    f.rep <- as.formula(paste0("~",region,"|",cluster))
    
    e.lmm0 <- do.call(partialCor, args = list(f.lmm0, data = data2, repetition = f.rep, structure = "CS", df = df.lmm))
    e.lmm <- do.call(partialCor, args = list(f.lmm, data = data2, repetition = f.rep, structure = "CS", df = df.lmm))
    if(lmmH){
        e.lmmH <- do.call(partialCor, args = list(f.lmm, data = data2, repetition = f.rep, structure = "HCS", df = df.lmm))
    }
    ## eRho.lmm <- coef(attr(e.lmm,"lmm"), effect = "correlation")
    ## (eRho.lmm["rho(1,2,dt=0)"]-eRho.lmm["rho(1,2,dt=1)"])/sqrt(prod(1-eRho.lmm[c("rho1","rho2")]))

    ## ** export
    out <- rbind(data.table(method = "averageSignal", type = "unclear",
                            estimate = cor.meanRegion$estimate, lower = cor.meanRegion$conf.int[1], upper = cor.meanRegion$conf.int[2], p.value = cor.meanRegion$p.value),
                 data.table(method = "averageCor", type = "marginal",
                            estimate = mean.regionalCor, lower = ci.regionalCor[1], upper = ci.regionalCor[2], p.value = pvalue.regionalCor),
                 data.table(method = "averageId", type = "conditional",
                            estimate = meanCor.individual, lower = ci.individual[1], upper = ci.individual[2], p.value = pvalue.individual),
                 data.table(method = "averageIdNorm", type = "conditional",
                            estimate = meanCor.individualNorm, lower = ci.individualNorm[1], upper = ci.individualNorm[2], p.value = pvalue.individualNorm),
                 data.table(method = "lmm0", type = "marginal", e.lmm0["marginal",c("estimate","lower","upper","p.value")]),
                 data.table(method = "lmm0", type = "conditional", e.lmm0["conditional",c("estimate","lower","upper","p.value")]),
                 data.table(method = "lmm0", type = "latent", e.lmm0["latent",c("estimate","lower","upper","p.value")]),
                 data.table(method = "lmm", type = "marginal", e.lmm["marginal",c("estimate","lower","upper","p.value")]),
                 data.table(method = "lmm", type = "conditional", e.lmm["conditional",c("estimate","lower","upper","p.value")]),
                 data.table(method = "lmm", type = "latent", e.lmm["latent",c("estimate","lower","upper","p.value")])
                 )
    if(lmmH){
        out <- rbind(out,
                     data.table(method = "lmmH", type = "marginal", e.lmmH["marginal",c("estimate","lower","upper","p.value")]),
                     data.table(method = "lmmH", type = "conditional", e.lmmH["conditional",c("estimate","lower","upper","p.value")]),
                     data.table(method = "lmmH", type = "latent", e.lmmH["latent",c("estimate","lower","upper","p.value")]))
    }
    return(out)
}

## * plotCor
##' @title Graphical representation of the estimated correlations.
##' @description Graphical representation of the estimated correlations.
##'
##' @param data [data.frame] dataset containing a column id, region, PET, and fMRI. 
##' @param type.norm [1 or 2] if \code{2}, the mean and variance of the PET, as well as the variance of the fMRI normalized data is set to the one of the non-normalized data (across all regions).
##' @param widths [numeric vector of length 2,>0] width of the scatterplot and boxplot panel. Passed to \code{ggarrange}.
##' @param round.dotplot [integer,>0] rounding used to group points in the points in the dotplot on top of the boxplot. 
##' @param xspace.dotplot [numeric,>0] horizontal space used for the points in the dotplot on top of the boxplot. 
##' @param size.dotplot [numeric vector,>0] size of the points in the dotplot on top of the boxplot. 
##' 
##' @examples
##' df <- simData(seed = 1, n.obs = 50,
##'              mu.PET = 1:10,
##'              mu.fMRI = 1:10,
##'              rho.PET = 0.8, rho.fMRI = 0.6, rho.marginal = 0.25, rho.conditional = 0.5,
##'              sigma.PET = 3, sigma.fMRI = 5)
##'
##' plotCor(df, type.norm = 1)
##' plotCor(df, type.norm = 2)
##' 
plotCor <- function(data, type.norm, n.id = 4, widths = c(1.5,0.5), round.dotplot = 2, xspace.dotplot = 0.3, size.dotplot = c(6,4)){

    require(data.table)
    require(ggplot2)
    require(ggpubr)

    ## *** normalize data
    if(type.norm %in% 1:2 == FALSE){
        stop("Argument \'type.norm\' should be 1 or 2. \n")
    }
    
    if(is.data.table(data)){
        data2 <- copy(data)
        by <- list(region = c("region"),
                   id = c("id","normalized"),
                   group = c("group"),
                   facet = ~id)
    }else if(is.data.frame(data)){
        data2 <- as.data.table(data)
        by <- list(region = c("region"),
                   id = c("id","normalized"),
                   group = c("group"),
                   facet = ~id)
    }else if(is.list(data)){
        data2 <- as.data.table(do.call(rbind,lapply(1:length(data), function(iD){
            cbind(scenario = names(data)[iD], data[[iD]])
        })))
        by <- list(region = c("scenario","region"),
                   id = c("scenario","id","normalized"),
                   group = c("group","scenario"),
                   facet = scenario~id)
    }else{
        stop("Unknown format for argument \'data\'. \n")
    }

    data2[,c("id","normalized") := .(paste0("Subject ",.SD$id),FALSE)]
    data2.norm <- data.table::copy(data2)
    data2.norm[, c("PET","fMRI","normalized") := list(as.double(scale(.SD$PET)), as.double(scale(.SD$fMRI)), TRUE), by = eval(by$region)]
    if(type.norm == 2){
        if("scenario" %in% names(data2.norm)){
            for(iS in unique(data2.norm$scenario)){
                data2.norm$PET[data2.norm$scenario == iS] <- data2.norm$PET[data2.norm$scenario == iS]*sd(data2$PET[data2$scenario == iS]) + mean(data2$PET[data2$scenario == iS])
                data2.norm$fMRI[data2.norm$scenario == iS] <- data2.norm$fMRI[data2.norm$scenario == iS]*sd(data2$fMRI[data2$scenario == iS])
            }
        }else{
            data2.norm$PET <- data2.norm$PET*sd(data2$PET) + mean(data2$PET)
            data2.norm$fMRI <- data2.norm$fMRI*sd(data2$fMRI) 
        }
    }
    data2.average <- data2[,.(id = "Average", PET = mean(PET), fMRI = mean(fMRI), normalized = unique(normalized)), by = eval(by$region)]
    data2.normaverage <- data2.norm[,.(id = "Average", PET = mean(PET), fMRI = mean(fMRI), normalized = unique(normalized)), by = eval(by$region)]
    data2.all <- rbind(data2,data2.norm,data2.average[,.SD,.SDcols = names(data2)],data2.normaverage[,.SD,.SDcols = names(data2)])
    data2.all[, id := factor(id, levels = unique(id))]

    data2.cor <- data2.all[id != "Average", .(cor=cor(PET,fMRI)) , by = eval(by$id)]
    data2.cor$label <- gsub("Subject ","",data2.cor$id)
    data2.cor$id <- "Correlation"
    round.factor <- 1+9*(round.dotplot %% 1)/10
    round.dotplot <- floor(round.dotplot)
    data2.cor[, group := interaction(data2.cor$normalized,round(round.factor*data2.cor$cor,round.dotplot))]
    data2.cor[, myx := as.numeric(normalized)+1 + seq(0, by = xspace.dotplot, length.out = .N) - (.N-1)*xspace.dotplot/2, by = eval(by$group)]

    ## *** display
    keep.id <- levels(data2.all$id)[c(1:(n.id-1),nlevels(data2.all$id)-(0:1))]
    data2.subset <- data2.all[id %in% keep.id]
    data2.subset$id <- factor(data2.subset$id, levels = unique(data2.subset$id), labels = c(unique(as.character(data2.subset$id))[1:(n.id-2)], "...", unique(as.character(data2.subset$id))[n.id+0:1]))
    ggScat <- ggplot(data2.subset, aes(x = PET, y = fMRI, group = normalized, color = normalized))
    ggScat <- ggScat + facet_grid(by$facet) + geom_point(aes(shape=region)) + geom_smooth(method = "lm", se = FALSE)
    ggScat <- ggScat + scale_shape_manual(values= 1:nlevels(data2.all$region))
    ggScat <- ggScat + guides(shape = guide_legend(nrow = 1))

    ggBox <- ggplot(data2.cor, aes(y = cor, fill = normalized))
    ggBox <- ggBox + geom_boxplot(alpha = 0.5, aes(x = normalized))
    ggBox <- ggBox + facet_grid(by$facet) + geom_point(aes(x = myx), size = size.dotplot[1])
    ggBox <- ggBox + geom_text(aes(label = label, x = myx), color = "white", size = size.dotplot[2])
    ggBox <- ggBox + ylab("")
    ggBox <- ggBox + theme(axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank())

    ggScat <- ggScat + theme(text = element_text(size=12), 
                           axis.line = element_line(linewidth = 1.25),
                           axis.ticks = element_line(linewidth = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.key.size = unit(1,"line"))
    ggBox <- ggBox + theme(text = element_text(size=12), 
                           axis.line = element_line(linewidth = 1.25),
                           axis.ticks = element_line(linewidth = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.key.size = unit(1,"line"))
    
    ggAll <- ggarrange(ggScat, ggBox, nrow = 1, widths = widths, legend = "bottom", common.legend = TRUE)


    ## *** export
    return(ggAll)
}

## * Collect simulation results
loadRes <- function(path, tempo.file = FALSE, type = NULL,
                    export.attribute = NULL, trace = 2, space = "     "){
    all.files <- list.files(path)
    file.tempo <- grep("(tempo)",all.files,value = TRUE)
    file.final <- setdiff(all.files, file.tempo)

    if(tempo.file){
        file.read <- file.tempo
    }else{
        file.read <- file.final
    }
    if(!is.null(type)){
        file.read <- grep(pattern=type,x=file.read,value=TRUE)
    }

    n.file <- length(file.read)
    if(n.file==0){
        if(trace>1){
            cat(space,"No file found in ",path,". \n\n",sep="")
        }
        return(NULL)
    }
    myApply <- switch(as.character(as.logical(trace)),
                      "TRUE" = pbapply::pblapply,
                      "FALSE" = lapply)

    ls.out <- do.call(myApply, args = list(X = 1:n.file, FUN = function(iFile){
        if(grepl("\\.rds$",file.read[iFile])){
            iRead <- try(readRDS(file = file.path(path,file.read[iFile])))
        }else if(grepl("\\.rda$", file.read[iFile])){
            iRead <- try(load(file = file.path(path,file.read[iFile])))
            if(!inherits(iRead,"try-error")){
                iRead <- eval(parse(text = iRead))
            }
        }else{
            return(NULL)
        }
        if(inherits(iRead,"try-error")){
            return(NULL)
        }else{
            iOut <- cbind(data.table::as.data.table(iRead),
                          file = file.read[iFile])
            return(iOut)
        }
    }))

    out <- do.call(rbind, ls.out)

    Urow <- unique(sapply(ls.out, NROW))
    if(length(Urow)==1){
        cat(space,length(ls.out)," files with ",Urow," lines \n\n",sep="")
    }else{
        cat(space,length(ls.out)," files, with from ",min(Urow)," to ",max(Urow)," lines \n\n",sep="")
    }
    return(out)
}

##----------------------------------------------------------------------
### FCT.R ends here
