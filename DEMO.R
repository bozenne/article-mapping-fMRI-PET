## DEMO.R --- 
##----------------------------------------------------------------------
## Author: Patrick Fisher
## Created: March 24 2025
## Version: 
## Last-Updated: 
##           By: 
##     Update #: 
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## *load data
df_CBF <- read.csv2('source/df_CBF_long_pruned__pvelab.csv')
df_demo <- read.csv2('source/MR45_MR001_ASL.csv')

## * unique ids
cimbi <- unique(df_CBF$CIMBI)

## * population age mean, sd, and range
age_mean <- mean(sapply(cimbi, function(i){return(unique(df_demo[df_demo$CIMBI==i, 'age']))}))
age_sd <- sd(sapply(cimbi, function(i){return(unique(df_demo[df_demo$CIMBI==i, 'age']))}))
age_range <- range(sapply(cimbi, function(i){return(unique(df_demo[df_demo$CIMBI==i, 'age']))}))
writeLines(paste0('age: ', signif(age_mean,3), ' +- ', signif(age_sd, 3), ' [', paste0(age_range, collapse = ';'), ']'))
## age: 32.5 +- 7.94 [24;58]


## * Sex table
sex_table <- table(sapply(cimbi, function(i){return(unique(df_demo[df_demo$CIMBI==i, 'sex']))}))
sex_table
## female   male 
## 10     14 


## * Scanner table
scanner_table <- table(sapply(cimbi, function(i){return(unique(tolower(df_demo[df_demo$CIMBI==i, 'scanner'])))}))
scanner_table
## mr001  mr45 
## 9    15 
