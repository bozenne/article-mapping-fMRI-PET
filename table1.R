### table1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 18 2025 (10:54) 
## Version: 
## Last-Updated: mar 18 2025 (12:22) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(officer)
library(flextable)

## * load data
strategy2.cor <- readRDS("results/strategy2.rds")
strategy3.cor <- readRDS("results/strategy3.rds")
region.name <- readRDS("results/regionName.rds")

## * prepare table
n.digit <- 3
strategy23.cor <- rbind(strategy2.cor,strategy3.cor["marginal",c("estimate", "lower", "upper", "p.value")])

df.table1 <- cbind(" "=c("average"="Average","marginal"="Linear Mixed Model",region.name)[rownames(strategy23.cor)],
                   "Estimated correlation" = formatC(round(strategy23.cor$estimate,n.digit), format = "f", digits = n.digit),
                   "Confidence interval" = paste0("[",formatC(round(strategy23.cor$lower,n.digit), format = "f", digits = n.digit),
                                                  "; ",formatC(round(strategy23.cor$upper,n.digit), format = "f", digits = n.digit),"]"),
                   "P-value" = formatC(round(strategy23.cor$p.value,n.digit), format = "f", digits = n.digit))
## formatC keeps trailing 0

table1 <- read_docx()
table1 <- body_add_flextable(table1,set_table_properties(autofit(flextable(as.data.frame(df.table1))), width = 0.9, layout = "autofit"))

## * export
print(table1, target = "tables/table1.docx")

##----------------------------------------------------------------------
### table1.R ends here
