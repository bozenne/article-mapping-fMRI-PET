### figure2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 18 2025 (10:26) 
## Version: 
## Last-Updated: mar 23 2025 (11:45) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(ggplot2)

## * load data
strategy2.cor <- readRDS("results/strategy2.rds")
strategy3.cor <- readRDS("results/strategy3H.rds")
region.name <- readRDS("results/regionName.rds")


## * prepare figure
df.figure2 <- cbind(region = c(rownames(strategy2.cor),"lmm"),
                    region.full = c(c("average" = "Average", region.name)[rownames(strategy2.cor)],"Linear mixed model"),
                    region.lmm = c(c("average" = "Average", region.name)[rownames(strategy2.cor)],"Linear mixed model"),
                    rbind(strategy2.cor,strategy3.cor["marginal",c("estimate", "lower", "upper", "p.value")])
                    )
df.figure2$region.full <- factor(df.figure2$region.full, levels = unique(df.figure2$region.full))
df.figure2$style <- df.figure2$region %in% c("average","lmm")


figure2 <- ggplot(df.figure2, aes(x = region.full, y =  estimate, ymin = lower, ymax = upper) )
figure2 <- figure2 + geom_point(aes(shape = style, size = style)) + geom_linerange(aes(linewidth = style))
figure2 <- figure2 + geom_hline(yintercept = 0, linetype = "dashed")
figure2 <- figure2 + coord_cartesian(ylim = c(-1,1))
figure2 <- figure2 + labs(x="",y="Pearson's correlation coefficient") + guides(linewidth = "none", shape = "none", size = "none")
figure2 <- figure2 + scale_size_manual(breaks = c(FALSE,TRUE), values = c(2,3))
figure2 <- figure2 + scale_shape_manual(breaks = c(FALSE,TRUE), values = c(19,15))
figure2 <- figure2 + scale_linewidth_manual(breaks = c(FALSE,TRUE), values = c(1,1.5))
figure2 <- figure2 + theme(text = element_text(size=15), ,
                           axis.text.x = element_text(angle = 35, size = 8, vjust = 1, hjust = 1),                                     
                           axis.line = element_line(linewidth = 1.25),
                           axis.ticks = element_line(linewidth = 2),
                           axis.ticks.length=unit(.25, "cm"))
figure2 <- figure2 + scale_x_discrete(breaks = levels(df.figure2$region.full),
                                      labels = c(setdiff(levels(df.figure2$region.full),c("Average","Linear mixed model")),
                                                 expression(bold("Average")),expression(bold("Linear Mixed Model"))))
figure2

## * export
ggsave(figure2, filename = file.path("figures","figure2.pdf"), width = 8, height = 5)
ggsave(figure2, filename = file.path("figures","figure2.png"), width = 8, height = 5)

##----------------------------------------------------------------------
### figure2.R ends here
