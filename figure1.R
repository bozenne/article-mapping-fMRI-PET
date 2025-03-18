### figure1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 18 2025 (10:16) 
## Version: 
## Last-Updated: mar 18 2025 (10:52) 
##           By: Brice Ozenne
##     Update #: 4
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(ggplot2)
library(ggrepel)

## * load data
strategy1.cor <- readRDS("results/strategy1.rds")
df.2cohorts <- readRDS("results/df2cohorts.rds")
region.name <- readRDS("results/regionName.rds")

## * prepare figure
nsig <- 2
sub.text.corr <- paste0("Pearson's rho: ", signif(strategy1.cor$estimate,nsig),
                        ', 95% CI: ', paste(signif(strategy1.cor$conf.int,nsig), collapse = '; '),
                        ', p = ', signif(strategy1.cor$p.value,nsig))

figure1 <- ggplot(cbind(df.2cohorts, region.full = region.name[df.2cohorts$region]), aes(x = pet, y = asl)) 
figure1 <- figure1 + geom_smooth(method = 'lm', alpha = .8, color = 'black') 
figure1 <- figure1 + geom_point(aes(color = region)) 
figure1 <- figure1 + theme_classic() 
figure1 <- figure1 + theme(axis.text = element_text(size=20),
                           axis.title = element_text(size=24,face="bold"), 
                           legend.position = 'none',
                           plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
                           plot.subtitle = element_text(hjust = 0.5, size = 20))
figure1 <- figure1 + labs(x = '[11C]Cimbi-36 BPND', y = 'Baseline CBF (ml/100g/min)', title = '5-HT2AR vs. CBF', subtitle = sub.text.corr)
figure1 <- figure1 + scale_x_continuous(breaks = seq(0.2,1.6,0.2))
figure1 <- figure1 + scale_y_continuous(limits = c(35,70), breaks = seq(35,70,5))
figure1 <- figure1 + geom_label_repel(aes(label = region.full, 
                                          fill = region), size = 5)
figure1

## * export
ggsave(figure1, filename = file.path("figures","figure1.pdf"), units = "in", width = 12.5, height = 7.5)
ggsave(figure1, filename = file.path("figures","figure1.png"), units = "in", width = 12.5, height = 7.5, dpi=300)



##----------------------------------------------------------------------
### figure1.R ends here
