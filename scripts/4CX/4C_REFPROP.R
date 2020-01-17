out_path = 'figs/lid_paper/4C4/';

data <- NULL

expt_name <- '4C4_FS_UP1'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))
dat.cnt.all$ref_prop <- 1/4
data <- rbind(data, dat.cnt.all)

expt_name <- '4C4_TR_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))
dat.cnt.all$ref_prop <- 1/4 * 3/4
data <- rbind(data, dat.cnt.all)

expt_name <- '4C4_TRBL_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))
dat.cnt.all$ref_prop <- 1/4 * 1/2
data <- rbind(data, dat.cnt.all)

expt_name <- '4C4_TRBLBR_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))
dat.cnt.all$ref_prop <- 1/4 * 1/4
data <- rbind(data, dat.cnt.all)

data$power[data$hours < data$cont_hrs] <- data$Deleterious_p[data$hours < data$cont_hrs]/910 * 100
data$power[data$hours > data$cont_hrs] <- data$Beneficial_p[data$hours > data$cont_hrs]/910 * 100

data$abs_cen <- abs(1-data$cen) * 100

ggplot(data[round(data$abs_cen) == 5,]) +
  # geom_boxplot(aes(x = rep, y = power, group = rep), fill = 'grey90') +
  geom_smooth(aes(x = rep, y = power, col = as.character(ref_prop*100)), method = 'loess', se = T, lwd = 1.2) +
  labs(title = 'Sensitivity in detecting 5% fitness effects',
       subtitle = 'with change in number of replicates and references') +
  scale_x_continuous(name = 'No. of Replicates',
                     breaks = seq(0,16,2),
                     minor_breaks = seq(0,16,1)) +
  scale_y_continuous(name = 'Sensitivity',
                     breaks = seq(0,100,10),
                     minor_breaks = seq(0,105,5)) +
  scale_color_discrete(name = 'Reference\nProportion') +
  coord_cartesian(ylim = c(0,100)) +
  # facet_wrap(.~ref_prop) +
  theme_linedraw()
ggsave(sprintf('%s4C4_REFREP.png',out_path),
       height = 7, width = 7,
       dpi = 300)
