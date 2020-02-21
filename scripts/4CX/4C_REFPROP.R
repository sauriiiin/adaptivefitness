

out_path = 'figs/lid_paper/4C4/';

refrep <- NULL

expt_name <- '4C4_FS_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))
dat.cnt.all$ref_prop <- 1/4
refrep <- rbind(refrep, dat.cnt.all)

expt_name <- '4C4_TR_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))
dat.cnt.all$ref_prop <- 1/4 * 3/4
refrep <- rbind(refrep, dat.cnt.all)

expt_name <- '4C4_TRBL_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))
dat.cnt.all$ref_prop <- 1/4 * 1/2
refrep <- rbind(refrep, dat.cnt.all)

expt_name <- '4C4_TRBLBR_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_DAT_CNT.RData", expt_name))
dat.cnt.all$ref_prop <- 1/4 * 1/4
refrep <- rbind(refrep, dat.cnt.all)

save(refrep,
     file = sprintf('%sREFREPDATA.RData',out_path))

refrep$power[refrep$hours < refrep$cont_hrs] <- refrep$Deleterious_p[refrep$hours < refrep$cont_hrs]/910 * 100
refrep$power[refrep$hours > refrep$cont_hrs] <- refrep$Beneficial_p[refrep$hours > refrep$cont_hrs]/910 * 100

refrep$abs_cen <- abs(1-refrep$cen) * 100

ggplot(refrep[round(refrep$abs_cen) <= 5,]) +
  # geom_boxplot(aes(x = rep, y = power, group = rep), fill = 'grey90') +
  geom_smooth(aes(x = rep, y = power, col = as.character(ref_prop*100)), method = 'loess', se = T, lwd = 1.2) +
  # labs(title = 'Sensitivity in detecting 5% fitness effects',
  #      subtitle = 'with change in number of replicates and references') +
  scale_x_continuous(name = 'No. of Replicates',
                     breaks = seq(0,16,2),
                     minor_breaks = seq(0,16,1)) +
  scale_y_continuous(name = 'Sensitivity',
                     breaks = seq(0,100,10),
                     minor_breaks = seq(0,105,5)) +
  scale_color_discrete(name = 'Reference\nProportion',
                       breaks = c('6.25','12.5','18.75','25'),
                       labels = c('06.25%','12.50%','18.75%','25.00%')) +
  coord_cartesian(ylim = c(0,100)) +
  # facet_wrap(.~ref_prop) +
  theme_linedraw() +
  # theme(legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 6)))
ggsave(sprintf('%s4C4_REFREP.png',out_path),
       height = 5, width = 5.5,
       dpi = 300)
