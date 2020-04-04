

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


###### SPECIFICITY AND SENSITIVITY
spe.data <- NULL
sen.data <- NULL

expt_name <- '4C4_FS_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SPECIFICITY.RData", expt_name))
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SENSITIVITY.RData", expt_name))
spe.data <- rbind(spe.data, spe.all)
sen.data <- rbind(sen.data, sen.all)

expt_name <- '4C4_TR_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SPECIFICITY.RData", expt_name))
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SENSITIVITY.RData", expt_name))
spe.data <- rbind(spe.data, spe.all)
sen.data <- rbind(sen.data, sen.all)

expt_name <- '4C4_TRBL_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SPECIFICITY.RData", expt_name))
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SENSITIVITY.RData", expt_name))
spe.data <- rbind(spe.data, spe.all)
sen.data <- rbind(sen.data, sen.all)

expt_name <- '4C4_TRBLBR_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SPECIFICITY.RData", expt_name))
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_SENSITIVITY.RData", expt_name))
spe.data <- rbind(spe.data, spe.all)
sen.data <- rbind(sen.data, sen.all)


save(sen.data,
     file = sprintf('%sSENSITIVITY.RData',out_path))
save(spe.data,
     file = sprintf('%sSPECIFICITY.RData',out_path))


##### FITNESS RESULTS
rmse.rr <- NULL

rmse <- NULL
expt_name <- '4C4_FS_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))
temp <- fit.all[fit.all$cont_hrs == fit.all$hours,]
i <- 1
for (r in unique(temp$rep)) {
  if (r %in% seq(2,16,2)) {
    dat.bg <- temp[temp$rep == r,]
    for (hr in unique(dat.bg$hours)) {
      # for (a in seq(1,dim(dat.bg[dat.bg$hours == hr,])[[1]]/(6144*length(unique(dat.bg$plate))))) {
        # rmse$atmpt[i] <- a
        rmse$rep[i] <- r
        rmse$hours[i] <- hr
        rmse$pix_avg[i] <- mean(dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0], na.rm = T)
        rmse$rmse[i] <- sqrt(mean((dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0] -
                                     dat.bg$bg[dat.bg$hours == hr & dat.bg$average > 0])^2,na.rm = T))
        
        # rmse$pix_avg[i] <- mean(dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)], na.rm = T)
        # rmse$rmse[i] <- sqrt(mean((dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)] -
        #                              dat.bg$bg[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)])^2,na.rm = T))
        i <- i + 1
      # }
    }
  }
}
rmse$ref <- 1/4
rmse <- data.frame(rmse)
rmse.rr <- rbind(rmse.rr, rmse)

rmse <- NULL
expt_name <- '4C4_TR_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))
temp <- fit.all[fit.all$cont_hrs == fit.all$hours,]
i <- 1
for (r in unique(temp$rep)) {
  if (r %in% seq(2,16,2)) {
    dat.bg <- temp[temp$rep == r,]
    for (hr in unique(dat.bg$hours)) {
      # for (a in seq(1,dim(dat.bg[dat.bg$hours == hr,])[[1]]/(6144*length(unique(dat.bg$plate))))) {
      # rmse$atmpt[i] <- a
      rmse$rep[i] <- r
      rmse$hours[i] <- hr
      rmse$pix_avg[i] <- mean(dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0], na.rm = T)
      rmse$rmse[i] <- sqrt(mean((dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0] -
                                   dat.bg$bg[dat.bg$hours == hr & dat.bg$average > 0])^2,na.rm = T))
      
      # rmse$pix_avg[i] <- mean(dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)], na.rm = T)
      # rmse$rmse[i] <- sqrt(mean((dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)] -
      #                              dat.bg$bg[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)])^2,na.rm = T))
      i <- i + 1
      # }
    }
  }
}
rmse$ref <- 1/4*3/4
rmse <- data.frame(rmse)
rmse.rr <- rbind(rmse.rr, rmse)

rmse <- NULL
expt_name <- '4C4_TRBL_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))
temp <- fit.all[fit.all$cont_hrs == fit.all$hours,]
i <- 1
for (r in unique(temp$rep)) {
  if (r %in% seq(2,16,2)) {
    dat.bg <- temp[temp$rep == r,]
    for (hr in unique(dat.bg$hours)) {
      # for (a in seq(1,dim(dat.bg[dat.bg$hours == hr,])[[1]]/(6144*length(unique(dat.bg$plate))))) {
      # rmse$atmpt[i] <- a
      rmse$rep[i] <- r
      rmse$hours[i] <- hr
      rmse$pix_avg[i] <- mean(dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0], na.rm = T)
      rmse$rmse[i] <- sqrt(mean((dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0] -
                                   dat.bg$bg[dat.bg$hours == hr & dat.bg$average > 0])^2,na.rm = T))
      
      # rmse$pix_avg[i] <- mean(dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)], na.rm = T)
      # rmse$rmse[i] <- sqrt(mean((dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)] -
      #                              dat.bg$bg[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)])^2,na.rm = T))
      i <- i + 1
      # }
    }
  }
}
rmse$ref <- 1/4*2/4
rmse <- data.frame(rmse)
rmse.rr <- rbind(rmse.rr, rmse)

rmse <- NULL
expt_name <- '4C4_TRBLBR_CC'
load(sprintf("/home/sbp29/R/Projects/adaptivefitness/output/%s_FITNESS.RData", expt_name))
temp <- fit.all[fit.all$cont_hrs == fit.all$hours,]
i <- 1
for (r in unique(temp$rep)) {
  if (r %in% seq(2,16,2)) {
    dat.bg <- temp[temp$rep == r,]
    for (hr in unique(dat.bg$hours)) {
      # for (a in seq(1,dim(dat.bg[dat.bg$hours == hr,])[[1]]/(6144*length(unique(dat.bg$plate))))) {
      # rmse$atmpt[i] <- a
      rmse$rep[i] <- r
      rmse$hours[i] <- hr
      rmse$pix_avg[i] <- mean(dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0], na.rm = T)
      rmse$rmse[i] <- sqrt(mean((dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0] -
                                   dat.bg$bg[dat.bg$hours == hr & dat.bg$average > 0])^2,na.rm = T))
      
      # rmse$pix_avg[i] <- mean(dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)], na.rm = T)
      # rmse$rmse[i] <- sqrt(mean((dat.bg$average[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)] -
      #                              dat.bg$bg[dat.bg$hours == hr & dat.bg$average > 0][seq(6144*length(unique(dat.bg$plate))*(a-1),6144*length(unique(dat.bg$plate))*a,1)])^2,na.rm = T))
      i <- i + 1
      # }
    }
  }
}
rmse$ref <- 1/4*1/4
rmse <- data.frame(rmse)
rmse.rr <- rbind(rmse.rr, rmse)
rmse.rr$per <- rmse.rr$rmse/rmse.rr$pix_avg * 100

save(rmse.rr,
     file = sprintf('%sRMSE_REFREP.RData',out_path))
