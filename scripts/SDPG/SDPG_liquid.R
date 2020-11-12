##### SDPG LIQUID GROWTH CURVES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/05/2020

##### INITIALIZATION
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(locfit)
library(growthrates)
library(stringr)
library(zoo)

path.out <- 'figs/SDPG/'
path.out.gc <- 'figs/SDPG/GC/'
expt_name <- 'SDPG' 

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### LOAD DATA AND CLEAN
df <- as.data.frame(readxl::read_excel("rawdata/SDPG/SDPG_LiquidGrowth.xlsx", sheet = 2))
df <- df[!(df$Arm == 'SDA' & df$Replicate %in% c('R1','R2')),]

map <- as.data.frame(readxl::read_excel("rawdata/SDPG/SDPG_96well.xlsx"))
map$orf_name[map$orf_name == 'YGR296W'] <- "YGR269W"
map$orf_name[map$orf_name == 'REF'] <- 'BF_control'

##### ANALYZE USING GROWTHRATES
growth_dat <- NULL
growth_res <- NULL
growth_exp <- NULL
growth_fit <- NULL

min.od <- 0.001
min.t <- 0
for (a in unique(df$Arm)) {
  for (r in unique(df$Replicate[df$Arm == a])) {
    temp <- df[df$Arm == a & df$Replicate == r,]
    temp$Time <- seq(15,15*dim(temp)[1],15)
    temp <- temp[,colSums(temp > 0.001) >= dim(temp)[1]*0.9]
    temp[temp <= 0] <- 0.001
    
    for (i in 4:dim(temp)[2]) {
      # lo <- loess.smooth(temp$Time, log(temp[,i]), span = 0.15, evaluation = length(temp[,i]))
      # fit0 <- fit_easylinear(lo$x, exp(lo$y), h=10, quota = 1);
      fit0 <- fit_easylinear(temp$Time, temp[,i], h = 20, quota = 1);
      AUC <- sum(diff(temp$Time)*rollmean(temp[,i],2))

      temp_res <- data.frame(arm = a, replicate = r, maxgr = coef(fit0)[[3]],
                             dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]],
                             auc = AUC, ood = temp[2,i], sod = temp[length(temp[,i]), i],
                             orf_name = as.character(map$orf_name[map$well == colnames(temp[i])]),
                             bio_rep = map$replicate[map$well == colnames(temp[i])])
      growth_res <- rbind(growth_res, temp_res)

      temp_dat <- slot(fit0,'obs')
      temp_dat$od <- temp[,i]
      temp_dat$arm <- a
      temp_dat$replicate <- r
      temp_dat$maxgr <- coef(fit0)[[3]]
      temp_dat$dtime <- log(2)/coef(fit0)[[3]]
      temp_dat$ltime <- coef(fit0)[[4]]
      temp_dat$auc <- AUC
      temp_dat$ood <- temp[2,i]
      temp_dat$sod <- temp[length(temp[,i]), i]
      temp_dat$orf_name <- as.character(map$orf_name[map$well == colnames(temp[i])])
      temp_dat$bio_rep <- map$replicate[map$well == colnames(temp[i])]
      growth_dat <- rbind(growth_dat, temp_dat)

      temp_exp <- slot(fit0,'obs')[slot(fit0,'ndx'),]
      temp_exp$arm <- a
      temp_exp$replicate <- r
      temp_exp$maxgr <- coef(fit0)[[3]]
      temp_exp$dtime <- log(2)/coef(fit0)[[3]]
      temp_exp$ltime <- coef(fit0)[[4]]
      temp_exp$auc <- AUC
      temp_exp$ood <- temp[2,i]
      temp_exp$sod <- temp[length(temp[,i]), i]
      temp_exp$orf_name <- as.character(map$orf_name[map$well == colnames(temp[i])])
      temp_exp$bio_rep <- map$replicate[map$well == colnames(temp[i])]
      growth_exp <- rbind(growth_exp, temp_exp)

      temp_fit <- data.frame(time = seq(0,4000,100),
                             y = predict(slot(fit0, 'fit'),
                                         newdata = data.frame(x = seq(0,4000,100))))
      temp_fit$arm <- a
      temp_fit$replicate <- r
      temp_fit$maxgr <- coef(fit0)[[3]]
      temp_fit$dtime <- log(2)/coef(fit0)[[3]]
      temp_fit$ltime <- coef(fit0)[[4]]
      temp_fit$auc <- AUC
      temp_fit$ood <- temp[2,i]
      temp_fit$sod <- temp[length(temp[,i]), i]
      temp_fit$orf_name <- as.character(map$orf_name[map$well == colnames(temp[i])])
      temp_fit$bio_rep <- map$replicate[map$well == colnames(temp[i])]
      growth_fit <- rbind(growth_fit, temp_fit)
    }
  }
}
growth_res$arm <- factor(growth_res$arm, levels = c('GLU','CAS','SDA'))
growth_dat$arm <- factor(growth_dat$arm, levels = c('GLU','CAS','SDA'))
growth_exp$arm <- factor(growth_exp$arm, levels = c('GLU','CAS','SDA'))
growth_fit$arm <- factor(growth_fit$arm, levels = c('GLU','CAS','SDA'))

growth_res <- growth_res[!(growth_res$orf_name %in% c('YHR021W-A','BLANK','BORDER')),]
growth_dat <- growth_dat[!(growth_dat$orf_name %in% c('YHR021W-A','BLANK','BORDER')),]
growth_exp <- growth_exp[!(growth_exp$orf_name %in% c('YHR021W-A','BLANK','BORDER')),]
growth_fit <- growth_fit[!(growth_fit$orf_name %in% c('YHR021W-A','BLANK','BORDER')),]

# head(growth_exp)
# head(growth_dat)
# head(growth_exp)
# head(growth_fit)

save(growth_res, growth_dat, growth_exp, growth_fit,
     file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/growthrates_results.RData')


# #### GROWTH CRUVES
# gc_all <- ggplot() +
#   geom_point(data = growth_dat, aes(x = time, y = od)) +
#   # geom_line(data = growth_dat, aes(x = time, y = y)) +
#   geom_point(data = growth_exp, aes(x = time, y = y), col = 'red') +
#   geom_line(data = growth_fit, aes(x = time, y = exp(y))) +
#   geom_text(data = growth_res, aes(x = 2000, y = 0.025, label = round(dtime,2)), size = 3) +
#   geom_text(data = growth_res, aes(x = 2000, y = 0.05, label = round(auc,2)), size = 3) +
#   scale_y_log10() +
#   coord_cartesian(ylim = c(0.01,6)) +
#   facet_wrap(.~arm*orf_name*replicate*bio_rep) +
#   theme_linedraw() +
#   theme(axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.key.size = unit(3, "mm"),
#         legend.position = "bottom",
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf('%sGROWTHCURVES_ALL.png',path.out), gc_all,
#        height = 1000, width = 1000, units = 'mm',
#        dpi = 300, limitsize = F)
# 
# 
# ##### AUC RESULTS
# ggplot(growth_res,
#        aes(x = orf_name, y = auc)) +
#   geom_boxplot() +
#   geom_point(aes(col = replicate)) +
#   coord_flip() +
#   labs(y = 'AUC') +
#   facet_wrap(~arm*replicate)
# 
# ##### DOUBLING TIME RESULTS
# ggplot(growth_res,
#        aes(x = orf_name, y = dtime)) +
#   geom_boxplot() +
#   geom_point(aes(col = replicate)) +
#   coord_flip() +
#   labs(y = 'Doubling Time (mins)') +
#   facet_wrap(~arm*replicate)
# 
