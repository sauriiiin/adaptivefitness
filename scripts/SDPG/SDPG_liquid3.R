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
library(growthcurver)
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
df$Time <- as.numeric(df$Time)

map <- as.data.frame(readxl::read_excel("rawdata/SDPG/SDPG_96well.xlsx"))
map$orf_name[map$orf_name == 'YGR296W'] <- "YGR269W"
map$orf_name[map$orf_name == 'REF'] <- 'BF_control'

##### ANALYZE USING GROWTHCURVER
growth_crv <- NULL
models_all <- NULL

min.od <- 0.001
min.t <- 0
for (a in unique(df$Arm)) {
  for (r in unique(df$Replicate[df$Arm == a])) {
    temp <- df[df$Arm == a & df$Replicate == r,]
    temp$Time <- seq(15,15*dim(temp)[1],15)
    df$Time[df$Arm == a & df$Replicate == r] <- temp$Time
    temp <- temp[,colSums(temp > 0.001) >= dim(temp)[1]*0.9]
    temp[temp <= 0] <- 0.001
    
    temp_time <- temp$Time
    temp_ods <- temp[4:ncol(temp)]
    colnames(temp_ods) <- c(paste(temp$Arm[1], temp$Replicate[1], colnames(temp_ods), sep = '_'))
    temp_mdl <- lapply(temp_ods, function(x) SummarizeGrowth(temp_time, x))
    models_all <- append(models_all, temp_mdl)
    
    temp_crv <- SummarizeGrowthByPlate(temp[,3:dim(temp)[2]])
    temp_crv$arm <- a
    temp_crv$replicate <- r
    growth_crv <- rbind(growth_crv, temp_crv)
  }
}
growth_crv$arm <- factor(growth_crv$arm, levels = c('GLU','CAS','SDA'))
growth_crv <- merge(growth_crv, map, by.x = 'sample', by.y = 'well', all.x = T)
growth_crv <- growth_crv[!(growth_crv$orf_name %in% c('YHR021W-A','BLANK','BORDER')),]
colnames(growth_crv) <- c("well", "k", "n0", "r", "t_mid", "t_gen", "auc_l", "auc_e",
                          "sigma", "note", "arm", "replicate", "orf_name", "bio_rep")

df_pred <- df
df_pred[,4:ncol(df_pred)] <- 0
for (a in unique(df_pred$Arm)) {
  for (r in unique(df_pred$Replicate[df_pred$Arm == a])) {
    for (w in names(df_pred[,4:ncol(df_pred)])) {
      if (!is.null(models_all[[c(paste(a, r, w, sep = '_'))]])) {
        df_pred[df_pred$Arm == a & df_pred$Replicate == r, w] <-
          predict(models_all[[c(paste(a, r, w, sep = '_'))]]$model)
        
      }
    }
  }
}
head(df_pred)

melt1 <- melt(df, id.vars = c("Arm","Replicate","Time"), variable.name = "well", value.name = "od")
melt2 <- melt(df_pred, id.vars = c("Arm","Replicate","Time"), variable.name = "well", value.name = "pred.od")
df_final <- cbind(melt1, pred.od=melt2[,5])
df_final <- merge(df_final, map, by = "well", all = T)

df_final <- df_final[!(df_final$orf_name %in% c('YHR021W-A','BLANK','BORDER')),]
df_final$Arm <- factor(df_final$Arm, levels = c('GLU','CAS','SDA'))
colnames(df_final) <- c("well", "arm", "replicate", "time", "od", "pred.od", "orf_name", "bio_rep")

head(df_final)
head(growth_crv)

##### SAVING ALL DATA
save(growth_crv, df_final, models_all,
     file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/growthcurver_results.RData')

##### GROWTH CURVES
gc_all <- ggplot() +
  geom_point(data = df_final, aes(x = time, y = od), col = '#303F9F') +
  geom_line(data = df_final, aes(x = time, y = pred.od), col = '#E91E63') +
  geom_text(data = growth_crv, aes(x = 2000, y = 0.25, label = round(t_gen,2)), size = 3) +
  geom_text(data = growth_crv, aes(x = 2000, y = 1, label = round(auc_l,2)), size = 3) +
  labs(x = 'Time',
       y = 'OD') +
  coord_cartesian(ylim = c(0.01,6)) +
  facet_wrap(.~arm*orf_name*replicate*bio_rep, ncol = 18) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sGROWTHCURVES_ALL2.png',path.out), gc_all,
       height = 1000, width = 1000, units = 'mm',
       dpi = 300, limitsize = F)

gc_all2 <- ggplot() +
  geom_line(data = df_final, aes(x = time, y = pred.od, col = as.factor(bio_rep))) +
  geom_text(data = growth_crv, aes(x = 2000, y = 0.25,
                                   label = round(t_gen,2),
                                   col = as.factor(bio_rep)),
            position = position_dodge(width = 1500), size = 3) +
  geom_text(data = growth_crv, aes(x = 2000, y = 1,
                                   label = round(auc_l,2),
                                   col = as.factor(bio_rep)),
            position = position_dodge(width = 1500), size = 3) +
  # scale_y_continuous(trans = 'pseudo_log') +
  labs(x = 'Time',
       y = 'OD') +
  coord_cartesian(ylim = c(0.01,6)) +
  facet_wrap(.~arm*orf_name*replicate, ncol = 9) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sGROWTHCURVES_ALL3.png',path.out), gc_all2,
       height = 1000, width = 1000, units = 'mm',
       dpi = 300, limitsize = F)


##### DOUBLING TIME RESULTS
liq_dtime <- ggplot(growth_crv,
                    aes(x = orf_name, y = t_gen)) +
  geom_boxplot() +
  geom_point(aes(col = replicate.x)) +
  coord_flip(ylim = c(50,350)) +
  scale_color_discrete(name = 'Replicate') +
  labs(y = 'Doubling Time (mins)',
       x = 'ORFs') +
  stat_compare_means(method = 't.test', ref.group = "BF_control",
                     label = "p.signif", label.y = 350,
                     size = 1.5, angle = 0, vjust = 0.5) +
  # stat_compare_means(aes(group = orf_name),
  #                    label = "p.signif", label.y = 320,
  #                    size = 1.5, hide.ns = T) +
  facet_wrap(~arm) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sLIQUID_DTIME.png',path.out), liq_dtime,
       height = 150, width = 300, units = 'mm',
       dpi = 300, limitsize = F)


