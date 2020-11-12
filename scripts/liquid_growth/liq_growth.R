##### LIQUID GROWTH CURVES CONTROL
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 11/11/2020

##### INITIALIZATION
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(growthcurver)
library(stringr)
library(zoo)

path.out <- 'figs/liq_growth/'

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### LOAD DATA AND CLEAN
df <- as.data.frame(readxl::read_excel("rawdata/liq_growth/liquidgrowth_control.xlsx", sheet = 1,
                                       col_types = c("text", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                     "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                     "numeric", "numeric", "numeric", "numeric", "numeric")))
df$Arm <- factor(df$Arm, levels = c('GLU','GAL'))
samples <- as.data.frame(readxl::read_excel("rawdata/liq_growth/liquidgrowth_control.xlsx", sheet = 2))
head(df)

##### PLOT RAW DATA
df2 <- melt(df, id.vars = c('Arm','Replicate','Time'), variable.name = 'sample', value.name = 'OD')
df2 <- merge(samples, df2, by = 'sample')

p <- ggplot(df2,
       aes(x = Time, y = OD, col = sample)) +
  geom_line() +
  facet_wrap(.~Arm) +
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
p
##### ANALYZE USING GROWTHCURVER
growth_crv <- NULL
# models_all <- NULL

min.od <- 0.001
min.t <- 0
for (a in unique(df$Arm)) {
  for (r in unique(df$Replicate[df$Arm == a])) {
    temp <- df[df$Arm == a & df$Replicate == r,]
    # temp$Time <- seq(15,15*dim(temp)[1],15)
    # df$Time[df$Arm == a & df$Replicate == r] <- temp$Time

    temp_time <- temp$Time
    temp_ods <- temp[4:ncol(temp)]
    colnames(temp_ods) <- c(paste(temp$Arm[1], temp$Replicate[1], colnames(temp_ods), sep = '_'))
    # temp_mdl <- lapply(temp_ods, function(x) SummarizeGrowth(temp_time, x))
    # models_all <- append(models_all, temp_mdl)

    temp_crv <- SummarizeGrowthByPlate(temp[,3:dim(temp)[2]])
    temp_crv$Arm <- a
    temp_crv$Replicate <- r
    growth_crv <- rbind(growth_crv, temp_crv)
  }
}
growth_crv$Arm <- factor(growth_crv$Arm, levels = c('GLU','GAL'))
growth_crv <- merge(samples, growth_crv, by = c('sample'))
head(growth_crv)

##### PLOT RESULTS
res_dtime <- ggplot(growth_crv,
       aes(x = strain, y = t_gen)) +
  geom_point(aes(col = as.factor(bio_rep), shape = tech_rep)) +
  scale_color_discrete(name = 'bioRep') +
  scale_shape_discrete(name = 'techRep') +
  labs(x = 'Strains',
       y = 'Doubling Time (mins)') +
  facet_wrap(.~Arm, ncol = 1) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = c(0.075,0.72),
        legend.background = element_rect(size=0.2, linetype="solid", 
                                         colour ="black"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

res_auc <- ggplot(growth_crv,
                    aes(x = strain, y = auc_l)) +
  geom_point(aes(col = as.factor(bio_rep), shape = tech_rep)) +
  scale_color_discrete(name = 'bioRep') +
  scale_shape_discrete(name = 'techRep') +
  labs(x = 'Strains',
       y = 'Area Under the Curve (AUC)') +
  facet_wrap(.~Arm, ncol = 1) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = c(0.075,0.72),
        legend.background = element_rect(size=0.2, linetype="solid", 
                                         colour ="black"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

crvs <- ggplot(df2,
       aes(x = Time, y = OD)) +
  geom_line(aes(col = as.factor(bio_rep), linetype = tech_rep)) +
  facet_wrap(.~Arm*strain, ncol = 4) +
  scale_color_discrete(name = 'bioRep') +
  scale_linetype_discrete(name = 'techRep') +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = c(0.075,0.72),
        legend.background = element_rect(size=0.2, linetype="solid", 
                                         colour ="black"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

res_all <- ggpubr::ggarrange(crvs, res_dtime, res_auc,
                  nrow = 1,
                  widths = c(3,1,1))
ggsave(sprintf('%sRESULTS.png',path.out), res_all,
       height = 200, width = 550, units = 'mm',
       dpi = 300, limitsize = F)

