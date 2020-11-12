##### SDPG LIQUID GROWTH CURVES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/05/2020

##### INITIALIZATION
library(tidyverse)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(stringr)
library(zoo)

path.out <- 'figs/SDPG/'
path.out.gc <- 'figs/SDPG/GC/'
expt_name <- 'SDPG' 

##### FUNCTIONS FROM ALEEZA GERSTEIN (Aleeza.Gerstein@umanitoba.ca)
# generic function to calculate growth rate 
nderiv <- function(fit, x, eps=1e-5) {
  (predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)
}

spline.slope <- function(x, y, n=300, eps=1e-5) {
  max(nderiv(loess(log(y) ~ x, degree=1, span=0.1),
             seq(min(x), max(x), length=n)), na.rm=TRUE)
}
  
spline.slope.dat <- function(d, ...) {
  sapply(d[-1], spline.slope, x=d$t, ...)
}
  
##### MY FUNCTION
reset_time <- function(data) {
  for (a in unique(data$Arm)) {
    for (r in unique(data$Replicate[data$Arm == a])) {
      data$Time <- as.character(data$Time)
      data$Time[data$Arm == a & data$Replicate == r] <-
        seq(15,15*dim(data[data$Arm == a & data$Replicate == r,])[1],15)
      
    }
  }
  data$Time <- as.numeric(data$Time)
  return(data)
}

growth_results <- function(data,map) {
  res <- NULL
  for (a in unique(data$Arm)) {
    for (r in unique(data$Replicate[data$Arm == a])) {
      temp <- data[data$Arm == a & data$Replicate == r,]
      temp <- temp[,c(3:dim(temp)[2])]
      # temp$Time <- seq(15,15*dim(temp)[1],15)
      # data$Time[data$Arm == a & data$Replicate == r] <- seq(15,15*dim(temp)[1],15)
      colnames(temp)[1] <- 't'
      
      temp <- temp[,colSums(temp > 0.001) >= dim(temp)[1]*0.9]
      temp[temp <= 0] <- 0.001
      lgr <- spline.slope.dat(temp)
      lgr <- data.frame(lgr)
      lgr$well <- rownames(lgr)  
      #add growth rate on to well information - use the name of your wells document here
      temp_res <- merge(map, lgr, by = 'well', all = T)
      temp_res$dtime <- log(2)/temp_res$lgr
      temp_res$Arm <- a
      temp_res$exp_rep <- r
      res <- rbind(res, temp_res)
    }
  }
  return(res)
}


##### LOADING OUR DATA
data_blk <- as.data.frame(readxl::read_excel("rawdata/SDPG/SDPG_LiquidGrowth.xlsx", sheet = 1))
data_plc <- as.data.frame(readxl::read_excel("rawdata/SDPG/SDPG_LiquidGrowth.xlsx", sheet = 2))
map <- as.data.frame(readxl::read_excel("rawdata/SDPG/SDPG_96well.xlsx"))
map$orf_name[map$orf_name == 'YGR296W'] <- "YGR269W"
map$orf_name[map$orf_name == 'REF'] <- 'BF_control'

data_blk$Arm <- factor(data_blk$Arm, levels = c('GLU','CAS','SDA'))
data_blk$Replicate <- factor(data_blk$Replicate, levels = c('R1','R2','R3'))
data_blk <- reset_time(data_blk)
data_plc$Arm <- factor(data_plc$Arm, levels = c('GLU','CAS','SDA'))
data_plc$Replicate <- factor(data_plc$Replicate, levels = c('R1','R2','R3'))
data_plc <- reset_time(data_plc)
# calculate the growth rate
head(data_plc)

res_plc <- growth_results(data_plc, map)
res_plc <- res_plc[res_plc$orf_name != 'BORDER' &  res_plc$orf_name != 'BLANK' & res_plc$orf_name != 'YHR021W-A',]
res_plc <- res_plc[!(res_plc$Arm == 'SDA' & res_plc$exp_rep %in% c('R1','R2')),]
res_plc$Arm <- factor(res_plc$Arm, levels = c('GLU','CAS','SDA'))

res_blk <- growth_results(data_blk, map)
res_blk <- res_blk[res_blk$orf_name != 'BORDER' &  res_blk$orf_name != 'BLANK' & res_blk$orf_name != 'YHR021W-A',]
res_blk <- res_blk[!(res_blk$Arm == 'SDA' & res_blk$exp_rep %in% c('R1','R2')),]
res_blk$Arm <- factor(res_blk$Arm, levels = c('GLU','CAS','SDA'))

##### SAVING THE RESULTS
save(res_plc, res_blk,
     file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/aleeza_results.RData')

##### STATS
res <- res_plc

# liq_stats <- data.frame(compare_means(dtime ~ orf_name, res,
#                          method = "t.test",
#                          p.adjust.method = 'bonferroni',
#                          group.by = c("Arm", "exp_rep"), ref.group = "REF"))
liq_stats <- NULL
i <- 1
for (a in unique(res$Arm)) {
  for (r in unique(res$exp_rep[res$Arm == a])) {
    ref_dtime <- res$dtime[res$orf_name == 'REF' & res$Arm == a & res$exp_rep == r]
    for (o in unique(res$orf_name[res$Arm == a & res$exp_rep == r & res$orf_name != 'REF'])) {
      orf_dtime <- res$dtime[res$orf_name == o & res$Arm == a & res$exp_rep == r]
      stat_res <- t.test(orf_dtime, ref_dtime)
      liq_stats$Arm[i] <- a
      liq_stats$Replicate[i] <- r
      liq_stats$orf_name[i] <- o
      liq_stats$p[i] <- stat_res$p.value
      liq_stats$tstat[i] <- stat_res$statistic[[1]]
      i <- i + 1
      # temp_cd <- cliff.delta(orf_dtime, ref_dtime, conf.level=.95,
      #                        use.unbiased=TRUE, use.normal=FALSE, return.dm=FALSE)
      # liq_stats$cliff.delta[liq_stats$group2 == o & liq_stats$Arm == a] <- temp_cd$estimate
      # liq_stats$magnitude[liq_stats$group2 == o & liq_stats$Arm == a] <- as.character(temp_cd$magnitude)
    }
  }
}
liq_stats <- data.frame(liq_stats)
head(liq_stats)

liq_stats$phenotype <- NULL
liq_stats$phenotype[liq_stats$p <= 0.05 & liq_stats$tstat < 0] <- 'Beneficial'
liq_stats$phenotype[liq_stats$p <= 0.05 & liq_stats$tstat > 0] <- 'Deleterious'
liq_stats$phenotype[is.na(liq_stats$phenotype)] <- 'Neutral'

liq_stats$Arm <- factor(liq_stats$Arm, levels = c('GLU','CAS','SDA'))

##### PLOT
liq_pheno_mat <- ggplot(res,
       aes(x = orf_name, y = dtime)) +
  # geom_point() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1) +
  labs(y = 'Doubling Time (mins)',
       x = 'ORFs') +
  # stat_compare_means(method = 't.test',
  #                    ref.group = 'REF',
  #                    label = "p.signif",
  #                    vjust = 0.5) +
  coord_flip() +
  facet_wrap(.~Arm*exp_rep, nrow = 7) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

liq_pheno_res <- ggplot(liq_stats) +
  geom_text(aes(y = orf_name, x = 1, label = phenotype, col = phenotype), size = 3) +
  facet_wrap(.~Arm*Replicate, nrow = 7) +
  scale_y_discrete(limits = sort(unique(res$orf_name))) +
  scale_color_manual(name = '',
                     breaks = c('Deleterious','Neutral','Beneficial'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107'),
                     guide = F) +
  labs(x = 'Liquid Phenotype',
       y = "ORFs") +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        panel.background = element_rect(fill = 'grey60'),
        axis.title = element_text(size = titles),
        axis.title.y = element_blank(),
        # axis.text = element_text(size = txt),
        axis.text = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

liq_pheno <- ggpubr::ggarrange(liq_pheno_mat, liq_pheno_res,
                                nrow = 1,
                                align = 'h',
                                widths = c(4,1),
                                common.legend = T,
                                legend = 'bottom')
ggsave(sprintf("%sLIQUID_PHENO.jpg",out_path),liq_pheno,
       height = 400, width = 100, units = 'mm',
       dpi = 300)

##### GROWTH CURVES
data_plc2 <- melt(data_plc, id.vars = c('Arm','Replicate','Time'))
data_plc2 <- merge(data_plc2, map, by.x = 'variable', by.y = 'well', all = T)

liq_gc <- ggplot(data_plc2[!(data_plc2$orf_name %in% c('BLANK','YHR021W-A','BORDER')) &
                   !(data_plc2$Arm == 'SDA' & data_plc2$Replicate == 'R2'),],
       aes(x = Time, y = value)) +
  geom_line(aes(col = Replicate, linetype = as.character(replicate))) +
  scale_color_discrete(name = 'Expt. Rep.') +
  scale_linetype_discrete(name = 'Tech. Rep.') +
  # scale_y_log10() +
  # geom_line(data = data_plc2[data_plc2$orf_name == 'REF' & data_plc2$Arm == 'GLU',],
  #           aes(x = Time, y = value), col = 'black') +
  facet_wrap(.~Arm*orf_name,ncol = 7) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        # axis.title.y = element_blank(),
        # axis.text = element_text(size = txt),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sLIQUID_GC.jpg",out_path),liq_gc,
       height = 300, width = 300, units = 'mm',
       dpi = 500)
