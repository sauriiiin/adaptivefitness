##### POWER DYNAMICS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/23/2019

##### INITIALIZE
library(ggplot2)
out_path = 'figs/lid_paper/';
dat.dir <- "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA1_LID/"

##### GET DATA
stats.files <- list.files(path = dat.dir,
                         pattern = "P.csv", recursive = TRUE)
fit.files <- list.files(path = dat.dir,
                        pattern = "S.csv", recursive = TRUE)
hours <- as.numeric(substr(stats.files,9,10))
dat.all <- NULL

##### PUTTING IT TOGETHER
for (i in 1:length(hours)) {
  hr <- hours[i]
  dat.stats <- read.csv(paste0(dat.dir,stats.files[i]),na.strings = "NaN")
  dat.stats$cont_hrs <- hours[i]
  dat.fit <- read.csv(paste0(dat.dir,fit.files[i]),na.strings = "NaN")
  cont.mean <- mean(dat.fit$fitness[dat.fit$hours == hr & dat.fit$orf_name == 'BF_control' & !is.na(dat.fit$fitness)])
  dat.stats$es <- dat.stats$cs_mean/cont.mean
  dat.stats$pthresh <- quantile(sort(dat.stats$p[dat.stats$hours == hr]),.05)
  dat.all <- rbind(dat.all,dat.stats)
}

# save(dat.all, file = "/home/sbp29/R/Projects/proto_plots/rawdata/4C3_GA1_CS/4C3_GA1_CS_DAT_ALL.RData")
effect_size <- sort(unique(round(dat.all$es,2)))
dat.pow <- NULL

##### POWER CALCULATIONS
for (es in effect_size) {
  N <- sum(dat.all$es > es-5e-03 & dat.all$es < es+5e-03)
  if (N > 0) {
    TP <- sum(dat.all$es > es-5e-03 & dat.all$es < es+5e-03 & dat.all$p <= dat.all$pthresh & dat.all$hours != dat.all$cont_hrs)
    FP <- sum(dat.all$es > es-5e-03 & dat.all$es < es+5e-03 & dat.all$p <= dat.all$pthresh & dat.all$hours == dat.all$cont_hrs)
    FPR <- FP/sum(dat.all$es > es-5e-03 & dat.all$es < es+5e-03 & dat.all$hours == dat.all$cont_hrs) * 100
    TN <- sum(dat.all$es > es-5e-03 & dat.all$es < es+5e-03 & dat.all$p > dat.all$pthresh & dat.all$hours == dat.all$cont_hrs)
    FN <- sum(dat.all$es > es-5e-03 & dat.all$es < es+5e-03 & dat.all$p > dat.all$pthresh & dat.all$hours != dat.all$cont_hrs)
    POW <- TP/N * 100
    SEN <- TP/sum(dat.all$es > es-5e-03 & dat.all$es < es+5e-03 & dat.all$hours != dat.all$cont_hrs) * 100
    SPE <- TN/sum(dat.all$es > es-5e-03 & dat.all$es < es+5e-03 & dat.all$hours == dat.all$cont_hrs) * 100
    ACC <- (TP + TN)/N * 100
    dat.pow <- rbind(dat.pow, c(es,N,TP,FP,FPR,TN,FN,POW,SEN,SPE,ACC))
    }
}

dat.pow <- data.frame(dat.pow)
colnames(dat.pow) <- c("es","N","TP","FP","FPR","TN","FN","pow","sen","spe","acc")

pd.lid <- ggplot(dat.pow) +
  geom_line(aes(x = es, y = pow),col = '#D32F2F', linetype = 'twodash', lwd = 2) +
  scale_x_continuous(breaks = seq(-2,2,0.05),
                     minor_breaks = seq(-2,2,0.01)) +
  scale_y_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,100,5)) +
  labs(title = "Power Dynamics",
       subtitle = "Using Cubic Spline @ 5% FPR",
       x = "Effect Size",
       y = "Power") +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.position = c(0.8,0.2),
        legend.background = element_blank(),
        legend.text = element_text(size=15),
        legend.title =  element_text(size=15),
        plot.title = element_text(size=25,hjust = 0),
        plot.subtitle = element_text(size=20,hjust = 0)) +
  coord_cartesian(xlim = c(0.8,1.2),
                  ylim = c(0,100))
ggsave(sprintf("%spower_cs.png",out_path),
       pd.lid,
       width = 10,height = 10)

ggplot(hello) +
  geom_line(aes(x = es, y = pow, col = norm), linetype = 'twodash', lwd = 1.3) +
  scale_x_continuous(breaks = seq(-2,2,0.05),
                     minor_breaks = seq(-2,2,0.01)) +
  scale_y_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,100,5)) +
  labs(title = "Power Dynamics",
       subtitle = "@ 5% FPR",
       x = "Effect Size",
       y = "Power") +
  scale_color_manual(name = "Normalization",
                     breaks = c("CS","LID"),
                     values = c("CS"="#512DA8","LID"="#C2185B"),
                     labels = c("Cubic Spline","LI Detector")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.position = c(0.85,0.2),
        legend.background = element_rect(color = 'grey60'),
        legend.text = element_text(size=15),
        legend.title =  element_text(size=15),
        plot.title = element_text(size=25,hjust = 0),
        plot.subtitle = element_text(size=20,hjust = 0)) +
  coord_cartesian(xlim = c(0.8,1.2),
                  ylim = c(0,100))
ggsave(sprintf("%spower_cslid.png",out_path),
       width = 10,height = 10)

