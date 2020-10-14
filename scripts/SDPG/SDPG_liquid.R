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
df <- as.data.frame(readxl::read_excel("rawdata/SDPG/PR_Results.xlsx"))
map <- as.data.frame(readxl::read_excel("rawdata/SDPG/SDPG_96well.xlsx"))
# colnames(map) <- c('')

##### INITIALIZE OUTPUT
span_parameter <- 0.15
outlier_threshold <- 1

out       = NULL
maxgr     = NULL
dtime     = NULL
ltime     = NULL
sod       = NULL
arm       = NULL
sample    = NULL
media     = NULL
replicate = NULL

ii <- 1
min.od <- 0.001
min.t <- 0
for (a in unique(df$Arm)) {
  temp <- df[df$Arm == a,]
  temp$Time <- seq(15,15*dim(temp)[1],15)
  for (i in 4:dim(temp)[2]) {
    if (sum(temp[[i]] > min.od, na.rm = T) > 30) {
      # temp[[i]] <- rollapplyr(temp[[i]], 5, mean, na.rm = TRUE, fill = 0, align = 'right', partial = T)
      # fit0 <- fit_easylinear(temp$Time[(temp[[i]] > min.od) & (temp$Time > min.t) & !is.na(temp[[i]])],
      #                        temp[(temp[[i]] > min.od) & (temp$Time > min.t) & !is.na(temp[[i]]),i], h=40, quota = 1);
      lo <- loess.smooth(temp$Time, log(temp[,i]),
                         span = span_parameter)
      fit0 <- fit_easylinear(lo$x, exp(lo$y), h=20, quota = 1);
      # fit0 <- fit_easylinear(temp$Time, temp[,i], h=7, quota = 0.95);
      arm[ii] = a
      maxgr[ii] = coef(fit0)[[3]]
      dtime[ii] = log(2)/coef(fit0)[[3]] #* 60
      ltime[ii] = coef(fit0)[[4]]
      sod[ii] = temp[[9,i]]
      sample[ii] = as.character(map$orf_name[map$well == colnames(temp[i])])
      replicate[ii] = map$replicate[map$well == colnames(temp[i])]
      media[ii] = temp$Media[1]
      ii <- ii + 1

      jpeg(sprintf('%s%s_%s_%s-%s_GC.png',
                   path.out.gc,
                   expt_name,
                   a,
                   map$orf_name[map$well == colnames(temp[i])],
                   as.character(map$replicate[map$well == colnames(temp[i])])),
      width=600, height=600)
      plot(fit0, log = 'y',
           main=sprintf('%s | %s | %s \n Doubling Time = %0.2f mins',
                        a, colnames(df[i]), map$orf_name[map$well == colnames(temp[i])],
                        log(2)/coef(fit0)[[3]]),
           ylim = c(0.01,6))
      dev.off()
    }
  }
  #df[df$Arm == a,] <- temp
}

out <- data.frame(Arm = arm,
                  Sample = sample,
                  Media = media,
                  Replicate = replicate,
                  MaxGR = maxgr,
                  DTime = dtime,
                  SOD = sod,
                  LagTime = ltime)

# write.csv(out, file = sprintf("%s%s_GC_PLC_Result.csv",path.out,expt_name))

##### DOUBLING TIME RESULTS
ggplot(out[out$Sample != 'YHR021W-A' & out$Arm != 'INIT_96',],
       aes(x = Sample, y = DTime)) +
  geom_boxplot() + geom_point() +
  coord_flip(ylim = c(60,120)) +
  labs(y = 'Doubling Time (mins)') +
  facet_wrap(~Arm)

ggsave(sprintf('%s%s_DTime_Box.png',path.out,expt_name),
       height = two.c, width = two.c, units = 'mm',
       dpi = 300)

##### GROWTH CRUVES

