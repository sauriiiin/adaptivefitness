##### LID PAPER FIGURES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/18/2019

##### INITIALIZE
library(readxl)
library(ggplot2)
library(tidyverse)
library(egg)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### INITIALIZE
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/lid_paper/';
density = 6144;

tablename_pval = sprintf('%s_%d_PVALUE',expt_name,density)
hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_pval))
pvals = seq(0,1,0.01)

##### FIGURE 4
data = dbGetQuery(conn, sprintf('select * from %s',tablename_pval))
data = data[data$orf_name != 'BFC100',]
pdata = data.frame()
i = 1

for (hr in unique(data$hours)) {
  temp = data[data$hours == hr,]
  n = dim(temp)[1]
  for (p in pvals) {
    pdata = rbind(pdata, c(hr, p, sum(temp$p <= p)/n))
  }
  if (i == 1) {
    colnames(pdata) = c('hours','p','fpr')
    g = ggplot() + 
      geom_line(data = pdata[pdata$hours == hr,], aes(x = p, y = fpr, col = as.character(hours)), lwd = 1.5)
  } else {
    g = g + geom_line(data = pdata[pdata$hours == hr,], aes(x = p, y = fpr, col = as.character(hours)), lwd = 1.5)
  }
  i = i + 1
}

##### A
fpr <- g + geom_line(data = pdata, aes(x = p, y = p, col = 'red'),
                     linetype = 'dashed', lwd = 1.1, alpha = 0.7) +
  labs(title = "A. False positive rate",
       x = "P-value cut-off",
       y = "False Positive Rate") +
  scale_x_continuous(breaks = seq(0,1,0.1),
                     minor_breaks = seq(0,1,0.05),
                     limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,0.1),
                     minor_breaks = seq(0,1,0.05),
                     limits = c(0,1)) +
  scale_color_manual(name = "Hours",
                     breaks=c("8","10","14","16","18"),
                     values=c("8"="#D32F2F","10"="#536DFE","14"="#388E3C","16"="#795548","18"="#00BCD4","red"="red",
                              "0"="transparent","9"="transparent","11"="transparent","13"="transparent","17"="transparent")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        # axis.title = element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = "bottom",
        plot.title = element_text(size=25,hjust = -0.1),
        plot.subtitle = element_text(size=20,hjust = -0.1)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

##### B
fpr.zoom <- g + geom_line(data = pdata, aes(x = p, y = p, col = 'red'),
                          linetype = 'dashed', lwd = 1.1, alpha = 0.7) +
  labs(title = "B. False positive rate for p < 0.1",
       x = "P-value cut-off",
       y = "False Positive Rate") +
  scale_x_continuous(breaks = seq(0,1,0.01),
                     minor_breaks = seq(0,1,0.005),
                     limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,0.01),
                     minor_breaks = seq(0,1,0.005),
                     limits = c(0,1)) +
  scale_color_manual(name = "Hours",
                     breaks=c("8","10","14","16","18"),
                     values=c("8"="#D32F2F","10"="#536DFE","14"="#388E3C","16"="#795548","18"="#00BCD4","red"="red",
                              "0"="transparent","9"="transparent","11"="transparent","13"="transparent","17"="transparent")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        # axis.title = element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = "bottom",
        plot.title = element_text(size=25,hjust = -0.16),
        plot.subtitle = element_text(size=20,hjust = -0.1)) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 0.1))

##### C
d <- read_excel("rawdata/cdata_ga1_clean.xlsx",col_types = "numeric")
d.smooth <- read_excel("rawdata/cdata_ga1_clean_smooth.xlsx",col_types = "numeric")

pow <- ggplot() +
  geom_line(data=d.smooth,aes(x=es, y=power, col='Smoothed Fit'),
            lwd = 1.2) +
  geom_point(data=d,aes(x=es, y=power, col='Data Point'),
             size = 4) +
  labs(title = "C. Change in power with effect size",
       x = "Effect Size",
       y = "Power") +
  scale_x_continuous(breaks = seq(0,10,0.05),
                     minor_breaks = seq(0,10,0.025)) +
  scale_y_continuous(breaks = seq(-20,140,10),
                     minor_breaks = seq(-20,140,5)) +
  scale_color_manual(name = "Result",
                     breaks=c("Data Point", "Smoothed Fit"),
                     values=c("Data Point"="#303F9F",
                              "Smoothed Fit"="#FF4081")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        # axis.title = element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = "bottom",
        plot.title = element_text(size=25,hjust = -0.17),
        plot.subtitle = element_text(size=20,hjust = -0.1)) +
  coord_cartesian(xlim = c(0.8, 1.2), ylim = c(0, 100))

##### D
data <- read.csv(file="rawdata/4C3_GA1_TEMP_P_ES_17.csv", header=TRUE, sep=",")
hr = 17
alldat <- data.frame()

data <- data[data$orf_name != 'BFC100',]
data$es <- abs(data$es)

for (es in c(0.05,0.1,0.15,0.2,0.5)) {
  data.temp <- data[data$es <= es,]
  fpdat <- data.temp[data.temp$hours == 17,]
  tpdat <- data.temp[data.temp$hours != 17,]
  fp <- NULL
  tp <- NULL
  for (i in 1:length(unique(data$p))) {
    p <- unique(data$p)[i]
    fp[i] <- dim(fpdat[fpdat$p <= p,])[1]/dim(fpdat)[1] * 100
    tp[i] <- dim(tpdat[tpdat$p <= p,])[1]/dim(tpdat)[1] * 100
  }
  plotdat <- data.frame(unique(data$p),fp,tp)
  colnames(plotdat) <- c('p','FalsePositive','TruePositive')
  plotdat$ES <- es*100
  
  alldat <- rbind(alldat,plotdat)
}

lw = 1.2

roc <- ggplot() +
  geom_line(data = alldat[alldat$ES==50,],
            aes(x = FalsePositive, y = p*100, col = 'p'),
            linetype = "dotdash", lwd = lw) +
  geom_line(data = alldat[alldat$ES==5,],
            aes(x = FalsePositive, y = TruePositive, col = '5%'),
            lwd = lw) +
  geom_line(data = alldat[alldat$ES==10,],
            aes(x = FalsePositive, y = TruePositive, col = '10%'),
            lwd = lw) +
  geom_line(data = alldat[alldat$ES==15,],
            aes(x = FalsePositive, y = TruePositive, col = '15%'),
            lwd = lw) +
  geom_line(data = alldat[alldat$ES==20,],
            aes(x = FalsePositive, y = TruePositive, col = '20%'),
            lwd = lw) +
  geom_line(data = alldat[alldat$ES==50,],
            aes(x = FalsePositive, y = TruePositive, col = '50%'),
            lwd = lw) +
  labs(title = "D. ROC curve",
       x = "False Positive Rate",
       y = "Power") +
  scale_x_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,100,5),
                     limits = c(0,100)) +
  scale_y_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,100,5),
                     limits = c(0,100),
                     sec.axis = sec_axis(~./100,
                                         breaks = seq(0,1,0.1),
                                         name='p-value')) +
  scale_colour_manual(name="Effect Size Threshold",
                      breaks=c("5%","10%","15%","20%","50%"),
                      values=c("5%"="#D32F2F","10%"="#536DFE","15%"="#388E3C",
                               "20%"="#795548","50%"="#00BCD4","p"="#FFA000")) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        axis.line.y.right = element_line(color = "#FFA000"),
        axis.text.y.right = element_text(size=15, color = "#FFA000"),
        axis.title.y.right = element_text(size=20, color = "#FFA000"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = "bottom",
        plot.title = element_text(size=25,hjust = -0.1))

##### FINAL FIG4
fig4 <- ggarrange(fpr, fpr.zoom, pow, roc,
                  nrow = 2)
ggsave(sprintf("%sfigure4.png",out_path),
       fig4,
       width = 20,height = 20)


##### FIGURE 3
tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_nfit = sprintf('%s_NIL_%d_FITNESS',substr(expt_name,1,6),density);
tablename_p2o = '4C3_pos2orf_name1';
tablename_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

p2c = dbGetQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
                               p2c_info[1],
                               p2c_info[2],
                               p2c_info[3],
                               p2c_info[4]))

n_plates = dbGetQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                    p2c_info[2],
                                    p2c_info[1],
                                    p2c_info[2]))
hr = hours[[1]][9]
pl = n_plates[[1]][1]

fitdat = dbGetQuery(conn, sprintf('select c.*, a.orf_name, a.hours, a.bg, a.average, a.fitness,
                                  b.bg nbg, b.average naverage, b.fitness nfitness
                                  from %s a, %s b, %s c
                                  where a.hours = %d and a.hours = b.hours
                                  and a.pos = b.pos and b.pos = c.pos
                                  and c.%s = %d order by c.%s, c.%s',
                                  tablename_fit,tablename_nfit,
                                  p2c_info[1],hr,p2c_info[2],
                                  pl,p2c_info[3],p2c_info[4]))
fitdat$bg[is.na(fitdat$average)] = NA
min = min(fitdat$average, na.rm=T)
max = max(fitdat$average, na.rm=T)

fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==1] = 'TL'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==1] = 'BL'
fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==0] = 'TR'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==0] = 'BR'

fitdat$colony[fitdat$orf_name == 'BF_control'] = 'Reference'
fitdat$colony[fitdat$orf_name != 'BF_control'] = 'Query'
fitdat$colony[is.na(fitdat$orf_name)] = 'Gap'

##### A
obs <- ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = average,
                 shape = colony,
                 alpha = colony),size = 7,na.rm = T) +
  labs(title = "A. Observed Colony Size (Pixel Count)",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1),limits = c(5,92)) +
  scale_y_continuous(breaks = seq(1,64,1),limits = c(60,5),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(min,max),
                        palette = "Set1") +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20,face="bold"),
        legend.position = "right",
        plot.title = element_text(size=25,hjust = 0),
        plot.subtitle = element_text(size=20,hjust = 0))

##### B
pre <- ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = bg,
                 shape = colony,
                 alpha = colony),size = 7,na.rm = T) +
  labs(title = "B. Predicted Colony Size (Pixel Count)",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1),limits = c(5,92)) +
  scale_y_continuous(breaks = seq(1,64,1),limits = c(60,5),trans = 'reverse') +
  scale_color_distiller(name = "PIX",
                        limits = c(min,max),
                        palette = "Set1") +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20,face="bold"),
        legend.position = "right",
        plot.title = element_text(size=25,hjust = 0),
        plot.subtitle = element_text(size=20,hjust = 0))

##### C
fit <- ggplot(data = fitdat, aes(x = `6144col`, y = `6144row`)) +
  geom_point(aes(x = `6144col`, y = `6144row`,col = fitness,
                 shape = colony,
                 alpha = colony),size = 7,na.rm = T) +
  labs(title = "C. Fitness Landscape",
       x = "",
       y = "") +
  scale_x_continuous(breaks = seq(1,96,1),limits = c(5,92)) +
  scale_y_continuous(breaks = seq(1,64,1),limits = c(60,5),trans = 'reverse') +
  scale_color_distiller(name = "FIT",
                        limits = c(0.7,1.3),
                        breaks = c(0.70,0.85,1.00,1.15,1.30),
                        palette = "Set1") +
  scale_shape_manual(name="Colony Kind",
                     values=c("Gap"=15,"Query"=15,"Reference"=15),
                     breaks=c("Reference","Query","Gap"),guide=F) +
  scale_alpha_manual(values=c("Gap"=1,"Query"=1,"Reference"=1),
                     guide=F) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20,face="bold"),
        legend.position = "right",
        plot.title = element_text(size=25,hjust = 0),
        plot.subtitle = element_text(size=20,hjust = 0))

##### D
raw <- ggplot(data = fitdat, aes(x=average, col = source)) +
  geom_density(lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                      breaks=c("TL","TR","BL","BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_x_continuous(breaks = seq(0,1000,50),
                     minor_breaks = seq(0,1000,25),
                     limits = c(min,max)) +
  scale_y_continuous(breaks = seq(0,1,0.002),
                     minor_breaks = seq(0,1,0.001),
                     limits = c(0,0.02)) +
  labs(title = 'D. Raw colony sizes',
       x = 'Observed Pixel Count', y = 'Density') +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = 'bottom',
        legend.background = element_rect(fill="lightblue", 
                                         size=0.5, linetype="solid"),
        plot.title = element_text(size=25,hjust = -0.15))

##### E
src.nrm <- ggplot(data = fitdat, aes(x=fitness, col = source)) +
  geom_density(lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                      breaks=c("TL","TR","BL","BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_x_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.025),
                     limits = c(0.7,1.3)) +
  scale_y_continuous(breaks = seq(0,15,1),
                     minor_breaks = seq(0,15,0.5),
                     limits = c(0,12)) +
  labs(title= 'E. With Source Normalization',
       x = 'Fitness', y = 'Density') +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = 'bottom',
        legend.background = element_rect(fill="lightblue", 
                                         size=0.5, linetype="solid"),
        plot.title = element_text(size=25,hjust = -0.13))

##### F
no.src.nrm <- ggplot(data = fitdat, aes(x=nfitness, col = source)) +
  geom_density(lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                      breaks=c("TL","TR","BL","BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_x_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.025),
                     limits = c(0.7,1.3)) +
  scale_y_continuous(breaks = seq(0,15,1),
                     minor_breaks = seq(0,15,0.5),
                     limits = c(0,12)) +
  labs(title= 'F. Without Source Normalization',
       x = 'Fitness', y = 'Density') +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.position = 'bottom',
        legend.background = element_rect(fill="lightblue", 
                                         size=0.5, linetype="solid"),
        plot.title = element_text(size=25,hjust = -0.13))

##### FINAL FIG 3
fig3.top <- ggarrange(obs, pre, fit,
                  nrow = 1)
ggsave(sprintf("%sfigure3_top.png",out_path),
       fig3.top,
       width = 30,height = 6.5)

fig3.bot <- ggarrange(raw, src.nrm, no.src.nrm,
                      nrow = 1)
ggsave(sprintf("%sfigure3_bot.png",out_path),
       fig3.bot,
       width = 30,height = 11)








