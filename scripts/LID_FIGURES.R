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
out_path = 'figs/';
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
ggsave(sprintf("%s%s_FIG4.png",
               out_path,expt_name),
       fig4,
       width = 20,height = 20)