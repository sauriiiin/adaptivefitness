##### ROC 4 Control Experiment
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/07/2019

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("ggplot2")
# library(RMariaDB)
library(ggplot2)
# source("R/functions/initialize.sql.R")
# conn <- initialize.sql("saurin_test")

data <- read.csv(file="rawdata/4C3_GA1_TEMP_P_ES_17.csv", header=TRUE, sep=",")

out_path = 'figs/';
expt_name = '4C3_GA1'
hr = 17
alldat <- data.frame()

#####
# query <- 'select a.orf_name, a.hours, b.p 
# from 4C3_GA1_TEMP_6144_FITNESS_STATS a, 4C3_GA1_TEMP_6144_PVALUE b 
# where a.orf_name = b.orf_name and a.hours = b.hours 
# order by p asc'
# query <- 'select a.orf_name, a.hours, a.cs_mean - 1.0009 es, b.p
# from 4C3_GA1_TEMP_6144_FITNESS_STATS a, 4C3_GA1_TEMP_6144_PVALUE b
# where a.orf_name = b.orf_name and a.hours = b.hours
# and abs(a.cs_mean - 1.0009) <= .5
# order by  p asc'
# data <- dbGetQuery(conn, query)
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

ggplot() +
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
  labs(title = "D. ROC Curve",
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
        plot.title = element_text(size=25,hjust = -0.115))
ggsave(sprintf("%s%s_ROC_%d.png",
               out_path,expt_name,hr),
       width = 10,height = 10.5)










