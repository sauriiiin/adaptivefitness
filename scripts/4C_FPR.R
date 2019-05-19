##### 4C FALSE POSITIVES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/10/2019

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("ggplot2")
library(RMariaDB)
library(ggplot2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/';
density = 6144;

tablename_pval = sprintf('%s_%d_PVALUE',expt_name,density)
hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_pval))
pvals = seq(0,1,0.01)

##### GET DATA
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
ggsave(sprintf("%s%s_FPR.png",
               out_path,expt_name),
       width = 10,height = 10.5)


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
ggsave(sprintf("%s%s_FPR_ZOOM.png",
               out_path,expt_name),
       width = 10,height = 10.5)



