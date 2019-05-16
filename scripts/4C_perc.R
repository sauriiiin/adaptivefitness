
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/16/2019

##### INITIALIZE
library(ggplot2)

##### GET/SET DATA
out_path = 'figs/';
expt_name = '4C3_GA1'

data <- read.csv(file="rawdata/perc.csv", header=TRUE, sep=",")
data$pl[data$i < 5] = "One"; data$pl[data$i > 4] = "Two"
data$source[data$i == 1 | data$i == 5] = 'TL'
data$source[data$i == 2 | data$i == 6] = 'TR'
data$source[data$i == 3 | data$i == 7] = 'BL'
data$source[data$i == 4 | data$i == 8] = 'BR'

data$a_diff <- data$avg_ul - data$avg_ll
data$f_diff <- data$f_ul - data$f_ll
data$m_diff <- data$mck_ul - data$mck_ll
data$mf_diff <- data$m_ul - data$m_ll

##### PERC PLOT
ggplot(data, aes(x = data$avg_sd, y = data$a_diff, col = pl, fill = source),size = 4) +
  geom_point(shape = 21, lwd = 2) +
  labs(title = 'Pixel Count',
       x = 'Standard Deviation',
       y = 'Upper Range - Lower Range') +
  scale_colour_manual(name="Plate",
                      values=c("One"="#212121","Two"="#FFA000")) +
  scale_fill_manual(name="Source",
                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
                      breaks=c("TL","TR","BL","BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_x_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,100,5),
                     limits = c(10,90)) +
  scale_y_continuous(breaks = seq(-100,100,25),
                     minor_breaks = seq(-100,100,12.5)) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=15),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        legend.position = 'bottom',
        plot.title = element_text(size=15,hjust = 0.5))
ggsave(sprintf("%s%s_PERC_PIX.png",
               out_path,expt_name),
       width = 10,height = 10.5)


ggplot(data, aes(x = source, y = mf_diff, col = pl, fill = source, size = mock_sd)) +
  geom_point(shape = 21) +
  labs(title = 'Mock Fitness',
       x = 'Source Plate',
       y = 'Upper Range - Lower Range') +
  scale_size(name = 'Pix Range Diff',range = c(2, 6)) +
  scale_colour_manual(name="Plate",
                      values=c("One"="#212121","Two"="#FFA000")) +
  scale_fill_manual(name="Source",
                    values=alpha(c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),0.2),
                    breaks=c("TL","TR","BL","BR"),
                    labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  # scale_x_continuous(breaks = seq(0,100,10),
  #                    minor_breaks = seq(0,100,5),
  #                    limits = c(10,90)) +
  # scale_y_continuous(breaks = seq(-100,100,25),
  #                    minor_breaks = seq(-100,100,12.5)) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=15),
        legend.text = element_text(size=10),
        legend.title = element_text(size=15),
        legend.position = 'right',
        plot.title = element_text(size=15,hjust = 0.5))
ggsave(sprintf("%s%s_PERC_sourceF.png",
               out_path,expt_name),
       width = 11,height = 10)
