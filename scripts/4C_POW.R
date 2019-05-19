##### 4C ES VS POWER
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/18/2019

##### INITIALIZE
library(readxl)
library(ggplot2)

d <- read_excel("rawdata/cdata_ga1_clean.xlsx",col_types = "numeric")
d.smooth <- read_excel("rawdata/cdata_ga1_clean_smooth.xlsx",col_types = "numeric")

##### SET VARIABLES
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/';
density = 6144;

# d.smooth <- data.frame(spline(d$es,d$power,n = 5*length(d$es),method='fmm'))

##### PLOT

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
ggsave(sprintf("%s%s_POW_FIG4C.png",
               out_path,expt_name),
       width = 10,height = 10.5)
