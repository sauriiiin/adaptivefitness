
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/16/2019

##### INITIALIZE
library(ggplot2)

##### GET/SET DATA
data <- read.csv(file="rawdata/perc.csv", header=TRUE, sep=",")
data$pl[data$i < 5] = 1; data$pl[data$i > 4] = 2
data$source[data$i == 1 | data$i == 5] = 'TL'
data$source[data$i == 2 | data$i == 6] = 'TR'
data$source[data$i == 3 | data$i == 7] = 'BL'
data$source[data$i == 4 | data$i == 8] = 'BR'

data$a_diff <- data$avg_ul - data$avg_ll
data$f_diff <- data$f_ul - data$f_ll
data$m_diff <- data$mck_ul - data$mck_ll
data$mf_diff <- data$m_ul - data$m_ll
