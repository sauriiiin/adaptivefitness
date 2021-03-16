library(ggplot2)
library(ggridges)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(RMariaDB)

source("R/functions/isoutlier.R")
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

data = dbGetQuery(conn, 'select * from trial_384_RAW a, F28FUR_pos2coor b
                  where a.pos = b.pos
                  order by hours, plate, col, row')


ggplot(data,
       aes(x = hours, y = pixels)) +
  # geom_point(aes(col = pos)) +
  geom_boxplot(aes(group = hours), outlier.shape = NA) +
  facet_grid(.~plate)


hi <- data[data$hours == 75,]
