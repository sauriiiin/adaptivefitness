##### JIAZHEN TUTORIAL
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 01/29/2021

##### INSTALL PACKAGES
# comment this section out once you have installed them.
# you do not need to run this section everytime.
install.packages("RMariaDB")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("ggExtra")
install.packages("grid")
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("stringr")

##### INITIALIZE
library(RMariaDB)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)

sql.usr <- 'sbp29'
sql.pwd <- 'Ku5hani@28'
sql.db <- 'Branden'

conn <- dbConnect(MariaDB(), dbname = sql.db,
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')

##### VIEW TABLES PRESENT IN SQL
## all tables
all.tables <- dbGetQuery(conn, 'show tables')
View(all.tables)

## find tables with a particular name
name.like <- 'respiration_plusMet_Glu'
all.tables[str_detect(all.tables$Tables_in_Branden, name.like),]

##### GETTING DATA FROM THAT TABLE
## Get RAW Colony Sizes
tblname_JPEG <- 'respiration_plusMet_Glu_384_JPEG'
raw.data <- dbGetQuery(conn, sprintf('select * from %s', tblname_JPEG))
View(raw.data)

## Getting RAW Data with Coordinates and Mutant Names
# name.like <- 'espiration'
# all.tables[str_detect(all.tables$Tables_in_Branden, name.like),]
tblname_p2o <- 'Respiration_pos2orf_name'
tblname_p2c <- 'Respiration_pos2coor'

cs.data <- dbGetQuery(conn, sprintf('select c.*, b.orf_name, a.hours, a.average
                                 from %s a, %s b, %s c
                                 where a.pos = b.pos and a.pos = c.pos
                                 order by a.hours, c.plate, c.col, c.row',
                                 tblname_JPEG, tblname_p2o, tblname_p2c))

##### PLOTTING DATA
ggplot(cs.data) +
  geom_point(aes(x = hours, y = average, col = orf_name)) +
  geom_smooth(aes(x = hours, y = average, col = orf_name),
              method = 'loess') +
  scale_x_continuous(breaks = seq(0,80,4)) +
  theme_linedraw() 
# ggsave("growth_plot.jpg",
#        width = 7, height = 7,
#        dpi = 300)

### EXERCISE
# Change x-axis title from 'hours' to 'Time (hrs)'
# Change x-axis so that ticks go from 0 to 170 at intervals of 10
# Change y-axis title from 'average' to 'Colony Size (pix.)'
# Change legend title from 'orf_name' to 'Mutants'

