##### ALL F28 DATA
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 04/07/2021

##### INITIALIZE
library(ggplot2)
library(ggridges)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(scales)
library(egg)
library(zoo)
library(ggrepel)
library(reshape2)
library(RMariaDB)

source("R/functions/empirical_p.R")
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/f28fu/repeat/'

##### LOAD DATA
load("/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/solid_results.RData")
load("/home/sbp29/R/Projects/adaptivefitness/figs/f28fu/repeat/f28fu_data.RData")

temp1 <- solid_fit[(solid_fit$density == 1536 & solid_fit$hours == 48 & solid_fit$arm != 'SDA') |
                       (solid_fit$density == 1536 & solid_fit$hours == 72 & solid_fit$arm == 'SDA') |
                       (solid_fit$density == 6144 & solid_fit$hours == 24 & solid_fit$arm != 'SDA') |
                       (solid_fit$density == 6144 & solid_fit$hours == 36 & solid_fit$arm == 'SDA'),]
head(temp1)
colnames(temp1)
temp1 <- temp1[,c(1,2,7,8,9,10,18,17)]
temp1$stage <- 'FS'
temp1$location <- "Pittsburgh"
temp1$plasmid <- "Yes"

temp2 <- foa.data[((foa.data$density == 1536 & foa.data$hours == 60) |
                     (foa.data$density == 6144 & foa.data$hours == 36 & foa.data$condition %in% c('GAL','SDA','CAS')) |
                     (foa.data$density == 6144 & foa.data$hours == 21 & foa.data$condition == 'GLU')) &
                      !(foa.data$orf_name %in% c('BOR','REF','YHR021W-A')),]
# temp2$condition <- factor(temp2$condition, levels = c('GLU','GAL','CAS','SDA'))
temp2 <- temp2[,c(1,2,6,7,14,15,10,11)]
temp2$stage <- 'FS'
temp2$location <- "Pittsburgh"
temp2$plasmid <- "No"

colnames(temp2) <- colnames(temp1)
head(temp2)


head(foa.sd.data)
foa.sd.data$replicate[foa.sd.data$exp_id == 997] <- "R1"
foa.sd.data$replicate[foa.sd.data$exp_id == 998] <- "R2"

foa.sd.data$density <- 6144
foa.sd.data$arm <- 'GAL'
foa.sd.data$stage <- 'FS'
foa.sd.data$location <- 'SanDiego'
foa.sd.data$plasmid <- 'No'

head(foa.sd.data)
temp3 <- foa.sd.data[,c(2,9,1,4,6,7,8,10,11,12,13)]
colnames(temp3) <- colnames(temp1)
head(temp3)

#####
##### F28FUR PS2 & S3
tblname_dat1 <- 'F28FUR_PS2_1536_FITNESS'
tblname_dat2 <- 'F28FUR_S3_384_FITNESS'
tblname_p2o <- 'F28FUR_pos2orf_name'
tblname_p2c <- 'F28FUR_pos2coor'

foa.ps2 <- dbGetQuery(conn, sprintf('select c.*, b.orf_name, a.hours, a.average, a.fitness
                                    from %s a, %s b, %s c
                                    where a.pos = b.pos and a.pos = c.pos
                                    order by density, hours, plate, col, row',
                                    tblname_dat1, tblname_p2o, tblname_p2c))
head(foa.ps2)

foa.ps2$arm[foa.ps2$plate == 1 & foa.ps2$hours != 30] <- 'GLU'
foa.ps2$arm[foa.ps2$plate == 2 & foa.ps2$hours != 30] <- 'GAL'
foa.ps2$arm[foa.ps2$plate == 3 & foa.ps2$hours != 30] <- 'CAS'
foa.ps2$arm[foa.ps2$plate == 4 & foa.ps2$hours != 30] <- 'SDA'
foa.ps2$arm[foa.ps2$plate == 1 & foa.ps2$hours == 30] <- 'GAL'
foa.ps2$arm[foa.ps2$plate == 2 & foa.ps2$hours == 30] <- 'CAS'
foa.ps2$arm[foa.ps2$plate == 3 & foa.ps2$hours == 30] <- 'SDA'
foa.ps2$arm[foa.ps2$plate == 4 & foa.ps2$hours == 30] <- 'GLU'

foa.ps2$replicate[foa.ps2$arm == 'GLU'] <- 'R1'
foa.ps2$replicate[foa.ps2$arm == 'GAL'] <- 'R1'
foa.ps2$replicate[foa.ps2$arm == 'CAS'] <- 'R2'
foa.ps2$replicate[foa.ps2$arm == 'SDA'] <- 'R3'


foa.ps2$stage <- 'PS2'
foa.ps2$location <- 'Pittsburgh'
foa.ps2$plasmid <- 'No'

foa.ps2 <- foa.ps2[,c(1,2,6:9,11,10,12:14)]
head(foa.ps2)


foa.s3 <- dbGetQuery(conn, sprintf('select c.*, b.orf_name, a.hours, a.average, a.fitness
                                    from %s a, %s b, %s c
                                    where a.pos = b.pos and a.pos = c.pos
                                    order by density, hours, plate, col, row',
                                    tblname_dat2, tblname_p2o, tblname_p2c))
head(foa.s3)
foa.s3$replicate <- 'R1'
foa.s3$arm <- '5FOA'
foa.s3$stage <- 'S3'
foa.s3$location <- 'Pittsburgh'
foa.s3$plasmid <- 'No'

foa.s3 <- foa.s3[,c(1,2,6:14)]
head(foa.s3)


f28.dat <- rbind(temp1, foa.s3, foa.ps2, temp2, temp3)

save(f28.dat, file = sprintf('%sf28alldata.RData',out_path))



