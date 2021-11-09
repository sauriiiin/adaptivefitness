library(dplyr)
library(ggplot2)
library(RMariaDB)

source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")


plate <- dbGetQuery(conn, 'select * from 384DENSITYCOOR')
plate$destination_position <- (plate$`384col`-1)*16 + plate$`384row`

## BORDERS
plate$strain_id[plate$`384row` == min(plate$`384row`) |
                  plate$`384row` == max(plate$`384row`) |
                  plate$`384col` == min(plate$`384col`) |
                  plate$`384col` == max(plate$`384col`)] <- -2

## GAPS
n_gaps <- 2
plate$strain_id[plate$destination_position %in%
  sample(plate$destination_position[is.na(plate$strain_id)],n_gaps)] <- 0

## STRAINS
n_strains <- 18
n_reps <- 17
for (i in seq(1,n_strains,1)) {
  plate$strain_id[plate$destination_position %in%
                    sample(plate$destination_position[is.na(plate$strain_id)],n_reps)] <- i
}
plate$destination_position[is.na(plate$strain_id)]

plate <- plate[order(plate$destination_position),]

## SAVING FILE
write.csv(plate, file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/condition_control_plate.csv')

## UPLOADING
plate.map <- read.csv('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/condition_control_map.csv')
dbWriteTable(conn, 'TRANS_CONTROL_MAP', plate.map)

