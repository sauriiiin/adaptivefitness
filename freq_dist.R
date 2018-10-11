##### FREQ_DIST
##### proto-gene frequency distribution plots
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/11/2018

library(RMariaDB)
source("R/functions/initialize.sql.R")

mydb <- initialize.sql("saurin_test")

dbListFields(mydb, '1536_PS2_previous')


dbDisconnect(mydb)
