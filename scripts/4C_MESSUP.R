##### 4C_MESSUP
##### MESSING UP THE FINAL PLATES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 06/26/2019

##### INITIALIZE
library(RMariaDB)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
alldat = dbGetQuery(conn, 'select b.*, a.orf_name, a.hours, a.bg, a.average, a.fitness
              from 4C3_GA1_6144_FITNESS a, 4C3_pos2coor6144 b
              where a.pos = b.pos
              order by hours, 6144plate, 6144col, 6144row')
alldat$rnd_hrs <- alldat$hours

for (i in 1:dim(alldat)[1]) {
  if (alldat$orf_name[i] != "BF_control" & !is.na(alldat$orf_name[i])) {
    temp <- sample_n(alldat[alldat$pos == alldat$pos[i],], 1)
    alldat$hours[i] <- temp$rnd_hrs
    alldat$average[i] <- temp$average
  }
}

dbWriteTable(conn, "4C3_GA1_RND_6144_FITNESS", alldat, overwrite = T)


