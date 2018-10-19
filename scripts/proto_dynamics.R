##### PROTO-DYNAMIC
##### proto-gene frequency distribution plots
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/19/2018

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("ggplot2")
library(RMariaDB)
library(ggplot2)
library(reshape2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### FETCH DATA
query = dbSendQuery(conn, ("select a.orf_name, a.cs_median cs_ps1,
b.cs_median cs_ps2,
c.cs_median cs_fs
from PT2_PGLU_PS1_6144_RES_eFDR a,
PT2_PGLU_PS2_6144_RES_eFDR b,
PT2_PGLU_FS_6144_RES_eFDR c
where a.protogene = 1
and a.hours = 12 and b.hours = 12 and c.hours = 12
and a.orf_name = b.orf_name and b.orf_name = c.orf_name"))
data = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

##### PLOTS
trial <- melt(data, id.vars = "orf_name")
ggplot(trial, aes(variable, value, group=factor(orf_name))) +
  geom_line(aes(color=factor(orf_name))) +
  theme(legend.position="none")


