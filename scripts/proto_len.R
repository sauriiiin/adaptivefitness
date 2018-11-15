##### PROTO_LEN
##### proto-gene length distribution in our collection
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 15/11/2018

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("ggplot2")
library(RMariaDB)
library(ggplot2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### FETCH DATA

query = dbSendQuery(conn, ("select a.orf_name, length(b.seq_nt_atg) nt_len
from BARFLEX_SPACE_AGAR a, ORFS_SEQUENCES b
where a.strain_id > 0 and a.orf_name = b.orf_name"))
data = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

query = dbSendQuery(conn, ("select a.orf_name, length(b.seq_nt_atg) nt_len
from BACTERIAL_SPACE_MAP a, ORFS_SEQUENCES b
where a.pos > 0 and a.orf_name = b.orf_name"))
data2 = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}




