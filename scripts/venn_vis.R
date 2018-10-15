##### VENN_VIS
##### venn diagram visualizations of proto data
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/13/2018

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("VennDiagram")
library(RMariaDB)
library(VennDiagram)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### FETCH DATA
query = dbSendQuery(conn, ("select distinct orf_name orfs
from PT_SA_CN_6144_FITNESS_STATS
where orf_name in
(select orf_name from PROTOGENES
where selected + longer + translated < 3)"))
proto = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

query = dbSendQuery(conn, ("select distinct orf_name orfs
from PT_SA_CN_6144_FITNESS_STATS
where orf_name not in
('NULL','BF_control')"))
total = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

#venn.diagram(list("Strains" = total$orfs,"Proto-genes" = proto$orfs), 
#             fill = c("red", "green"),
#             alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,
#             filename = "trial.png")

overlap <- calculate.overlap(
  x = list(
    "Strains" = total$orfs,
    "Proto-genes" = proto$orfs
  )
);


##### END
dbDisconnect(conn)
