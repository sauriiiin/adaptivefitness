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