##### CONTROL DISTRIBUTION (4C)
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 04/26/2019

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("ggplot2")
library(RMariaDB)
library(ggplot2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### GET/SET DATA
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/';
density = 6144;

tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_p2o = '4C3_pos2orf_name1';
tablename_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

cont.name  = 'BF_control'
dbFetch(query, n=-1)

query = dbSendQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_fit))
hours = dbFetch(query,n=-1)

query = dbSendQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
                                  p2c_info[1],
                                  p2c_info[2],
                                  p2c_info[3],
                                  p2c_info[4]))
p2c = dbFetch(query,n=-1)

query = dbSendQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                  p2c_info[2],
                                  p2c_info[1],
                                  p2c_info[2]))
n_plates = dbFetch(query, n=-1) 
