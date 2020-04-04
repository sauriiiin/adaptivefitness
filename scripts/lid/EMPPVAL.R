##### EMPIRICAL P VALUE
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 03/23/2020

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
source("R/functions/initialize.sql.R")
source("R/functions/isoutlier.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/paper/'

expt_name <- '4C4'
stage <- 'FS_CC'
density <- 6144
cont.name <- 'BF_control'
tablename_p2o <- sprintf('%s_pos2orf_name', expt_name)
tablename_bpos <- sprintf('%s_borderpos', expt_name)

tablename_fit <- sprintf('%s_%s_%d_FITNESS',
                         expt_name, stage, density)

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### 
contpos = dbGetQuery(conn, sprintf('select pos from %s where orf_name = "%s" and pos < 100000
                                   and pos not in (select pos from %s)',
                                   tablename_p2o,cont.name,tablename_bpos))

temp1536 <- NULL
contreps <- NULL
i <- 1
if (density == 1536) {
  repeat{contreps <- cbind(contreps, contpos$pos + 10000*i); i <- i + 1; if (i == 9) {break}}
} else {
  repeat{temp1536 <- cbind(temp1536, contpos$pos + 10000*i); i <- i + 1; if (i == 9) {break}}
  i <- 1
  repeat{contreps <- cbind(contreps, temp1536 + 100000*i); i <- i + 1; if (i == 9) {break}}
}
contreps <- data.frame(contreps)

hours <- 11.04
contfit <- NULL
for (ii in 1:dim(contreps)[1]) {
  cp <- paste(contreps[ii,],',',collapse = '')
  temp <- dbGetQuery(conn, sprintf('select fitness from %s
                             where hours = %0.2f and pos in (%s)
                             and fitness is not null',
                             tablename_fit,hours,
                             strtrim(cp,nchar(cp)-1)))
  if (length(temp$fitness) > 0) {
    outlier <- isoutlier(temp$fitness)
    temp$fitness[outlier] <- NA
    contfit <- c(contfit, mean(temp$fitness, na.rm = T))
  }
}
contmean <- mean(contfit, na.rm = T)
contstd <- sd(contfit, na.rm = T)

contfit <- data.frame(fitness = contfit)
save(contfit,
     file = sprintf('%sCONTFIT.RData',out_path))

