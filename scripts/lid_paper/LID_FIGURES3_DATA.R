##### LID PAPER FIGURES #3 DATA
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 02/14/2020

##### INITIALIZE
library(RMariaDB)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/paper/'

##### GETTING ALL THE SOURCENORMALIZATION DATA
dat.avg <- dbGetQuery(conn, 'select b.*, orf_name, hours, fitness,
                      bg, average
                      from
                      4C4_FS_NOSN_6144_FITNESS a, 4C4_pos2coor b
                      where a.pos = b.pos
                      order by plate, col, row')
dat.avg$name <- '1. Raw Colony Size'
dat.avg$method <- 1

dat.nosn <- dbGetQuery(conn, 'select b.*, orf_name, hours, fitness,
                       bg, average
                       from
                       4C4_FS_NOSN_6144_FITNESS a, 4C4_pos2coor b
                       where a.pos = b.pos
                       and fitness between 0 and 2
                       order by plate, col, row')
dat.nosn$name <- '2. LID-SN'
dat.nosn$method <- 2

dat.nocc <- dbGetQuery(conn, 'select b.*, orf_name, hours, fitness,
                       bg, average
                       from
                       4C4_FS_6144_FITNESS a, 4C4_pos2coor b
                       where a.pos = b.pos
                       and fitness between 0 and 2
                       order by plate, col, row')
dat.nocc$name <- '3. LID-AC'
dat.nocc$method <- 3

dat.lid <- dbGetQuery(conn, 'select b.*, orf_name, hours, fitness,
                      bg, average
                      from
                      4C4_FS_CC_6144_FITNESS a, 4C4_pos2coor b
                      where a.pos = b.pos
                      and fitness between 0 and 2
                      order by plate, col, row')
dat.lid$name <- '4. LID'
dat.lid$method <- 4

dat.bean <- dbGetQuery(conn, 'select b.*, orf_name, hours, fitness,
                       bg, average
                       from
                       4C4_FS_BEAN_6144_FITNESS a, 4C4_pos2coor b
                       where a.pos = b.pos
                       and fitness between 0 and 2
                       order by plate, col, row')
dat.bean$name <- '5. MCAT'
dat.bean$method <- 5

dat.sn <- rbind(dat.avg[dat.avg$hours == 11.04,],
                dat.nosn[dat.nosn$hours == 11.04,],
                dat.nocc[dat.nocc$hours == 11.04,],
                dat.lid[dat.lid$hours == 11.04,],
                dat.bean[dat.bean$hours == 11.04,])
dat.sn$source[dat.sn$row%%2==1 & dat.sn$col%%2==1] = '1TL'
dat.sn$source[dat.sn$row%%2==0 & dat.sn$col%%2==1] = '3BL'
dat.sn$source[dat.sn$row%%2==1 & dat.sn$col%%2==0] = '2TR'
dat.sn$source[dat.sn$row%%2==0 & dat.sn$col%%2==0] = '4BR'

save(dat.sn,
     file = sprintf('%sSOURCENORMALIZATIONDATA.RData',out_path))

##### BACKGROUND PREDICTION DATA
dat.nocc <- dbGetQuery(conn, 'select b.*, orf_name, hours, fitness,
                      bg, average
                       from
                       4C4_FS_6144_FITNESS a, 4C4_pos2coor b
                       where a.pos = b.pos
                       and fitness between 0 and 2
                       order by plate, col, row')
dat.nocc$name <- '5. LID-AC'
dat.nocc$method <- 5

dat.rnd <- dbGetQuery(conn, 'select b.*, orf_name, hours, fitness,
                      bg, average
                      from
                      4C4_FS_6144_FITNESS a, 4C4_pos2coor b
                      where a.pos = b.pos
                      order by plate, col, row')
dat.rnd$name <- '6. RND'
dat.rnd$method <- 6

for (hr in unique(dat.rnd$hours)) {
  n <- length(dat.rnd$average[dat.rnd$hours == hr])
  m <- mean(dat.rnd$average[dat.rnd$hours == hr], na.rm = T)
  s <- sd(dat.rnd$average[dat.rnd$hours == hr], na.rm = T)
  dat.rnd$bg[dat.rnd$hours == hr] <- rnorm(n,m,s)
  dat.rnd$fitness[dat.rnd$hours == hr] <- dat.rnd$average[dat.rnd$hours == hr]/
    dat.rnd$bg[dat.rnd$hours == hr]
}

dat.bg <- rbind(dat.rnd,dat.nosn,dat.nocc,dat.lid,dat.bean)
save(dat.bg,
     file = sprintf('%sBACKGROUNDDATA.RData',out_path))

rmse <- NULL
i = 1
for (m in unique(dat.bg$method)) {
  for (hr in unique(dat.bg$hours[dat.bg$method == m])) {
    rmse$method[i] <- m
    rmse$name[i] <- unique(dat.bg$name[dat.bg$method == m])
    rmse$hour[i] <- hr
    rmse$avg[i] <- mean(dat.bg$average[dat.bg$hours == hr & dat.bg$method == m & dat.bg$average > 0], na.rm = T)
    rmse$rmse[i] <- sqrt(mean((dat.bg$average[dat.bg$hours == hr & dat.bg$method == m & dat.bg$average > 0] -
                                 dat.bg$bg[dat.bg$hours == hr & dat.bg$method == m & dat.bg$average > 0])^2,na.rm = T))
    i = i + 1
  }
}
rmse <- data.frame(rmse)
rmse$per <- rmse$rmse/rmse$avg * 100

save(rmse,
     file = sprintf('%sRMSEDATA.RData',out_path))


##### CHANGE IN SENSITIVITY WITH METHOD STAGES
dat.method <- NULL
i = 1
dat.method$method[i] <- 'NO-NORM'
dat.method$es[i] <- .95
dat.method$sen[i] <- 10
dat.method$method[i+1] <- 'NO-NORM'
dat.method$es[i+1] <- 1.05
dat.method$sen[i+1] <- 7
i = i + 2
dat.method$method[i] <- 'LID-SN'
dat.method$es[i] <- .95
dat.method$sen[i] <- 55
dat.method$method[i + 1] <- 'LID-SN'
dat.method$es[i + 1] <- 1.05
dat.method$sen[i + 1] <- 26
i = i + 2
dat.method$method[i] <- 'LID-AC'
dat.method$es[i] <- .95
dat.method$sen[i] <- 70
dat.method$method[i+1] <- 'LID-AC'
dat.method$es[i+1] <- 1.05
dat.method$sen[i+1] <- 66
i = i + 2
dat.method$method[i] <- 'LID'
dat.method$es[i] <- .95
dat.method$sen[i] <- 70
dat.method$method[i+1] <- 'LID'
dat.method$es[i+1] <- 1.05
dat.method$sen[i+1] <- 70

dat.method$rep <- 8

i = i + 2
dat.method$method[i] <- 'NO-NORM'
dat.method$es[i] <- .95
dat.method$sen[i] <- 22
dat.method$method[i+1] <- 'NO-NORM'
dat.method$es[i+1] <- 1.05
dat.method$sen[i+1] <- 8
i = i + 2
dat.method$method[i] <- 'LID-SN'
dat.method$es[i] <- .95
dat.method$sen[i] <- 82
dat.method$method[i + 1] <- 'LID-SN'
dat.method$es[i + 1] <- 1.05
dat.method$sen[i + 1] <- 79
i = i + 2
dat.method$method[i] <- 'LID-AC'
dat.method$es[i] <- .95
dat.method$sen[i] <- 96
dat.method$method[i+1] <- 'LID-AC'
dat.method$es[i+1] <- 1.05
dat.method$sen[i+1] <- 91
i = i + 2
dat.method$method[i] <- 'LID'
dat.method$es[i] <- .95
dat.method$sen[i] <- 98
dat.method$method[i+1] <- 'LID'
dat.method$es[i+1] <- 1.05
dat.method$sen[i+1] <- 96

dat.method <- data.frame(dat.method)
dat.method$rep[9:16] <- 16

save(dat.method,
     file = sprintf('%sMETHODSENDATA.RData',out_path))


##### FIGURE 6 : VIRTUAL PLATE 2 RESULTS
load("output/4C4_rnd_fit_data.RData")
load("output/4C4_rnd_ref_data.RData")
load("output/4C4_rnd_sta_data.RData")

sta.data$res <- sta.data$effect_p
sta.data$res[sta.data$effect_p == 'Beneficial' & sta.data$from == 'more'] <- '5TrueBeneficial'
sta.data$res[sta.data$effect_p == 'Beneficial' & sta.data$from == 'less'] <- '3FalseBeneficial'
sta.data$res[sta.data$effect_p == 'Deleterious' & sta.data$from == 'more'] <- '4FalseDeleterious'
sta.data$res[sta.data$effect_p == 'Deleterious' & sta.data$from == 'less'] <- '2TrueDeleterious'
sta.data$res[sta.data$effect_p == 'Neutral' & sta.data$from == 'more'] <- '6BenNeutral'
sta.data$res[sta.data$effect_p == 'Neutral' & sta.data$from == 'less'] <- '1DelNeutral'

sta.data <- sta.data[!sta.data$hours == 0,]
sta.data$method[sta.data$method == "BEAN"] <- "MCAT"

pie.dat <- data.frame(rbind(cbind(hours = sta.data$hours[sta.data$method == 'MCAT'],
                                  value = sta.data$from[sta.data$method == 'MCAT'],
                                  method = 'TRUTH'),
                            cbind(hours = sta.data$hours[sta.data$method == 'MCAT'],
                                  value = sta.data$res[sta.data$method == 'MCAT'],
                                  method = 'MCAT'),
                            cbind(hours = sta.data$hours[sta.data$method == 'LID'],
                                  value = sta.data$res[sta.data$method == 'LID'],
                                  method = 'LID')))

pie.dat$value[pie.dat$value == 'less'] <- '2TrueDeleterious'
pie.dat$value[pie.dat$value == 'more'] <- '5TrueBeneficial'

pie.dat <- arrange(transform(pie.dat,
                             method=factor(method,levels=c('LID','TRUTH','MCAT'))),
                   method)

pie.per <- plyr::count(pie.dat, vars = c('method','value'))
pie.per$per <- pie.per$freq/(sum(pie.per$freq)/3) * 100

save(pie.per,
     file = sprintf('%sVP2PIEDATA.RData',out_path))


rnd_data <- dbGetQuery(conn, 'select * from 4C4_FS_RND_6144_FITNESS a, 4C4_pos2coor b
                       where a.pos = b.pos')
i <- 1
for (h in unique(rnd_data$hours)) {
  rnd_data$plate[rnd_data$hours == h] <- i
  i <- i + 1
}
rnd_data$colony[rnd_data$orf_name == 'BF_control'] <- 'Reference'
rnd_data$colony[rnd_data$orf_name != 'BF_control'] <- 'Query'

save(rnd_data,
     file = sprintf('%sVP2RNDDATA.RData',out_path))

##### FIGURE S2: VIRTUAL PLATE 1 EXAMPLE
load("/home/sbp29/R/Projects/adaptivefitness/output/4C4_FS_NONORM_FITNESS.RData")
fit.all$colony[fit.all$orf_name == 'BF_control'] <- 'Reference'
fit.all$colony[fit.all$orf_name != 'BF_control'] <- 'Query'

vir.sml <- fit.all[fit.all$cont_hrs == 2.9 & fit.all$hours == 2.9 & fit.all$plate == 1,]
vir.sml$kind <- '1. Time = 2.9 hr'
vir.big <- fit.all[fit.all$cont_hrs == 11.04 & fit.all$hours == 11.04 & fit.all$plate == 1,]
vir.big$kind <- '2. Time = 11.04 hr'
vir.mix <- fit.all[fit.all$cont_hrs == 11.04 & fit.all$hours == 2.9 & fit.all$plate == 1,]
vir.mix$kind <- '3. Virtual Plate'

vir1.eg <- rbind(vir.sml,vir.big,vir.mix)

save(vir1.eg,
     file = sprintf('%sVP1EGDATA.RData',out_path))


##### FIGURE SX: Prescreens Plate Effect
stage <- 'PS1'
data.ps1 <- dbGetQuery(conn, sprintf('select * from 4C4_%s_1536_JPEG a, 4C4_pos2coor b
                                     where a.pos = b.pos
                                     order by hours, b.plate, b.col, b.row', stage))
data.ps1$stage <- "Upscale Plates (#1)"

stage <- 'PS2'
data.ps2 <- dbGetQuery(conn, sprintf('select * from 4C4_%s_1536_JPEG a, 4C4_pos2coor b
                                     where a.pos = b.pos
                                     order by hours, b.plate, b.col, b.row', stage))
data.ps2$stage <- "Transition Plates (#2)"

data.ps <- rbind(data.ps1, data.ps2[data.ps2$hours == 21,])
data.ps$source[data.ps$row%%2==1 & data.ps$col%%2==1] = '1TL'
data.ps$source[data.ps$row%%2==0 & data.ps$col%%2==1] = '3BL'
data.ps$source[data.ps$row%%2==1 & data.ps$col%%2==0] = '2TR'
data.ps$source[data.ps$row%%2==0 & data.ps$col%%2==0] = '4BR'

data.ps$stage <- factor(data.ps$stage, levels = c("Upscale Plates (#1)", "Transition Plates (#2)"))

for (s in unique(data.ps$stage)) {
  data.ps$mean_cs[data.ps$stage == s] <- mean(data.ps$average[data.ps$stage == s], na.rm = T)
  data.ps$median_cs[data.ps$stage == s] <- median(data.ps$average[data.ps$stage == s], na.rm = T)
}

save(data.ps,
     file = sprintf('%s4C4PSDATA.RData',out_path))


##### PLATE SURFACE
dat.raw <- dbGetQuery(conn, 'select * from 4C4_FS_6144_RAW a, 4C4_pos2coor b
                     where hours = 11.04 and a.pos = b.pos
                      order by hours, plate, col, row')

dat.raw$source[dat.raw$row%%2==1 & dat.raw$col%%2==1] = 'Top Left'
dat.raw$source[dat.raw$row%%2==0 & dat.raw$col%%2==1] = 'Bottom Left'
dat.raw$source[dat.raw$row%%2==1 & dat.raw$col%%2==0] = 'Top Right'
dat.raw$source[dat.raw$row%%2==0 & dat.raw$col%%2==0] = 'Bottom Right'

save(dat.raw,
     file = sprintf('%s4C4SURFACE.RData',out_path))


###### RND ES DATA
rnd.es <- dbGetQuery(conn, 'select c.*, d.rnd_hrs from
                    (select a.orf_name, a.hours, a.p lid_p, a.stat lid_stat, a.es lid_es, a.pix_es lid_pix_es,
                     b.p mcat_p, b.stat mcat_stat, b.es mcat_es, b.pix_es mcat_pix_es
                     from 4C4_FS_RND_6144_PVALUE a, 4C4_FS_RND_BEAN_6144_PVALUE b
                     where a.orf_name = b.orf_name and a.hours = b.hours) c
                     left join
                     (select orf_name, max(hours) hours, max(rnd_hrs) rnd_hrs
                     from 4C4_FS_RND_6144_DATA
                     where orf_name is not NULL
                     group by hours, orf_name) d
                     on c.orf_name = d.orf_name and c.hours = d.hours')

rnd.es$lid_effect[rnd.es$lid_p <= 0.05 & rnd.es$lid_stat > 0] <- 'Beneficial'
rnd.es$lid_effect[rnd.es$lid_p <= 0.05 & rnd.es$lid_stat < 0] <- 'Deleterious'
rnd.es$lid_effect[is.na(rnd.es$lid_effect)] <- 'Neutral'

rnd.es$mcat_effect[rnd.es$mcat_p <= 0.05 & rnd.es$mcat_stat > 0] <- 'Beneficial'
rnd.es$mcat_effect[rnd.es$mcat_p <= 0.05 & rnd.es$mcat_stat < 0] <- 'Deleterious'
rnd.es$mcat_effect[is.na(rnd.es$mcat_effect)] <- 'Neutral'

rnd.es$truth[rnd.es$hours < rnd.es$rnd_hrs] <- 'Beneficial'
rnd.es$truth[rnd.es$hours > rnd.es$rnd_hrs] <- 'Deleterious'

rnd.es$lid_result <- paste(strtrim(rnd.es$lid_effect,3), strtrim(rnd.es$truth,3), sep = '/')
rnd.es$mcat_result <- paste(strtrim(rnd.es$mcat_effect,3), strtrim(rnd.es$truth,3), sep = '/')

save(rnd.es,
     file = sprintf('%sVP2PIXES.RData',out_path))


