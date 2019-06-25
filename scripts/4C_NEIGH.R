##### NEIGHBOR EFFECT
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/16/2019

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggExtra)
library(gridExtra)
source("R/functions/initialize.sql.R")

##### GET/SET DATA
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/neigh/';
density = 6144;

##### CHECK POSITION WISE VARIABILITY
conn <- initialize.sql("saurin_test")

tablename_fit = sprintf('%s_%d_FITNESS',expt_name,density);
tablename_p2o = '4C3_pos2orf_name1';
tablename_bpos = '4C3_borderpos';

p2c_info = NULL
p2c_info[1] = '4C3_pos2coor6144'
p2c_info[2] = '6144plate'
p2c_info[3] = '6144col'
p2c_info[4] = '6144row'

hours = dbGetQuery(conn, sprintf('select distinct hours from %s order by hours asc', tablename_fit))
p2c = dbGetQuery(conn, sprintf('select * from %s a order by a.%s, a.%s, a.%s',
                               p2c_info[1],
                               p2c_info[2],
                               p2c_info[3],
                               p2c_info[4]))
n_plates = dbGetQuery(conn, sprintf('select distinct %s from %s a order by %s asc',
                                    p2c_info[2],
                                    p2c_info[1],
                                    p2c_info[2]))

for (hr in hours[[1]][9:length(hours[[1]])]) {
  for (pl in n_plates[[1]]) {
    alldat = dbGetQuery(conn, sprintf('select a.*, b.*
                                      from %s a, %s b
                                      where a.hours = %d
                                      and a.pos = b.pos
                                      and b.%s = %d
                                      order by b.%s, b.%s',
                                      tablename_fit,
                                      p2c_info[1],hr,p2c_info[2],
                                      pl,
                                      p2c_info[3],p2c_info[4]))
    
    alldat$source[alldat$`6144row`%%2==1 & alldat$`6144col`%%2==1] = 'TL'
    alldat$source[alldat$`6144row`%%2==0 & alldat$`6144col`%%2==1] = 'BL'
    alldat$source[alldat$`6144row`%%2==1 & alldat$`6144col`%%2==0] = 'TR'
    alldat$source[alldat$`6144row`%%2==0 & alldat$`6144col`%%2==0] = 'BR'
    
    alldat$colony[alldat$orf_name == 'BF_control'] = 'Reference'
    alldat$colony[alldat$orf_name != 'BF_control'] = 'Query'
    alldat$colony[is.na(alldat$orf_name)] = 'Gap'

    # cnt = 1
    for (sr in unique(alldat$source)) {
      temp <- alldat[alldat$source == sr,]
      for (i in seq(2,length(unique(temp$`6144col`)))) {
        col <- unique(temp$`6144col`)[i]
        lf <- tail(temp$`6144col`[temp$`6144col` < col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
        rt <- temp$`6144col`[temp$`6144col` > col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
        for (ii in seq(2,length(unique(temp$`6144row`[temp$`6144col` == col])))) {
          row <- unique(temp$`6144row`[temp$`6144col` == col])[ii]
          up <- tail(temp$`6144row`[temp$`6144row` < row & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
          dw <- temp$`6144row`[temp$`6144row` > row & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
          if (!is.na(alldat$average[alldat$`6144col` == col & alldat$`6144row` == row])) {
            a <- alldat$average[alldat$`6144col` == col & alldat$`6144row` == row]
            u <- alldat$average[alldat$`6144col` == col & alldat$`6144row` == up]
            d <- alldat$average[alldat$`6144col` == col & alldat$`6144row` == dw]
            l <- alldat$average[alldat$`6144col` == lf & alldat$`6144row` == row]
            r <- alldat$average[alldat$`6144col` == rt & alldat$`6144row` == row]
            alldat$var[alldat$`6144col` == col & alldat$`6144row` == row] <- sd(c(a,u,d,l,r),na.rm = T)/mean(c(a,u,d,l,r),na.rm = T)
            alldat$neigh[alldat$`6144col` == col & alldat$`6144row` == row] <-  mean(c(u,d,l,r),na.rm = T)
            alldat$diff[alldat$`6144col` == col & alldat$`6144row` == row] <- a - mean(c(u,d,l,r),na.rm = T)
            }
          # cnt <- cnt + 1
        }
      }
      diff_std <- sd(alldat$diff[alldat$source == sr],na.rm = T)
      diff_mean <- mean(alldat$diff[alldat$source == sr],na.rm = T)
      alldat$outlier[alldat$source == sr & !is.na(alldat$average) & alldat$diff > (diff_mean + 3*diff_std)] = 'Bigger'
      alldat$outlier[alldat$source == sr & !is.na(alldat$average) & alldat$diff < (diff_mean - 3*diff_std)] = 'Smaller'
      alldat$outlier[is.na(alldat$outlier)] = 'Normal'
      
      # ngh.sca <- ggplot(alldat[alldat$source == sr,]) +
      #   geom_abline(linetype = 2, col = 'red', lwd = 1.2) +
      #   geom_abline(linetype = 2, col = 'blue', lwd = 1,
      #               intercept = 3*diff_std) +
      #   geom_abline(linetype = 2, col = 'blue', lwd = 1,
      #               intercept = -3*diff_std) +
      #   geom_point(aes(x=average, y=neigh,
      #                  col = outlier,
      #                  shape = colony),
      #              size = 3, alpha = 0.7) +
      #   scale_color_manual(name="wrt Neigh Ref",
      #                      breaks=c("Smaller","Normal","Bigger"),
      #                      values=c("Smaller"="#FFA000","Bigger"="#009688","Normal"="grey50"),
      #                      guide = F) +
      #   scale_shape_manual(name = 'Colony Type',
      #                      breaks=c('Reference','Query','Gap'),
      #                      values=c('Reference' = 18,'Query' = 15, 'Gap' =1),
      #                      guide = F) +
      #   labs(title = "Comparison With Neighboring References :",
      #        x = "Colony Pixel Count",
      #        y = "Neighbor Pixel Count Average") +
      #   scale_x_continuous(breaks = seq(0,1000,100),
      #                      minor_breaks = seq(0,1000,25)) +
      #   scale_y_continuous(breaks = seq(0,1000,100),
      #                      minor_breaks = seq(0,1000,25)) +
      #   theme_linedraw() +
      #   theme(axis.text.x = element_text(size=10),
      #         axis.title.x = element_text(size=15),
      #         axis.text.y = element_text(size=10),
      #         axis.title.y = element_text(size=15),
      #         legend.position = c(0.8,0.2),
      #         legend.background = element_rect(fill="gray90",
      #                                          size=.5,
      #                                          linetype="dotted"),
      #         legend.text = element_text(size=10),
      #         legend.title =  element_text(size=15),
      #         plot.title = element_text(size=20,hjust = 0.5),
      #         plot.subtitle = element_text(size=13,hjust = 0.5)) +
      #   coord_cartesian(xlim = c(200,600),
      #                   ylim = c(200,600))
      # 
      # ngh.plt <- ggplot(alldat[alldat$source == sr,]) +
      #   geom_point(aes(x=`6144col`, y=`6144row`,
      #                  col = outlier,
      #                  shape = colony,
      #                  alpha = sr),
      #              size = 3) +
      #   geom_point(data = alldat,
      #              aes(x=`6144col`, y=`6144row`,
      #                  col = outlier,
      #                  shape = colony,
      #                  alpha = 'REST'),
      #              size = 3) +
      #   scale_color_manual(name="wrt Neigh Ref",
      #                      breaks=c("Smaller","Normal","Bigger"),
      #                      values=c("Smaller"="#FFA000","Bigger"="#009688","Normal"="grey50")) +
      #   scale_shape_manual(name = 'Colony Type',
      #                      breaks=c('Reference','Query','Gap'),
      #                      values=c('Reference' = 18,'Query' = 15, 'Gap' =1)) +
      #   scale_alpha_manual(name = 'Source',
      #                      breaks=c(sr,'REST'),
      #                      values=c(sr = 1,'REST'=0.4)) +
      #   scale_x_continuous(breaks = seq(0,96,2),limits = c(1,96)) +
      #   scale_y_continuous(breaks = seq(0,64,2),limits = c(64,1),trans = 'reverse') +
      #   labs(title = sprintf("%s | %d hours | Plate %d | %s",
      #                        expt, hr, pl, sr),
      #        x = "Column",
      #        y = "Row") +
      #   theme_linedraw() +
      #   theme(axis.text.x = element_text(size=10),
      #         axis.title.x = element_text(size=15),
      #         axis.text.y = element_text(size=10),
      #         axis.title.y = element_text(size=15),
      #         legend.position = 'right',
      #         legend.text = element_text(size=10),
      #         legend.title =  element_text(size=15),
      #         plot.title = element_text(size=20,hjust = 0.5),
      #         plot.subtitle = element_text(size=13,hjust = 0.5)) +
      #   guides(color = guide_legend(override.aes = list(size=3, alpha = 1)),
      #          shape = guide_legend(override.aes = list(size=3)),
      #          alpha = guide_legend(override.aes = list(size=3, alpha = c(1,0.3))))
      # 
      # fig<- ggarrange(ngh.sca, ngh.plt,
      #                 widths = c(1,1.5),
      #                 nrow = 1)
      # 
      # ggsave(sprintf("%s%s_NEIGH_OUTLIERS_%d_%d_%s.png",
      #                out_path,expt_name,hr,pl,sr),
      #        fig,
      #        width = 25,height = 10)
      # 
      # ggplot(alldat[alldat$source == sr,]) +
      #   geom_density(aes(x = average, col = source),lwd = 1.2) +
      #   geom_point(aes(x = average, y = 0, fill = outlier),
      #              col = 'transparent',
      #              alpha = 0.7,
      #              shape = 21,
      #              size = 3) +
      #   scale_fill_manual(name="wrt Neigh Ref",
      #                     breaks=c("Smaller","Normal","Bigger"),
      #                     values=c("Smaller"="#FFA000","Bigger"="#009688","Normal"="grey50")) +
      #   scale_color_manual(name="Source",
      #                      values=c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548"),
      #                      breaks=c("TL","TR","BL","BR"),
      #                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
      #   labs(title = "Outliers on Density Plot",
      #        subtitle = sprintf("%s | %d hours | Plate %d | %s",
      #                           expt, hr, pl, sr),
      #        x = 'Observed Pixel Count',
      #        y = 'Density') +
      #   theme_linedraw() +
      #   theme(axis.text.x = element_text(size=10),
      #         axis.title.x = element_text(size=15),
      #         axis.text.y = element_text(size=10),
      #         axis.title.y = element_text(size=15),
      #         legend.position = 'right',
      #         legend.text = element_text(size=10),
      #         legend.title =  element_text(size=15),
      #         plot.title = element_text(size=20,hjust = 0.5),
      #         plot.subtitle = element_text(size=13,hjust = 0.5)) +
      #   guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
      #   coord_cartesian(xlim = c(200,600),
      #                   ylim = c(0,0.016))
      # 
      # ggsave(sprintf("%s%s_NEIGH_OUT_DEN_%d_%d_%s.png",
      #                out_path,expt_name,hr,pl,sr),
      #        width = 10,height = 10)
    }
    
    for (o in alldat$pos[alldat$outlier == 'Bigger']) {
      c = alldat$`6144col`[alldat$pos == o]
      r = alldat$`6144row`[alldat$pos == o]
      num.gaps <- sum(alldat$colony[alldat$`6144row` == r - 1 & alldat$`6144col` == c |
                                      alldat$`6144row` == r + 1 & alldat$`6144col` == c |
                                      alldat$`6144row` == r & alldat$`6144col` == c - 1 |
                                      alldat$`6144row` == r & alldat$`6144col` == c + 1 |
                                      alldat$`6144row` == r - 1 & alldat$`6144col` == c - 1 |
                                      alldat$`6144row` == r + 1 & alldat$`6144col` == c + 1|
                                      alldat$`6144row` == r - 1 & alldat$`6144col` == c + 1 |
                                      alldat$`6144row` == r + 1 & alldat$`6144col` == c - 1] == 'Gap')
      if (num.gaps > 0) {
        alldat$gaps[alldat$pos == o] = num.gaps
      }
    }

    num.out <- sum(alldat$outlier != 'Normal')
    num.big <- sum(alldat$outlier == 'Bigger')
    big.gap <- sum(alldat$gaps > 0, na.rm = T)

    # ggplot(alldat[alldat$source == sr,]) +
    #   geom_point(data = alldat,
    #              aes(x=`6144col`, y=`6144row`,
    #                  col = outlier,
    #                  shape = colony,
    #                  alpha = outlier),
    #              size = 3) +
    #   scale_color_manual(name="wrt Neigh Ref",
    #                      breaks=c("Smaller","Normal","Bigger"),
    #                      values=c("Smaller"="#FFA000","Bigger"="#009688","Normal"="grey50")) +
    #   scale_shape_manual(name = 'Colony Type',
    #                      breaks=c('Reference','Query','Gap'),
    #                      values=c('Reference' = 18,'Query' = 15, 'Gap' =1)) +
    #   scale_alpha_manual(name="wrt Neigh Ref",
    #                      breaks=c("Smaller","Normal","Bigger"),
    #                      values=c("Smaller"=1,"Bigger"=1,"Normal"=0.7),
    #                      guide = F) +
    #   scale_x_continuous(breaks = seq(0,96,2),limits = c(1,96)) +
    #   scale_y_continuous(breaks = seq(0,64,2),limits = c(64,1),trans = 'reverse') +
    #   labs(title = "Comparison With Neighboring References",
    #        subtitle = sprintf("%s | %d hours | Plate %d | Outliers = %d | Bigger = %d | Bigger + Gap = %d",
    #                        expt, hr, pl, num.out, num.big, big.gap),
    #        x = "Column",
    #        y = "Row") +
    #   theme_linedraw() +
    #   theme(axis.text.x = element_text(size=10),
    #         axis.title.x = element_text(size=15),
    #         axis.text.y = element_text(size=10),
    #         axis.title.y = element_text(size=15),
    #         legend.position = 'right',
    #         legend.text = element_text(size=10),
    #         legend.title =  element_text(size=15),
    #         plot.title = element_text(size=20,hjust = 0.5),
    #         plot.subtitle = element_text(size=13,hjust = 0.5)) +
    #   guides(color = guide_legend(override.aes = list(size=3, alpha = 1)),
    #          shape = guide_legend(override.aes = list(size=3)))
    # ggsave(sprintf("%s%s_NEIGH_OUTLIERS_%d_%d.png",
    #                out_path,expt_name,hr,pl),
    #        width = 15,height = 10)
  }
}

for (o in alldat$pos) {
  c = alldat$`6144col`[alldat$pos == o]
  r = alldat$`6144row`[alldat$pos == o]
  num.gaps <- sum(alldat$colony[alldat$`6144row` == r - 1 & alldat$`6144col` == c |
                                  alldat$`6144row` == r + 1 & alldat$`6144col` == c |
                                  alldat$`6144row` == r & alldat$`6144col` == c - 1 |
                                  alldat$`6144row` == r & alldat$`6144col` == c + 1 |
                                  alldat$`6144row` == r - 1 & alldat$`6144col` == c - 1 |
                                  alldat$`6144row` == r + 1 & alldat$`6144col` == c + 1|
                                  alldat$`6144row` == r - 1 & alldat$`6144col` == c + 1 |
                                  alldat$`6144row` == r + 1 & alldat$`6144col` == c - 1] == 'Gap')
  if (num.gaps > 0) {
    alldat$gaps[alldat$pos == o] = num.gaps
  }
}

alldat$Gap[is.na(alldat$gaps)] = "NO"
alldat$Gap[!is.na(alldat$gaps)] = "YES"

# ggplot() +
#   geom_point(data = alldat[alldat$Gap == "NO",], aes(x = average, y = fitness, col = source, alpha = source)) +
#   geom_point(data = alldat[!is.na(alldat$gaps),], aes(x = average, y = fitness, col = "Gap", alpha = "Gap")) +
#   geom_smooth(data = alldat, aes(x = average, y = fitness), method = "lm", color = "black") +
#   scale_alpha_manual(values = c("Gap" = 1,"TL" = 0.3,"TR" = 0.3,"BL" = 0.3,"BR" = 0.3), guide = F) +
#   labs(x = "Colony Size",
#        y = "Fitness") +
#   scale_x_continuous(breaks = seq(0,800,100),
#                      minor_breaks = seq(0,800,50)) +
#   scale_y_continuous(breaks = seq(0,2,0.1),
#                      minor_breaks = seq(0,2,0.05)) +
#   scale_color_manual(values = c("TL"="#D32F2F","TR"="#536DFE","BL"="#388E3C","BR"="#795548",
#                                 "Gap" = "black")) +
#   theme_linedraw() +
#   coord_cartesian(xlim = c(200,700),
#                   ylim = c(0.8,1.2))

  
avg.den.ori <- ggplot() +
  geom_point(data = alldat, aes(x = average, col = Gap), stat = "density")


median(alldat$average[alldat$Gap ==  "YES"], na.rm = T)
median(alldat$average[alldat$Gap ==  "NO"], na.rm = T)

alldat$average_mcg <- alldat$average

alldat$average_mcg[alldat$Gap ==  "YES"] = alldat$average_mcg[alldat$Gap ==  "YES"] *
  median(alldat$average_mcg[alldat$Gap ==  "NO"], na.rm = T)/median(alldat$average_mcg[alldat$Gap ==  "YES"], na.rm = T)

avg.den.mcg <- ggplot() +
  geom_point(data = alldat, aes(x = average_mcg, col = Gap), stat = "density")

avg.ul <- mean(alldat$average, na.rm = T) + 2*sd(alldat$average, na.rm = T)
avg.ll <- mean(alldat$average, na.rm = T) - 2*sd(alldat$average, na.rm = T)

alldat$extreme[alldat$average > avg.ul] <- 'Big'

avg.mcg.ul <- mean(alldat$average_mcg, na.rm = T) + 2*sd(alldat$average_mcg, na.rm = T)
avg.mcg.ll <- mean(alldat$average_mcg, na.rm = T) - 2*sd(alldat$average_mcg, na.rm = T)

alldat$extreme_mcg[alldat$average_mcg > avg.mcg.ul] <- 'Big'

fit.ul <- mean(alldat$fitness, na.rm = T) + 2*sd(alldat$fitness, na.rm = T)
fit.ll <- mean(alldat$fitness, na.rm = T) - 2*sd(alldat$fitness, na.rm = T)

alldat$extreme_fit[alldat$fitness > fit.ul] <- 'Big'
alldat$extreme_fit[alldat$fitness < fit.ll] <- 'Small'

ggplot(alldat) +
  geom_point(data = alldat,
             aes(x=`6144col`, y=`6144row`,
                 col = extreme_fit,
                 shape = colony),
             size = 3)
  

##### CHANGING DATA ON SQL USING GAP MEDIAN CORRECTION
fitdat = dbGetQuery(conn, sprintf('select a.*, b.*
                                      from %s a, %s b
                                      where a.pos = b.pos
                                      order by b.%s, b.%s',
                                  tablename_fit,
                                  p2c_info[1],
                                  p2c_info[3],p2c_info[4]))

fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==1] = 'TL'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==1] = 'BL'
fitdat$source[fitdat$`6144row`%%2==1 & fitdat$`6144col`%%2==0] = 'TR'
fitdat$source[fitdat$`6144row`%%2==0 & fitdat$`6144col`%%2==0] = 'BR'

fitdat$colony[fitdat$orf_name == 'BF_control'] = 'Reference'
fitdat$colony[fitdat$orf_name != 'BF_control'] = 'Query'
fitdat$colony[is.na(fitdat$orf_name)] = 'Gap'

for (hr in unique(fitdat$hours)) {
  for (pl in unique(fitdat$`6144plate`[fitdat$hours == hr])) {
    temp <- fitdat[fitdat$hours == hr & fitdat$`6144plate` == pl,]
    for (o in temp$pos) {
      c = temp$`6144col`[temp$pos == o]
      r = temp$`6144row`[temp$pos == o]
      num.gaps <- sum(temp$colony[temp$`6144row` == r - 1 & temp$`6144col` == c |
                                      temp$`6144row` == r + 1 & temp$`6144col` == c |
                                      temp$`6144row` == r & temp$`6144col` == c - 1 |
                                      temp$`6144row` == r & temp$`6144col` == c + 1 |
                                      temp$`6144row` == r - 1 & temp$`6144col` == c - 1 |
                                      temp$`6144row` == r + 1 & temp$`6144col` == c + 1|
                                      temp$`6144row` == r - 1 & temp$`6144col` == c + 1 |
                                      temp$`6144row` == r + 1 & temp$`6144col` == c - 1] == 'Gap')
      if (num.gaps > 0) {
        temp$gaps[temp$pos == o] = num.gaps
        # fitdat$gaps[fitdat$hours == hr & fitdat$`6144plate` == pl & fitdat$pos == o] = num.gaps
      }
    }
    temp$NearGap[is.na(temp$gaps)] = "NO"
    temp$NearGap[!is.na(temp$gaps)] = "YES"
    temp$average[temp$NearGap ==  "YES"] = temp$average[temp$NearGap ==  "YES"] *
      median(temp$average[temp$NearGap ==  "NO"], na.rm = T)/median(temp$average[temp$NearGap ==  "YES"], na.rm = T)
    fitdat$mca[fitdat$hours == hr & fitdat$`6144plate` == pl] = temp$average
  }
}

ggplot(temp) +
  geom_point(aes(x = mca, col = NearGap), stat = "density")

# dbWriteTable(conn, "4C3_GA1_MCG_6144_FITNESS", fitdat[1:6], overwrite = T)

##### MEDIAN CORRECTED FITNESS V/S ORIGINAL
mcgdat = dbGetQuery(conn, sprintf('select a.*, b.*
                                      from %s a, %s b
                                      where a.pos = b.pos
                                      order by b.%s, b.%s',
                                  tablename_fit,
                                  p2c_info[1],
                                  p2c_info[3],p2c_info[4]))
