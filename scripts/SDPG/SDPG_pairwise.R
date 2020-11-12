library(effsize)
load("/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/results.RDATA")
# lidres <- mcatres
##### CLEAN DATA
lidres <- lidres[lidres$orf_name != 'YHR021W-A' &
                   lidres$orf_name != 'BOR' &
                   !is.na(lidres$orf_name),]
# Removing plates with bad colony grids
lidres <- lidres[!(lidres$density == 6144 & lidres$arm == 'GLU' & lidres$rep != 'R3' & lidres$hours >= 20),]
lidres <- lidres[!(lidres$density == 6144 & lidres$arm == 'CAS' & lidres$rep == 'R2' & lidres$hours >= 20),]

# Saturation information
lidres$time <- NA
lidres$time[lidres$density == 1536 & lidres$hours == 48 & lidres$arm != 'SDA'] <- "Saturated"
lidres$time[lidres$density == 1536 & lidres$hours %in% c(72, 78) & lidres$arm == 'SDA'] <- "Saturated"
lidres$time[lidres$density == 6144 & lidres$hours == 24 & lidres$arm != 'SDA'] <- "Saturated"
lidres$time[lidres$density == 6144 & lidres$hours == 36 & lidres$arm == 'SDA'] <- "Saturated"
lidres$time[is.na(lidres$time)] <- "Other"

lidres <- lidres[lidres$time == "Saturated",]

# Removing outliers
for (d in unique(lidres$density)) {
  for (a in unique(lidres$arm[lidres$density == d])) {
    for (r in unique(lidres$rep[lidres$density == d & 
                                lidres$arm == a])) {
      for (h in unique(lidres$hours[lidres$density == d & 
                                    lidres$arm == a & 
                                    lidres$rep == r])) {
        for (o in unique(lidres$orf_name[lidres$density == d & 
                                         lidres$arm == a & 
                                         lidres$rep == r &
                                         lidres$hours == h])) {
          temp <- lidres$fitness[lidres$density == d & 
                                   lidres$arm == a & 
                                   lidres$rep == r &
                                   lidres$hours == h &
                                   lidres$orf_name == o]
          m <- median(temp, na.rm = T)
          madev <- mad(temp)
          ul <- m + 3*madev
          ll <- m - 3*madev
          
          lidres$fitness[lidres$density == d & 
                           lidres$arm == a & 
                           lidres$rep == r &
                           lidres$hours == h &
                           lidres$orf_name == o &
                           lidres$fitness > ul] <- NA
          lidres$fitness[lidres$density == d & 
                           lidres$arm == a & 
                           lidres$rep == r &
                           lidres$hours == h &
                           lidres$orf_name == o &
                           lidres$fitness < ll] <- NA
          
          
        }
        cont_median <- median(lidres$average[lidres$density == d & 
                                               lidres$arm == a & 
                                               lidres$rep == r &
                                               lidres$hours == h &
                                               lidres$orf_name == 'BF_control'], na.rm = T)
        lidres$ccs[lidres$density == d & 
                     lidres$arm == a & 
                     lidres$rep == r &
                     lidres$hours == h] <- lidres$cs_median[lidres$density == d & 
                                                              lidres$arm == a & 
                                                              lidres$rep == r &
                                                              lidres$hours == h] * cont_median
      }
    }
  }
}

# Factorize the Experimental Arms
lidres$arm <- factor(lidres$arm, levels = c('GLU','CAS','SDA'))

comparisons <- NULL
for (d in unique(lidres$density)) {
  for (a in unique(lidres$arm[lidres$density == d])) {
    for (r in unique(lidres$rep[lidres$density == d & lidres$arm == a])) {
      for (o in unique(sort(lidres$orf_name[lidres$density == d & lidres$arm == a & lidres$rep == r])))
        comparisons <- rbind(comparisons, c(d,a,r,o))
    }
  }
}

pair_comp <- NULL
for (i in 1:dim(comparisons)[1]) {
  g1_fit <- lidres$fitness[lidres$density == comparisons[i,1] & 
                              lidres$arm == comparisons[i,2] & 
                              lidres$rep == comparisons[i,3] & 
                              lidres$orf_name == comparisons[i,4]]
  for (ii in 1:dim(comparisons)[1]) {
    g2_fit <- lidres$fitness[lidres$density == comparisons[ii,1] & 
                                lidres$arm == comparisons[ii,2] & 
                                lidres$rep == comparisons[ii,3] & 
                                lidres$orf_name == comparisons[ii,4] &
                                !is.na(lidres$fitness)]
    sig <- wilcox.test(g1_fit, g2_fit, conf.level = 0.95)
    es <- cliff.delta(g2_fit, g1_fit, conf.level = 0.95,
                      use.unbiased=TRUE, use.normal=FALSE, return.dm=FALSE)
    
    # pair_comp <- rbind(pair_comp, c(comparisons[i,],
    #                                 comparisons[ii,],
    #                                 sig$p.value,
    #                                 es$estimate,
    #                                 as.character(es$magnitude)))
    
    pair_comp <- rbind(pair_comp, cbind(g1_density = comparisons[i,1],
                                        g1_arm = comparisons[i,2], 
                                        g1_rep = comparisons[i,3], 
                                        g1_orf_name = comparisons[i,4],
                                        g2_density = comparisons[ii,1],
                                        g2_arm = comparisons[ii,2], 
                                        g2_rep = comparisons[ii,3], 
                                        g2_orf_name = comparisons[ii,4],
                                        pvalue = sig$p.value,
                                        cliffdelta = es$estimate,
                                        magnitude = as.character(es$magnitude)))
  }
}

save(pair_comp, file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/pair_comp.RDATA')
