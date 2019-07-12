
ggplot(alldat) + 
  geom_histogram(aes(x = diff)) +
  geom_line(aes(x = diff), stat = 'density') +
  geom_vline(xintercept = median(alldat$diff, na.rm = T))

ggplot(alldat) + 
  geom_histogram(aes(x = coef)) +
  geom_line(aes(x = coef), stat = 'density') +
  geom_vline(xintercept = c(quantile(alldat$coef, .50, na.rm = T),
                            quantile(alldat$coef, .05, na.rm = T),
                            quantile(alldat$coef, .95, na.rm = T)))

ggplot(alldat) +
  # geom_histogram(aes(x = average)) +
  geom_line(aes(x = average), stat = 'density') +
  geom_vline(xintercept = mean(alldat$average, na.rm = T)) +
  geom_line(aes(x = mca), stat = 'density', col = 'red') +
  geom_vline(xintercept = mean(alldat$mca, na.rm = T), col = 'red')

dim(alldat[alldat$coef[!is.na(alldat$coef)] > quantile(alldat$coef, .95, na.rm = T),])
    
hello <- alldat[alldat$nearSick != 'N' & alldat$nearBig == 'B',]
hello2 <- alldat[alldat$nearBig != 'N' & alldat$nearSick == 'S',]
hello3 <- alldat[alldat$nearBig != 'N' & alldat$nearBig != 'B2' &
                   alldat$nearSick != 'N' & alldat$nearSick != 'S2',]


#### AFTER MCA
alldat2 <- alldat
temp <- alldat2
# cnt <- 0
for (i in seq(1,dim(temp)[1])) {
  col <- temp$`6144col`[i]
  row <- temp$`6144row`[i]
  
  lf <- tail(temp$`6144col`[temp$`6144row` == row & temp$`6144col` < col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
  rt <- temp$`6144col`[temp$`6144row` == row & temp$`6144col` > col & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
  if (length(lf) == 0) {
    up = NaN
  }
  if (length(rt) == 0) {
    up = NaN
  }
  
  up <- tail(temp$`6144row`[temp$`6144col` == col & temp$`6144row` < row & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)],1)
  dw <- temp$`6144row`[temp$`6144col` == col & temp$`6144row` > row & temp$orf_name == 'BF_control' & !is.na(temp$orf_name)][1]
  if (length(up) == 0) {
    up = NaN
  }
  if (length(dw) == 0) {
    up = NaN
  }
  
  a <- alldat2$mca[alldat2$`6144col` == col & alldat2$`6144row` == row]
  u <- mean(alldat2$mca[alldat2$`6144col` == col & alldat2$`6144row` == up], na.rm = T)
  d <- mean(alldat2$mca[alldat2$`6144col` == col & alldat2$`6144row` == dw], na.rm = T)
  l <- mean(alldat2$mca[alldat2$`6144col` == lf & alldat2$`6144row` == row], na.rm = T)
  r <- mean(alldat2$mca[alldat2$`6144col` == rt & alldat2$`6144row` == row], na.rm = T)
  alldat2$var[alldat2$`6144col` == col & alldat2$`6144row` == row] <- sd(c(a,u,d,l,r),na.rm = T)/median(c(a,u,d,l,r),na.rm = T)
  alldat2$neigh[alldat2$`6144col` == col & alldat2$`6144row` == row] <-  median(c(u,d,l,r),na.rm = T)
  # alldat2$wt_neigh[alldat2$`6144col` == col & alldat2$`6144row` == row] <-
  #   weighted.mean(c(u,d,l,r), c(1/abs(row-up),1/abs(row-dw),1/abs(col-lf),1/abs(col-rt)),na.rm = T)
  alldat2$diff[alldat2$`6144col` == col & alldat2$`6144row` == row] <- a - median(c(u,d,l,r),na.rm = T)
  
}

diff_std <- sd(alldat2$diff,na.rm = T)
diff_mean <- mean(alldat2$diff,na.rm = T)
alldat2$coef <- alldat2$mca/alldat2$neigh

