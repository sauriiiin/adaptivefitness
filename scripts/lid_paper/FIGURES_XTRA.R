rnd.es[(rnd.es$lid_pix_es > 0 & rnd.es$truth == 'Deleterious'),]


fit.den.dat <- dbGetQuery(conn, 'select * from 4C4_FS_CC_6144_FITNESS
                          where hours > 0')

fit.den.dat$colony[fit.den.dat$orf_name == 'BF_control'] <- 'Reference'
fit.den.dat$colony[fit.den.dat$orf_name != 'BF_control'] <- 'Mutant'

ggplot(fit.den.dat,
       aes(x = fitness)) +
  geom_line(aes(col = colony), stat = 'density', trim = T) +
  facet_wrap(.~hours, scales = 'free')

lid.rnd.sen.data <- NULL
mcat.rnd.sen.data <- NULL
i <- 1
for (ref.h in unique(rnd.es$hours)) {
  for (mut.h in unique(rnd.es$rnd_hrs[rnd.es$hours == ref.h])) {
    rnd.es$es[rnd.es$hours == ref.h & rnd.es$rnd_hrs == mut.h] <-
      virP1$es[virP1$ref_hrs == ref.h & virP1$que_hrs == mut.h]
    # if (ref.h > 2 & mut.h > 2) {
      # rnd.es$cen[rnd.es$hours == ref.h & rnd.es$rnd_hrs == mut.h] <-
      #   lid.sen.data$cen[lid.sen.data$hours == mut.h & lid.sen.data$cont_hrs == ref.h]
      lid.rnd.sen.data$hours[i] <- ref.h
      lid.rnd.sen.data$rnd_hrs[i] <- mut.h
      # lid.rnd.sen.data$cen[i] <- lid.sen.data$cen[lid.sen.data$hours == mut.h & lid.sen.data$cont_hrs == ref.h]
      lid.rnd.sen.data$cen[i] <- virP1$es[virP1$ref_hrs == ref.h & virP1$que_hrs == mut.h]
      lid.rnd.sen.data$Deleterious[i] <- sum(rnd.es$lid_effect[rnd.es$hours == ref.h & rnd.es$rnd_hrs == mut.h] == 'Deleterious')
      lid.rnd.sen.data$Neutral[i] <- sum(rnd.es$lid_effect[rnd.es$hours == ref.h & rnd.es$rnd_hrs == mut.h] == 'Neutral')
      lid.rnd.sen.data$Beneficial[i] <- sum(rnd.es$lid_effect[rnd.es$hours == ref.h & rnd.es$rnd_hrs == mut.h] == 'Beneficial')
      
      mcat.rnd.sen.data$hours[i] <- ref.h
      mcat.rnd.sen.data$rnd_hrs[i] <- mut.h
      # mcat.rnd.sen.data$cen[i] <- lid.sen.data$cen[lid.sen.data$hours == mut.h & lid.sen.data$cont_hrs == ref.h]
      mcat.rnd.sen.data$cen[i] <- virP1$es[virP1$ref_hrs == ref.h & virP1$que_hrs == mut.h]
      mcat.rnd.sen.data$Deleterious[i] <- sum(rnd.es$mcat_effect[rnd.es$hours == ref.h & rnd.es$rnd_hrs == mut.h] == 'Deleterious')
      mcat.rnd.sen.data$Neutral[i] <- sum(rnd.es$mcat_effect[rnd.es$hours == ref.h & rnd.es$rnd_hrs == mut.h] == 'Neutral')
      mcat.rnd.sen.data$Beneficial[i] <- sum(rnd.es$mcat_effect[rnd.es$hours == ref.h & rnd.es$rnd_hrs == mut.h] == 'Beneficial')
      
      i <- i + 1
    # }
  }
}

lid.rnd.sen.data <- data.frame(lid.rnd.sen.data)
mcat.rnd.sen.data <- data.frame(mcat.rnd.sen.data)

lid.rnd.sen.data[,4:6] <- lid.rnd.sen.data[,4:6]/rowSums(lid.rnd.sen.data[,4:6]) * 100
mcat.rnd.sen.data[,4:6] <- mcat.rnd.sen.data[,4:6]/rowSums(mcat.rnd.sen.data[,4:6]) * 100

ggplot(lid.rnd.sen.data) +
  geom_point(aes(x = (cen-1)*100, y = Deleterious, col = 'Deleterious')) +
  geom_point(aes(x = (cen-1)*100, y = Beneficial, col = 'Beneficial')) +
  # geom_point(aes(x = (cen-1)*100, y = Neutral, col = 'Neutral')) +
  coord_cartesian(xlim = c(-50,50),
                  ylim = c(0,100))

ggplot(mcat.rnd.sen.data,
       aes(x = (cen-1)*100, y = Deleterious, col = 'Deleterious')) +
  geom_point() +
  geom_point(aes(x = (cen-1)*100, y = Beneficial, col = 'Beneficial')) +
  # geom_point(aes(x = (cen-1)*100, y = Neutral, col = 'Neutral')) +
  coord_cartesian(xlim = c(-50,50),
                  ylim = c(0,100))
  


#####


hello <- plyr::count(rnd.es, vars = c('hours','rnd_hrs','es','lid_effect'))
ggplot(hello) +
  geom_boxplot(aes(x = (es - 1) * 100, col = lid_effect, group = es))


ggplot(rnd.es,
       aes(x = (es - 1) * 100, fill = lid_result)) +
  geom_histogram(binwidth = 10) +
  labs(title = 'LID',
       x = 'Fitness Effect',
       y = '') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Neu/Ben',
                               'Ben/Ben',
                               'Del/Ben',
                               'Ben/Del',
                               'Del/Del',
                               'Neu/Del'),
                    values = c('Neu/Del'='#536DFE',
                               'Del/Del'='#303F9F',
                               'Ben/Del'='#C5CAE9',
                               'Del/Ben'='#FFECB3',
                               'Ben/Ben'='#FFA000',
                               'Neu/Ben'='#FFC107'),
                    labels = c('Neu/Del'='False-Neutral Deleterious',
                               'Del/Del'='True Deleterious',
                               'Ben/Del'='False-Beneficial Deleterious',
                               'Del/Ben'='False-Deleterious Beneficial',
                               'Ben/Ben'='True Beneficial',
                               'Neu/Ben'='False-Neutral Beneficial'),
                    drop = F) +
  theme_linedraw() +
  # facet_wrap(.~hours) +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))


ggplot(rnd.es,
       aes(x = (es - 1) * 100, fill = mcat_result)) +
  geom_histogram(binwidth = 10) +
  labs(title = 'MCAT',
       x = 'Fitness Effect',
       y = '') +
  scale_x_continuous(breaks = seq(-1000,1000,100),
                     minor_breaks = seq(-1000,1000,10),
                     labels = paste(seq(-1000,1000,100),'%', sep = '')) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Neu/Ben',
                               'Ben/Ben',
                               'Del/Ben',
                               'Ben/Del',
                               'Del/Del',
                               'Neu/Del'),
                    values = c('Neu/Del'='#536DFE',
                               'Del/Del'='#303F9F',
                               'Ben/Del'='#C5CAE9',
                               'Del/Ben'='#FFECB3',
                               'Ben/Ben'='#FFA000',
                               'Neu/Ben'='#FFC107'),
                    labels = c('Neu/Del'='False-Neutral Deleterious',
                               'Del/Del'='True Deleterious',
                               'Ben/Del'='False-Beneficial Deleterious',
                               'Del/Ben'='False-Deleterious Beneficial',
                               'Ben/Ben'='True Beneficial',
                               'Neu/Ben'='False-Neutral Beneficial'),
                    drop = F) +
  theme_linedraw() +
  # facet_wrap(.~hours) +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

######
lid.rnd.sen.data$Beneficial <- lid.rnd.sen.data$Beneficial + lid.rnd.sen.data$Deleterious
lid.rnd.sen.data$Neutral <- lid.rnd.sen.data$Neutral + lid.rnd.sen.data$Beneficial

ggplot(lid.rnd.sen.data) +
  geom_area(aes(x = (cen-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(100 * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', lwd = 1, linetype = 'dashed') +
  geom_text(x=0, y=100*0.9, label="LID", col = 'white', size = 3) +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(100 * seq(0,1,0.1)),
                     minor_breaks = c(100 * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100,
                     labels = paste0(round(seq(-2,2,0.05)*100), '%', sep = '')) +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107'),
                     guide = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(4, "mm"),
        legend.position = c(0.87,0.2),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm")) +
  coord_cartesian(xlim = c(-20,20),
                  ylim = c(0,100))



  
