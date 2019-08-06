##### FIGURE 4
##### 4A
raw <- ggplot(data = fitdat, aes(x=average, col = source)) +
  geom_line(stat = 'density', lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("1TL"="#D32F2F","2TR"="#536DFE","3BL"="#388E3C","4BR"="#795548"),
                      breaks=c("1TL","2TR","3BL","4BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_x_continuous(breaks = seq(0,1000,50),
                     minor_breaks = seq(0,1000,25)) +
  scale_y_continuous(breaks = seq(0,1,0.002),
                     minor_breaks = seq(0,1,0.001),
                     limits = c(0,0.02)) +
  labs(title = 'Raw colony sizes',
       x = 'Pixel Count', y = 'Density') +
  theme_linedraw()

##### 4B
src.nrm <- ggplot(data = fitdat, aes(x=fitness, col = source)) +
  geom_line(stat = 'density', lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("1TL"="#D32F2F","2TR"="#536DFE","3BL"="#388E3C","4BR"="#795548"),
                      breaks=c("1TL","2TR","3BL","4BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_x_continuous(breaks = seq(0,2,0.1),
                     minor_breaks = seq(0,2,0.025),
                     limits = c(0.7,1.3)) +
  scale_y_continuous(breaks = seq(0,15,2),
                     minor_breaks = seq(0,15,0.5),
                     limits = c(0,13)) +
  labs(title= 'LID W Source Normalization',
       x = 'Fitness', y = '') +
  theme_linedraw()

##### 4C
no.src.nrm <- ggplot(data = fitdat, aes(x=nfitness, col = source)) +
  geom_line(stat = 'density', lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("1TL"="#D32F2F","2TR"="#536DFE","3BL"="#388E3C","4BR"="#795548"),
                      breaks=c("1TL","2TR","3BL","4BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_x_continuous(breaks = seq(0,2,0.1),
                     minor_breaks = seq(0,2,0.025),
                     limits = c(0.7,1.3)) +
  scale_y_continuous(breaks = seq(0,15,2),
                     minor_breaks = seq(0,15,0.5),
                     limits = c(0,13)) +
  labs(title= 'LID W/O Source Normalization',
       x = 'Fitness', y = '') +
  theme_linedraw()

##### BEAN
bean <- ggplot(data = alldat.b, aes(x=fitness, col = source)) +
  geom_line(stat = 'density', lwd = 1.2) + 
  scale_colour_manual(name="Source",
                      values=c("1TL"="#D32F2F","2TR"="#536DFE","3BL"="#388E3C","4BR"="#795548"),
                      breaks=c("1TL","2TR","3BL","4BR"),
                      labels=c("Top Left","Top Right","Bottom Left","Bottom Right")) +
  scale_x_continuous(breaks = seq(0,2,0.1),
                     minor_breaks = seq(0,2,0.025),
                     limits = c(0.7,1.3)) +
  scale_y_continuous(breaks = seq(0,15,2),
                     minor_breaks = seq(0,15,0.5),
                     limits = c(0,13)) +
  labs(title= 'BEAN',
       x = 'Fitness', y = '') +
  theme_linedraw()

##### FINAL FIGURE 4
fig4 <- ggarrange(raw, no.src.nrm, src.nrm, bean,
                  nrow = 1,
                  common.legend =  T,
                  legend = 'bottom')
ggsave(sprintf("%sSOURCENORM.jpg",out_path),
       fig4,
       width = 16,height = 5)


#####
plt.pow.p <- ggplot(dat.cnt2) +
  geom_area(aes(x = cen, y = Beneficial_p, fill = 'Beneficial'), alpha = 0.6) +
  # geom_point(aes(x = cen, y = Beneficial, fill = 'Beneficial'), shape = 21, col = 'black') +
  geom_area(aes(x = cen, y = Deleterious_p, fill = 'Deleterious'), alpha = 0.6) +
  # geom_point(aes(x = cen, y = Deleterious, fill = 'Deleterious'), shape = 21, col = 'black') +
  geom_area(aes(x = cen, y = Neutral_p, fill = 'Neutral'), alpha = 0.6) +
  # geom_point(aes(x = cen, y = Neutral, fill = 'Neutral'), shape = 21, col = 'black') +
  labs(title = "Whats the Power to Detect Small Effects?",
       subtitle = 'at p <= 0.05',
       x = 'Effect Size',
       y = 'Percent Data') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.01)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.01)) +
  scale_color_manual(name = '',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#D32F2F',
                                'Neutral'='#303F9F',
                                'Beneficial'='#4CAF50')) +
  scale_fill_manual(name = '',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#D32F2F',
                               'Neutral'='#303F9F',
                               'Beneficial'='#4CAF50')) +
  theme_linedraw() +
  theme(legend.position = 'bottom') +
  coord_cartesian(xlim = c(0.8,1.2),
                  ylim = c(0,t))

stats.tmp$hours <- as.character(stats.tmp$hours)
stats.tmp$cont_hrs <- as.character(stats.tmp$cont_hrs)

plt.fpr.z <- ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,]) +
  geom_abline(col = 'red', linetype = 'dashed', lwd = 1) +
  geom_line(aes(x = p, y = fpr, col = hours), lwd = 1.2) +
  labs(title = "Does FPR Follow Random Expectation?",
       subtitle = "",
       x = "p-value",
       y = "False Positive Rate") +
  scale_x_continuous(breaks = seq(0,1,0.025),
                     minor_breaks = seq(0,1,0.00625),
                     limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,0.025),
                     minor_breaks = seq(0,1,0.00625),
                     limits = c(0,1)) +
  scale_color_manual(name = "Hours",
                     breaks=c("13","14","16","17","18"),
                     values=c("13"="#D32F2F","14"="#536DFE","16"="#388E3C","17"="#795548","18"="#00BCD4",
                              "0"="transparent","8"="transparent","9"="transparent","10"="transparent","11"="transparent")) +
  theme_linedraw() +
  theme(legend.position = 'bottom') +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 0.1))

edis <- annotate_figure(ggarrange(plt.fpr.z, plt.pow.p,
                  nrow = 1, ncol = 2),
                  top = 'LID')

ggsave(sprintf("%s%s_EDIS_%d.jpg",out_path,expt_name,rep),
       height = 6, width = 10,
       dpi = 300)

##### MESS UP PLATES
oridata$from[is.na(oridata$from)] = 'ref'

ggplot(oridata[oridata$hours == 14 & oridata$`6144plate` == 1,]) +
  geom_point(aes(x = `6144col`, y = `6144row`, col = from)) +
  scale_x_continuous(breaks = seq(1,96,2),limits = c(1,96)) +
  scale_y_continuous(breaks = seq(1,64,2),limits = c(64,1),trans = 'reverse') +
  scale_color_manual(name = '',
                     breaks = c('ref', 'less', 'more'),
                     values = c('ref' = '#FFC107', 'less' = '#FF5252', 'more' = '#388E3C'),
                     labels = c('Ref', '<14 Hrs', '>14 Hrs')) +
  labs(title = 'Example Messed Up Plate',
       subtitle = 'T (R) = 14 Hours | Plate 1',
       x = 'Columns',
       y = 'Rows') +
  theme_linedraw()
ggsave(sprintf("%sMESSEDUP.jpg",out_path),
       height = 7, width = 10,
       dpi = 300)
