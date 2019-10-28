hello <- dbGetQuery(conn, 'select a.*, c.p, c.q, c.stat, b.mean, b.median, d.p, d.q, d.stat
from STATS_v2_CND1_SC a, STATS_v2_CND1_SC b, QVALUES_v2_CND1_SC c, QVALUES_v2_CND1_SC d
where a.exp_id = 997 and b.exp_id = 998
and c.exp_id = a.exp_id and d.exp_id = b.exp_id
and a.orf_name = b.orf_name and c.orf_name = b.orf_name and d.orf_name = c.orf_name')

ggplot(hello) +
  geom_abline(col = 'red') +
  geom_point(aes(x = median, y = median..10)) +
  labs(title = 'Archaeology Folder Set1 Vs Set2',
       subtitle = sprintf('Comparing ORFwise Fitness | r = %0.3f',cor.test(hello$median[], hello$median..10[], method="pearson")[[4]]),
       x = 'Set1',
       y = 'Set2') +
  theme_linedraw()


ggplot(hello) +
  geom_vline(xintercept = 0.01, col = 'red') +
  geom_histogram(aes(x = q , y = (..count.. * 100)/4965, fill = 'raw'), binwidth = 0.001, alpha = 0.8) +
  geom_histogram(data = hello[hello$median > 1.075692 | hello$median < 0.956022,], 
                 aes(x = q , y = (..count.. * 100)/4965, fill = 'w perc'), binwidth = 0.001, alpha = 0.6) +
  labs(title = 'Impact of PERC limits',
       subtitle = 'Set1',
       x = 'q - value',
       y = 'False Positive Rate (%)') +
  scale_fill_manual(name = '',
                    breaks = c('raw','w perc'),
                    values = c('raw' = 'grey', 'w perc' = 'blue'),
                    drop = F)

ggplot(hello) +
  geom_vline(xintercept = 0.01, col = 'red') +
  geom_histogram(aes(x = q..12 , y = (..count.. * 100)/4965, fill = 'raw'), binwidth = 0.001, alpha = 0.8) +
  geom_histogram(data = hello[hello$median..10 > 1.077052 | hello$median..10 < 0.946367,], 
                 aes(x = q..12, y = (..count.. * 100)/4965, fill = 'w perc'), binwidth = 0.001, alpha = 0.6) +
  labs(title = 'Impact of PERC limits',
       subtitle = 'Set2',
       x = 'q - value',
       y = 'False Positive Rate (%)') +
  scale_fill_manual(name = '',
                    breaks = c('raw','w perc'),
                    values = c('raw' = 'grey', 'w perc' = 'blue'),
                    drop = F)

fit.s1 <- ggplot(hello) +
  geom_histogram(aes(x = median, fill = 'ALL'), binwidth = 0.01, alpha = 0.8) +
  geom_histogram(data = hello[hello$q < 0.01,],
                 aes(x = median, fill = 'q < 0.01'), binwidth = 0.01, alpha = 0.8) +
  geom_vline(xintercept = c(1.075692, 0.956022), col = 'red', linetype = 'dashed') +
  # scale_y_continuous(trans = 'log10') +
  scale_fill_discrete(name = '') +
  labs(title = 'Fitness Distribution',
       subtitle = 'Set1',
       x = 'Fitness',
       y = 'Frequency') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0, 1.25),
                  ylim = c(0, 1000))

fit.s2 <- ggplot(hello) +
  geom_histogram(aes(x = median..10, fill = 'ALL'), binwidth = 0.01, alpha = 0.8) +
  geom_histogram(data = hello[hello$q..12 < 0.01,],
                 aes(x = median..10, fill = 'q < 0.01'), binwidth = 0.01, alpha = 0.8) +
  geom_vline(xintercept = c(1.077052, 0.946367), col = 'red', linetype = 'dashed') +
  # scale_y_continuous(trans = 'log10') +
  scale_fill_discrete(name = '') +
  labs(title = '',
       subtitle = 'Set2',
       x = 'Fitness',
       y = 'Frequency') +
  theme_linedraw() +
  coord_cartesian(xlim = c(0, 1.25),
                  ylim = c(0, 1000))

ggarrange(fit.s1, fit.s2,
          common.legend = T,
          legend = 'right')
