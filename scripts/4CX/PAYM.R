
out_path <- 'figs/PAYM/';

##### SPATIAL BIAS CV
load(sprintf('figs/paper/SPATIAL.RData'))

spatial$plate <- as.factor(spatial$plate)
spatial$hours <- as.factor(spatial$hours)


my_comparisons <- list(c("RAW", "LID") )

ggplot(spatial[spatial$hours == 11.04 &
                 spatial$plate == 999 &
                 (spatial$name == 'RAW' | spatial$name == 'LID'),],
       aes(x = name, y = cv, fill = plate)) +
  geom_hline(yintercept = median(spatial$cv[spatial$name == "LID" & spatial$plate == 999], na.rm = T),
             col = 'red', linetype = 'dashed') +
  geom_boxplot(outlier.colour = "black", outlier.shape = 4, outlier.size = 0.8)  +
  labs(x = '',
       y = 'CV%\n(SD/Mean * 100)') +
  scale_fill_discrete(name = '',
                      breaks = c(1,2,3,4,999),
                      labels = c('Plate 1', 'Plate 2', 'Plate 3',
                                 'Plate 4', 'All Plates'),
                      guide = F) +
  stat_compare_means(comparisons = my_comparisons) +
  # facet_wrap(.~plate) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt))

ggsave(sprintf("%sSPATIAL.jpg",out_path),
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)
