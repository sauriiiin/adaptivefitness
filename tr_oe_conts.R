

kapitzky <- read.csv(file = 'rawdata/translatome/kapitzky_data.csv')
head(kapitzky)

fig.conts.expt <- merge(merge(data.mad, melt(controls, id.vars = c('strain_id','standard_name','orf_name'),
                     variable.name = 'condition', value.name = 'control_type'), by = c('condition','strain_id','orf_name')) %>%
  filter(hours %in% c(141, 36, 20, 125, 26), control_type != '') %>%
  group_by(condition, strain_id, orf_name, standard_name, control_type, fitness_ll, fitness_ul) %>%
  summarize(fitness.median = median(fitness.median, na.rm = T), .groups = 'keep') %>%
  data.frame(), melt(kapitzky,  id.vars = c('strain_id','standard_name','orf_name'),
                     variable.name = 'condition', value.name = 'fitness_expt'), 
  by = c('condition','strain_id','orf_name','standard_name')) %>%
  ggplot(aes(x = fitness.median, y = fitness_expt)) +
  geom_vline(aes(xintercept = fitness_ll), size = 0.5, linetype = 'dashed', col = 'red') +
  geom_vline(aes(xintercept = fitness_ul), size = 0.5, linetype = 'dashed', col = 'red') +
  geom_point(aes(col = control_type)) +
  geom_text_repel(aes(label = standard_name), size = 2) +
  labs(x = 'Fitness Observed',
       y = 'Fitness Expected') +
  scale_color_manual(name = 'Control Type',
                     values = c('Resistant' = '#FFA000',
                                'Sensitive' = '#303F9F')) +
  facet_wrap(.~condition, nrow = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/CONDITIONCONTROLS_EXPT.jpg",fig_path), fig.conts.expt,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)
