
head(fdata)

fdata$orf_type[str_detect(fdata$orf_name, 'smor')] <- 'unanno pg'
fdata$orf_type[is.na(fdata$orf_type) & fdata$pg_2012 == 1] <- 'anno pg'
fdata$orf_type[is.na(fdata$orf_type)] <- 'gene'

head(fdata)

plyr::count(fdata, vars = c('orf_type','oesp1_pheno'))
plyr::count(fdata, vars = c('orf_type','oesp2_pheno'))


oesp1_pheno <- ggplot(fdata,
       aes(x = orf_type, fill = oesp1_pheno)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..),
            col = '#F5F5F5', size = 2,
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  scale_x_discrete(limits = c('unanno pg','anno pg','gene'),
                   labels = c('smorfs','annotated\nproto-genes','genes')) +
  labs(title = 'Pilot #1') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

oesp2_pheno <- ggplot(fdata,
                      aes(x = orf_type, fill = oesp2_pheno)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..),
            col = '#F5F5F5', size = 2,
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  scale_x_discrete(limits = c('unanno pg','anno pg','gene'),
                   labels = c('smorfs','annotated\nproto-genes','genes')) +
  labs(title = 'Pilot #2') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

oesp_pheno <- ggpubr::ggarrange(oesp1_pheno, oesp2_pheno,
                                nrow = 1, ncol = 2,
                                common.legend = T, legend = 'bottom')
ggsave(sprintf("%sOESP_PHENO.jpg",out_path), oesp_pheno,
       height = one.c, width = two.c, units = 'mm',
       dpi = 300)


#####
oesp1_fitness <- ggplot(fdata,
                      aes(x = orf_type, y = oesp1_cs_median)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), fill = 'transparent') +
  scale_x_discrete(limits = c('unanno pg','anno pg','gene'),
                   labels = c('smorfs','annotated\nproto-genes','genes')) +
  labs(title = 'Pilot #1',
       y = 'fitness') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm")) +
  coord_cartesian(ylim = c(0,1.2))

oesp2_fitness <- ggplot(fdata,
                        aes(x = orf_type, y = oesp2_cs_median)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), fill = 'transparent') +
  scale_x_discrete(limits = c('unanno pg','anno pg','gene'),
                   labels = c('smorfs','annotated\nproto-genes','genes')) +
  labs(title = 'Pilot #2',
       y = 'fitness') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm")) +
  coord_cartesian(ylim = c(0,1.2))

oesp_fitness <- ggpubr::ggarrange(oesp1_fitness, oesp2_fitness,
                                nrow = 1, ncol = 2,
                                common.legend = T, legend = 'bottom')
ggsave(sprintf("%sOESP_FITNESS.jpg",out_path), oesp_fitness,
       height = one.c, width = two.c, units = 'mm',
       dpi = 300)


cbind(quantile(fdata$oesp1_cs_median[fdata$orf_type == 'unanno pg'],
               c(0,0.05,0.25,0.5,0.75,0.95,1)),
       quantile(fdata$oesp2_cs_median[fdata$orf_type == 'unanno pg'],
                c(0,0.05,0.25,0.5,0.75,0.95,1)))

sort(fdata$oesp2_cs_median[fdata$orf_type == 'unanno pg' & fdata$oesp2_pheno == 'Deleterious'], decreasing = T)
