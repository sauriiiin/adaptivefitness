
load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/solid_results.RData')
load("/home/sbp29/R/Projects/adaptivefitness/figs/f28fu/repeat/f28fu_data.RData")

f28.res <- solid_fit[(solid_fit$density == 1536 & solid_fit$hours == 48 & solid_fit$arm != 'SDA') |
                       (solid_fit$density == 1536 & solid_fit$hours == 72 & solid_fit$arm == 'SDA') |
                       (solid_fit$density == 6144 & solid_fit$hours == 24 & solid_fit$arm != 'SDA') |
                       (solid_fit$density == 6144 & solid_fit$hours == 36 & solid_fit$arm == 'SDA'),]
head(f28.res)
foa.res <- foa.data[((foa.data$density == 1536 & foa.data$hours == 60) |
                       (foa.data$density == 6144 & foa.data$hours == 36)) &
                      !(foa.data$orf_name %in% c('BOR','REF','YHR021W-A')),]
foa.res$condition <- factor(foa.res$condition, levels = c('GLU','GAL','CAS','SDA'))

f28.res.plot <- ggplot(f28.res) +
  geom_boxplot(aes(x = orf_name, y = fitness, fill = orf_name),
               outlier.shape = NA) +
  facet_wrap(.~density*arm) +
  scale_y_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.01)) +
  coord_cartesian(ylim = c(0.8,1.2)) +
  scale_fill_discrete(guide = F) +
  labs(title = 'F28 Results',
       x = 'Strains', y = 'Fitness') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))


foa.res.plot <- ggplot(foa.res[foa.res$condition!= 'GAL' & !is.na(foa.res$orf_name),]) +
  geom_boxplot(aes(x = orf_name, y = fitness_clean, fill = orf_name),
               outlier.shape = NA) +
  facet_wrap(.~density*condition) +
  scale_y_continuous(breaks = seq(0,2,0.05),
                     minor_breaks = seq(0,2,0.01)) +
  coord_cartesian(ylim = c(0.8,1.2)) +
  scale_fill_discrete(guide = F) +
  labs(title = '5FOA Results',
       x = 'Strains', y = 'Fitness') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))


res.res.f28foa <- ggarrange(f28.res.plot, foa.res.plot, nrow = 1)
ggsave(sprintf("%sBOXPLOT_F28_FOA.jpg",out_path), res.res.f28foa,
       height = two.c, width = two.c*2, units = 'mm',
       dpi = 300)




##### SCATTER
foa.res.m <- foa.res[!is.na(foa.res$orf_name),] %>% 
  group_by(orf_name, condition, density) %>%
  summarise(fitness = median(fitness_clean, na.rm = T)) %>%
  data.frame()

f28.res.m <- f28.res %>%
  group_by(orf_name, arm, density) %>%
  summarise(fitness = median(fitness, na.rm = T)) %>%
  data.frame

hello <- merge(foa.res.m, f28.res.m, by.x = c("orf_name","density","condition"),
      by.y = c("orf_name","density","arm"))

ggplot(hello, aes(x = fitness.x, y = fitness.y)) +
  # geom_abline() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman') +
  geom_point() +
  labs(x = 'Fitness without Plasmid', y = 'Fitness with Plasmid') +
  facet_wrap(.~density*condition) +
  coord_cartesian(xlim = c(0.85, 1.10),
                  ylim = c(0.85, 1.10)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sSCATTER_F28_FOA.jpg",out_path),
       height = two.c*2/3, width = two.c, units = 'mm',
       dpi = 300)

##### FOA DENSITY COMPARE
hello2

merge(foa.res.m[foa.res.m$density == 1536,], foa.res.m[foa.res.m$density == 6144,],
      by = c('orf_name','condition')) %>%
  ggplot(aes(x = fitness.x, y = fitness.y)) +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman') +
  geom_point() +
  labs(x = 'Fitness at 1546', y = 'Fitness at 6144') +
  facet_wrap(.~condition) +
  coord_cartesian(xlim = c(0.9, 1.10),
                  ylim = c(0.9, 1.10)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sSCATTER_1536_6144.jpg",out_path),
       height = two.c*2/3, width = two.c*2/3, units = 'mm',
       dpi = 300)

