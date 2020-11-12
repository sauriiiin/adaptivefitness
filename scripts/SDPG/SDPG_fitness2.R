
fit_stat2 <- data.frame(compare_means(fitness ~ orf_name, lidres[lidres$hours >= 8,],
                                     method = "wilcox.test", paired = FALSE,
                                     p.adjust.method = 'bonferroni',
                                     group.by = c("density", "arm", "rep", "hours"), ref.group = "BF_control"))
head(fit_stat2)

for (a in unique(fit_stat2$arm)) {
  for (d in unique(fit_stat2$density[fit_stat2$arm == a])) {
    for (r in unique(fit_stat2$rep[fit_stat2$density == d & fit_stat2$arm == a])) {
      for (h in unique(fit_stat2$hours[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r])) {
        cont_fit <- lidres$fitness[lidres$density == d & lidres$arm == a & lidres$rep == r &
                                     lidres$hours == h & lidres$orf_name == 'BF_control']
        cont_mean <- median(lidres$fitness[lidres$density == d & lidres$arm == a & lidres$rep == r &
                                             lidres$hours == h & lidres$orf_name == 'BF_control'], na.rm = T)
        for (o in unique(fit_stat2$group2[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r &
                                         fit_stat2$hours == h & fit_stat2$group2 != 'BF_control'])) {
          orf_fit <- lidres$fitness[lidres$density == d & lidres$arm == a & lidres$rep == r &
                                      lidres$hours == h & lidres$orf_name == o]
          orf_mean <- median(lidres$fitness[lidres$density == d & lidres$arm == a & lidres$rep == r &
                                              lidres$hours == h & lidres$orf_name == o], na.rm = T)
          temp_cd <- cliff.delta(orf_fit, cont_fit, conf.level=.95,
                                 use.unbiased=TRUE, use.normal=FALSE, return.dm=FALSE)
          fit_stat2$empirical_p[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r &
                                 fit_stat2$hours == h & fit_stat2$group2 == o] <- lidres$p[lidres$density == d & lidres$arm == a & lidres$rep == r &
                                                                                           lidres$hours == h & lidres$orf_name == o][1]
          fit_stat2$cliff.delta[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r &
                                 fit_stat2$hours == h & fit_stat2$group2 == o] <- temp_cd$estimate
          fit_stat2$magnitude[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r &
                               fit_stat2$hours == h & fit_stat2$group2 == o] <- as.character(temp_cd$magnitude)
          fit_stat2$effect_size[fit_stat2$density == d & fit_stat2$arm == a & fit_stat2$rep == r &
                                 fit_stat2$hours == h & fit_stat2$group2 == o] <- (orf_mean - cont_mean)/cont_mean * 100
          
        }
      }
    }
  }
}

fit_stat2$phenotype[fit_stat2$p.signif == 'ns'] <- 'Neutral'
fit_stat2$phenotype[fit_stat2$p.signif != 'ns' & fit_stat2$cliff.delta < 0] <- 'Deleterious'
fit_stat2$phenotype[fit_stat2$p.signif != 'ns' & fit_stat2$cliff.delta > 0] <- 'Beneficial'

fit_stat2$empirical_phenotype[fit_stat2$empirical_p > 0.05] <- 'Neutral'
fit_stat2$empirical_phenotype[fit_stat2$empirical_p <= 0.05 & fit_stat2$cliff.delta < 0] <- 'Deleterious'
fit_stat2$empirical_phenotype[fit_stat2$empirical_p <= 0.05 & fit_stat2$cliff.delta > 0] <- 'Beneficial'

fit_stat2$saturation <- NA
fit_stat2$saturation[fit_stat2$density == 1536 & fit_stat2$hours == 48 & fit_stat2$arm != 'SDA'] <- "Saturated"
fit_stat2$saturation[fit_stat2$density == 1536 & fit_stat2$hours %in% c(72, 78) & fit_stat2$arm == 'SDA'] <- "Saturated"
fit_stat2$saturation[fit_stat2$density == 6144 & fit_stat2$hours == 24 & fit_stat2$arm != 'SDA'] <- "Saturated"
fit_stat2$saturation[fit_stat2$density == 6144 & fit_stat2$hours == 36 & fit_stat2$arm == 'SDA'] <- "Saturated"
fit_stat2$saturation[is.na(fit_stat2$saturation)] <- "Other"

fit_stat2$arm <- factor(fit_stat2$arm , levels = c('GLU','CAS','SDA'))

head(fit_stat2)

fit_stat2 <- fit_stat2[fit_stat2$hours > 10,]

#### PHENO 2: MATRIX
pheno_mat <- ggplot(fit_stat2,
                    aes(x = group2, y = name, fill = phenotype), col = 'black') +
  geom_tile() +
  labs(x = "ORFs",
       y = "") +
  scale_fill_manual(name = '',
                    breaks = c("1536","6144",
                               "GLU","CAS","SDA",
                               "R1","R2","R3",
                               'Deleterious','Neutral','Beneficial'),
                    values = c("1536" = "#607D8B", "6144" = "#FFA000",
                               "GLU" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F",
                               "R1" = "#E1BEE7", "R2" = "#E040FB", "R3" = "#7B1FA2",
                               'Deleterious'='#3F51B5','Neutral'='#212121','Beneficial'='#FFC107'),
                    drop = F) +
  scale_y_discrete(limits = unique(fit_stat2$name[order(fit_stat2$arm)])) +
  theme_linedraw() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = txt, angle = 90),
        axis.title = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"))
pheno_leg <- get_legend(pheno_mat)

axis_y <- ggplot(fit_stat2) +
  geom_tile(aes(y = name, x = 'Condition', fill = arm), col = 'black') +
  geom_tile(aes(y = name, x = 'Replicate', fill = rep), col = 'black') +
  geom_tile(aes(y = name, x = 'Density', fill = as.character(density)), col = 'black') +
  geom_point(aes(y = name, x = 'Hours', col = as.character(hours))) +
  # geom_text(aes(y = name, x = 'Replicate', label = rep), col = 'black', size = 3, angle = 90) +
  # geom_text(aes(y = name, x = 'Condition', label = arm), col = 'black', size = 3, angle = 90) +
  # geom_text(aes(y = name, x = 'Density', label = density), col = 'black', size = 3, angle = 90) +
  labs(x = '', y = '') +
  scale_x_discrete(breaks = c('Condition','Density','Replicate','Hours'),
                   limits = c('Condition','Density','Replicate','Hours')) +
  scale_y_discrete(limits = unique(fit_stat2$name[order(fit_stat2$arm)])) +
  scale_fill_manual(name = '',
                    breaks = c("1536","6144",
                               "GLU","CAS","SDA",
                               "R1","R2","R3",
                               'Deleterious','Neutral','Beneficial'),
                    values = c("1536" = "#607D8B", "6144" = "#FFA000",
                               "GLU" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F",
                               "R1" = "#E1BEE7", "R2" = "#E040FB", "R3" = "#7B1FA2",
                               'Deleterious'='#3F51B5','Neutral'='#212121','Beneficial'='#FFC107'),
                    drop = F) +
  scale_color_discrete(name = '',
                       breaks = as.character(sort(unique(fit_stat2$hours)))) +
  theme_linedraw() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = txt, angle = 90),
        axis.title = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm")) +
  guides(col = guide_legend(ncol=2))
axis_leg <- get_legend(axis_y)

fit_pheno2 <- ggpubr::ggarrange(ggpubr::ggarrange(axis_y, pheno_mat, nrow = 1, widths = c(1,6),
                                                  align = 'hv', common.legend = T,legend = 'none'),
                                ggpubr::ggarrange(pheno_leg, axis_leg, heights = c(1,3), nrow = 2),
                                nrow = 1, widths = c(5,1))
ggsave(sprintf("%sFITNESS_PHENO2.jpg",out_path), fit_pheno2,
       height = 175, width = 200, units = 'mm',
       dpi = 300)


#### PHENO 3: DYNAMICS
fit_pheno3 <- ggplot(fit_stat2,
       aes(x = hours, y = cliff.delta)) +
  geom_hline(yintercept = c(-0.474,0.474), col = '#212121', linetype = 'dashed', lwd = 0.3) +
  geom_hline(yintercept = c(-0.33,0.33), col = '#757575', linetype = 'dashed', lwd = 0.3) +
  geom_hline(yintercept = c(-0.147,0.147), col = '#BDBDBD', linetype = 'dashed', lwd = 0.3) +
  geom_line(aes(group = 1, col = rep)) +
  geom_point(aes(col = phenotype)) +
  facet_wrap(.~arm*density*group2*rep,nrow = 30,
             scale = 'free_x') +
  scale_color_manual(name = '',
                    breaks = c('Deleterious','Neutral','Beneficial',
                               "R1","R2","R3"),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107',
                               "R1" = "#E1BEE7", "R2" = "#E040FB", "R3" = "#7B1FA2")) +
  labs(x = 'Time (hours)',
       y = "Cliff's Delta") +
  # scale_y_discrete(limits = c('Deleterious','Neutral','Beneficial')) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFITNESS_PHENO3.jpg",out_path),fit_pheno3,
       height = 1500, width = 1000, units = 'mm', limitsize = F,
       dpi = 300)


##### PHENO 4
sand_pheno <- dbGetQuery(conn, "select orf_name, exp_1 GLU, exp_91 CAS, exp_93 SDA
                         from SAND_PHENOTYPE
                         where (exp_28 = 'beneficial'
                         or exp_31 = 'beneficial'
                         or exp_33 = 'beneficial'
                         or exp_91 = 'beneficial'
                         or exp_93 = 'beneficial')
                         and orf_name in
                         (select orf_name from PROTOGENES
                         where pg_2012 = 1)
                         and orf_name != 'YHR021W-A'")
sand_pheno <- melt(sand_pheno, id.vars = 'orf_name', variable.name = 'arm', value.name = 'phenotype')
sand_pheno$density <- 6144
sand_pheno$phenotype[sand_pheno$arm == 'GLU'] <- "Neutral"
sand_pheno$arm <- factor(sand_pheno$arm, levels = c('GLU','CAS','SDA'))

fit_pheno4_mat <- ggplot(fit_stat2[fit_stat2$saturation == 'Saturated',]) +
  geom_tile(aes(x = rep, y = group2, fill = phenotype)) +
  labs(x = 'Replicates',
       y = "ORFs") +
  facet_wrap(.~arm*density, nrow = 3) +
  scale_fill_manual(name = '',
                     breaks = c('Deleterious','Neutral','Beneficial'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

fit_pheno4_cnt <- ggplot(plyr::count(fit_stat2[fit_stat2$saturation == 'Saturated',],
                   vars = c('group2', 'arm', 'density', 'phenotype'))) +
  geom_text(aes(x = phenotype, y = group2, label = freq, col = phenotype), size = 3) +
  labs(x = 'Phenotype',
       y = "ORFs") +
  scale_x_discrete(limits = c('Deleterious','Neutral','Beneficial')) +
  facet_wrap(.~arm*density, nrow = 3) +
  scale_color_manual(name = '',
                    breaks = c('Deleterious','Neutral','Beneficial'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    guide = F) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        panel.background = element_rect(fill = 'grey60'),
        axis.title = element_text(size = titles),
        axis.title.y = element_blank(),
        # axis.text = element_text(size = txt),
        axis.text = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

fit_pheno4_sand <- ggplot(sand_pheno) +
  geom_text(aes(x = 1, y = orf_name, label = value, col = value), size = 3) +
  facet_wrap(.~variable*density, nrow = 3) +
  scale_color_manual(name = '',
                     breaks = c('deleterious','neutral','beneficial'),
                     values = c('deleterious'='#3F51B5',
                                'neutral'='#212121',
                                'beneficial'='#FFC107'),
                     guide = F) +
  labs(x = 'Expected Phenotype',
       y = "ORFs") +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        panel.background = element_rect(fill = 'grey60'),
        axis.title = element_text(size = titles),
        axis.title.y = element_blank(),
        # axis.text = element_text(size = txt),
        axis.text = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
  
fit_pheno4 <- ggpubr::ggarrange(fit_pheno4_mat, fit_pheno4_cnt, fit_pheno4_sand,
                                nrow = 1,
                                align = 'h',
                                widths = c(4,1,0.5),
                                common.legend = T,
                                legend = 'bottom')
ggsave(sprintf("%sFITNESS_PHENO4.jpg",out_path),fit_pheno4,
       height = 250, width = 250, units = 'mm',
       dpi = 300)


#####
col_pheno <- plyr::count(fit_stat2[fit_stat2$saturation == 'Saturated',], vars = c('group2', 'arm' ,'phenotype'))

phenos <- NULL
i <- 1
for (a in unique(col_pheno$arm)) {
  for (o in unique(col_pheno$group2[col_pheno$arm == a])) {
    if (dim(col_pheno[col_pheno$arm == a & col_pheno$group2 == o & col_pheno$phenotype != 'Neutral',])[1] > 0) {
      phenos$arm[i] <- a
      phenos$orf_name[i] <- o
      phenos$col_pheno[i] <- col_pheno$phenotype[col_pheno$arm == a & col_pheno$group2 == o &
                                                   col_pheno$phenotype != 'Neutral' & sort(col_pheno$freq)][1]
      phenos$liq_pheno[i] <- liq_stats$phenotype[liq_stats$Arm == a & liq_stats$group2 == o]
      i <- i + 1
    } else {
      phenos$arm[i] <- a
      phenos$orf_name[i] <- o
      phenos$col_pheno[i] <- 'Neutral'
      phenos$liq_pheno[i] <- liq_stats$phenotype[liq_stats$Arm == a & liq_stats$group2 == o]
      i <- i + 1
    }
  }
}
phenos <- data.frame(phenos)
phenos$arm <- factor(phenos$arm, levels = c('GLU','CAS','SDA'))

fit_pheno5 <- ggplot(phenos) +
  geom_text(aes(x = 'Solid', y = orf_name, col = col_pheno, label = col_pheno), size = 3) +
  geom_text(aes(x = 'Liquid', y = orf_name, col = liq_pheno, label = liq_pheno), size = 3) +
  facet_wrap(.~arm, nrow = 3) +
  labs(x = 'Media',
       y = 'ORFs') +
  scale_x_discrete(limits = c('Solid','Liquid')) +
  scale_color_manual(name = '',
                     breaks = c('Deleterious','Neutral','Beneficial'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107'),
                     guide = F) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        panel.background = element_rect(fill = 'grey60'),
        axis.title = element_text(size = titles),
        axis.title.y = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%sFITNESS_PHENO5.jpg",out_path), fit_pheno5,
#        height = 250, width = 100, units = 'mm',
#        dpi = 300)

#####
fit_pheno6 <- ggpubr::ggarrange(fit_pheno4_mat, fit_pheno4_cnt, fit_pheno5, fit_pheno4_sand,
                                nrow = 1,
                                align = 'h',
                                widths = c(3,1,1,0.5),
                                common.legend = T,
                                legend = 'bottom')
ggsave(sprintf("%sFITNESS_PHENO6.jpg",out_path),fit_pheno6,
       height = 300, width = 300, units = 'mm',
       dpi = 300)

##### COUNTS
fit_pheno7 <- ggplot(phenos) +
  geom_tile(aes(x = arm, y = orf_name, fill = col_pheno), col = 'black') +
  labs(title = 'Solid Phenotype',
       x = 'Condition',
       y = 'ORFs') +
  scale_fill_manual(name = '',
                     breaks = c('Deleterious','Neutral','Beneficial'),
                     values = c('Deleterious'='#3F51B5',
                                'Neutral'='#212121',
                                'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFITNESS_PHENO7.jpg",out_path),fit_pheno7,
       height = 150, width = 150, units = 'mm',
       dpi = 300)

liq_pheno2 <- ggplot(phenos) +
  geom_tile(aes(x = arm, y = orf_name, fill = liq_pheno), col = 'black') +
  labs(title = 'Liquid Phenotype',
       x = 'Condition',
       y = 'ORFs') +
  scale_fill_manual(name = '',
                    breaks = c('Deleterious','Neutral','Beneficial'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sLIQUID_PHENO2.jpg",out_path),liq_pheno2,
       height = 150, width = 150, units = 'mm',
       dpi = 300)

#####
common_pheno <- matrix(c(18, 9, 16, 11), nrow = 2,
                       dimnames = list(c("Same", "Different"),
                                       c("Pitt", "SD")))
fisher.test(common_pheno)

##### SAVING SOLID RESULTS
solid_fit <- lidres
colnames(solid_fit) <- c("pos", "density", "plate", "row", "col", "strain_id", "orf_name",
                         "hours", "average", "fitness", "N", "fitness_mean", "fitness_median", "fitness_std",
                         "empirical_p", "es",
                         "arm", "replicate", "normalized_cs", "saturation")
head(solid_fit)
solid_stats <- fit_stat2
colnames(solid_stats) <- c("density", "arm", "replicate", "hours", "variable", "reference", "orf_name",
                           "p", "p.adj", "p.format", "p.signif", "method",
                           "empirical_p",
                           "cliff.delta", "magnitude",
                           "es", "phenotype", "empirical_phenotype", "saturation")
solid_stats <- merge(solid_stats, sand_pheno[,1:3], by = c("orf_name","arm","density"), all = T)
colnames(solid_stats) <- c("orf_name", "arm", "density", "replicate", "hours", "variable", "reference",
                           "p", "p.adj", "p.format", "p.signif", "method",
                           "empirical_p",
                           "cliff.delta", "magnitude",
                           "es", "pitt_phenotype", "empirical_phenotype", "saturation",
                           "sandiego_phenotype")

head(solid_stats[,-6])
solid_stats <- solid_stats[,-6]

solid_stats$sandiego_phenotype[solid_stats$sandiego_phenotype == 'beneficial' &
                                 !is.na(solid_stats$sandiego_phenotype)] <- 'Beneficial'
solid_stats$sandiego_phenotype[solid_stats$sandiego_phenotype == 'neutral' &
                                 !is.na(solid_stats$sandiego_phenotype)] <- 'Neutral'
solid_stats$sandiego_phenotype[solid_stats$sandiego_phenotype == 'deleterious' &
                                 !is.na(solid_stats$sandiego_phenotype)] <- 'Deleterious'

save(solid_fit, solid_stats,
     file = 'figs/SDPG/solid_results.RData')
