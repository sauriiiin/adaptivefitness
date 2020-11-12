##### F28 RESULTS
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 11/12/2020

##### INITIALIZE
library(ggplot2)
library(ggridges)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(scales)
library(egg)
library(zoo)
library(ggrepel)
library(reshape2)

out_path = 'figs/SDPG/ALL/';
out_path_gc <- 'figs/SDPG/GC/'

##### FIGURE TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### LOAD DATA
load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/solid_results.RData')
load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/liquid_results.RData')
liquid_results <- all_results
liquid_summary <- all_results2
colnames(liquid_summary) <- c("arm", "replicate", "orf_name", "bio_rep", "method", "value", "liquid_phenotype")

##### COMBINED PHENOTYPES
liquid_phenotype <- liquid_results[,c(1,2,3,seq(1,ncol(liquid_results))[str_detect(names(liquid_results), 'phenotype')])]
liquid_phenotype <- data.frame(liquid_phenotype %>%
                                 group_by(condition, exp_rep, orf_name) %>%
                                 summarise(aleeza_phenotype = max(aleeza_phenotype, na.rm = T),
                                           aleeza_plc_phenotype = max(alleza_plc_phenotype, na.rm = T),
                                           growthrates_phenotype = max(growthrates_phenotype, na.rm = T),
                                           growthrates_plc_phenotype = max(growthrates_plc_phenotype, na.rm = T),
                                           growthcurver_phenotype = max(growthcurver_phenotype, na.rm = T),
                                           growthcurver_plc_phenotype = max(growthcurver_plc_phenotype, na.rm = T),
                                           growthcurver_auc_phenotype = max(growthcurver_auc_phenotype, na.rm = T),
                                           growthcurver_plc_auc_phenotype = max(growthcurver_plc_auc_phenotype, na.rm = T),
                                           growthcurver_auc2_phenotype = max(growthcurver_auc2_phenotype, na.rm = T),
                                           growthcurver_plc_auc2_phenotype = max(growthcurver_plc_auc2_phenotype, na.rm = T)))
colnames(liquid_phenotype) <- c("arm", "replicate", "orf_name", "aleeza_phenotype", "aleeza_plc_phenotype",
                                "growthrates_phenotype", "growthrates_plc_phenotype",
                                "growthcurver_phenotype", "growthcurver_plc_phenotype",
                                "growthcurver_auc_phenotype", "growthcurver_plc_auc_phenotype",
                                "growthcurver_auc2_phenotype", "growthcurver_plc_auc2_phenotype")
head(liquid_phenotype)
solid_1536_phenotype <- solid_stats[solid_stats$density == 1536 & solid_stats$saturation == 'Saturated',c(1,2,4,16,19)]
solid_6144_phenotype <- solid_stats[solid_stats$density == 6144 & solid_stats$saturation == 'Saturated',c(1,2,4,16)]
solid_phenotype <- merge(solid_1536_phenotype, solid_6144_phenotype, by = c('orf_name','arm','replicate'),
                         suffixes = c('_1536','_6144'), all = T)

all_phenotype <- merge(liquid_phenotype, solid_phenotype, by = c('orf_name','arm','replicate'), all = T)
all_phenotype <- all_phenotype[,c(1:14,16,15)]
head(all_phenotype)

pheno_mat <- ggplot(melt(all_phenotype, id.vars = c('orf_name','arm','replicate'),
                         variable.name = 'method', value.name = 'phenotype'),
                    aes(x = method, y = orf_name)) +
  geom_tile(aes(fill = phenotype), col = 'black') +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Deleterious','Neutral','Beneficial'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    na.value = 'grey50') +
  scale_x_discrete(limits = names(all_phenotype[4:ncol(all_phenotype)]),
                   labels = str_remove_all(names(all_phenotype[4:ncol(all_phenotype)]), '_phenotype')) +
  labs(x = 'Method',
       y = 'ORFs') +
  facet_wrap(.~arm*replicate) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sALL_PHENOTYPE_MATRIX.png',out_path), pheno_mat,
       height = 300, width = 300, units = 'mm',
       dpi = 300, limitsize = F)


##### SOLID V LIQUID

most_common_phenotype <- NULL
i <- 1
for (a in unique(all_phenotype$arm)) {
  for (o in unique(all_phenotype$orf_name[all_phenotype$arm == a])) {
    for (n in names(all_phenotype[all_phenotype$arm == a & all_phenotype$orf_name == o,])[c(11,14:16)]){
      most_common_phenotype$arm[i] <- a
      most_common_phenotype$orf_name[i] <- o
      most_common_phenotype$method[i] <- str_remove(n, '_phenotype')
      if (length(names(which.max(table(all_phenotype[all_phenotype$arm == a & all_phenotype$orf_name == o,n])))) > 0) {
        most_common_phenotype$phenotype[i] <- names(which.max(table(all_phenotype[all_phenotype$arm == a & all_phenotype$orf_name == o,n])))
        i <- i + 1
      } else {
        i <- i + 1
      }
    }
  }
}
most_common_phenotype <- data.frame(most_common_phenotype)
most_common_phenotype$arm <- factor(most_common_phenotype$arm, levels = c('GLU','CAS','SDA'))
common_pheno_mat <- ggplot(most_common_phenotype,
       aes(x = method, y = orf_name)) +
  geom_tile(aes(fill = phenotype), col = 'black') +
  facet_wrap(.~arm) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Deleterious','Neutral','Beneficial'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    na.value = 'grey50') +
  labs(x = 'Method',
       y = 'ORFs') +
  facet_wrap(.~arm, nrow = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sCOMMON_PHENOTYPE_MATRIX.png',out_path), common_pheno_mat,
       height = 120, width = 300, units = 'mm',
       dpi = 300, limitsize = F)


###### DENSITY V DENSITY
solid_1536_stats <- solid_stats[solid_stats$density == 1536 & solid_stats$saturation == 'Saturated',c(1,2,4,12,13)]
solid_6144_stats <- solid_stats[solid_stats$density == 6144 & solid_stats$saturation == 'Saturated',c(1,2,4,12,13)]
solid_cliff <- merge(solid_1536_stats, solid_6144_stats, by = c('arm','replicate','orf_name'),
                     suffixes = c('_1536','_6144'))

cor.test(solid_cliff$cliff.delta_1536, solid_cliff$cliff.delta_6144, method = 'pearson')

dd_corr <- ggplot(solid_cliff,
       aes(x = cliff.delta_6144, y = cliff.delta_1536)) +
  geom_abline() +
  # geom_vline(xintercept = c(-0.474,-0.33,-0.147,0.147,0.33,0.474)) +
  # geom_hline(yintercept = c(-0.474,-0.33,-0.147,0.147,0.33,0.474)) +
  geom_point(aes(col = arm, shape = replicate)) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = 'spearman') +
  labs(x = "Cliff's Delta at 6144 Density",
       y = "Cliff's Delta at 1536 Density") +
  coord_cartesian(xlim = c(-1,1),
                  ylim = c(-1,1)) +
  scale_color_discrete(name = 'Condition') +
  scale_shape_discrete(name = 'expRep') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sDENSITYDENSITY_CLIFF_CORR.png',out_path), dd_corr,
       height = 120, width = 130, units = 'mm',
       dpi = 300, limitsize = F)


dd_corr2 <- ggplot(solid_cliff,
                  aes(x = cliff.delta_6144, y = cliff.delta_1536)) +
  geom_abline() +
  # geom_vline(xintercept = c(-0.474,-0.33,-0.147,0.147,0.33,0.474)) +
  # geom_hline(yintercept = c(-0.474,-0.33,-0.147,0.147,0.33,0.474)) +
  geom_point() +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = 'spearman') +
  labs(x = "Cliff's Delta at 6144 Density",
       y = "Cliff's Delta at 1536 Density") +
  coord_cartesian(xlim = c(-1,1),
                  ylim = c(-1,1)) +
  scale_color_discrete(name = 'Condition') +
  scale_shape_discrete(name = 'expRep') +
  theme_linedraw() +
  facet_wrap(.~arm*replicate) +
  theme(plot.title = element_text(size = titles + 2,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sDENSITYDENSITY_CLIFF_CORR2.png',out_path), dd_corr2,
       height = 120, width = 175, units = 'mm',
       dpi = 300, limitsize = F)

##### REPLICATE V REPLICATE
solid_r1_cliff <- solid_stats[solid_stats$replicate == 'R1' & solid_stats$saturation == 'Saturated',c(1:3,13)]
solid_r2_cliff <- solid_stats[solid_stats$replicate == 'R2' & solid_stats$saturation == 'Saturated',c(1:3,13)]
solid_r3_cliff <- solid_stats[solid_stats$replicate == 'R3' & solid_stats$saturation == 'Saturated',c(1:3,13)]

solid_rep_cliff1 <- merge(solid_r1_cliff, solid_r2_cliff, by = c('arm','density','orf_name'),
                     suffixes = c('_R1','_R2'), all = T)
solid_rep_cliff <- merge(solid_rep_cliff1, solid_r3_cliff, by = c('arm','density','orf_name'), all = T)

colnames(solid_rep_cliff) <- c('arm','density','orf_name','R1','R2','R3')

rr12 <- ggplot(solid_rep_cliff,
       aes(x = R1, y = R2)) +
  geom_abline() +
  geom_point() +
  facet_wrap(.~arm*density, nrow = 1) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = 'spearman') +
  labs(x = "Cliff's Delta in R1",
       y = "Cliff's Delta in R2") +
  coord_cartesian(xlim = c(-1,1),
                  ylim = c(-1,1)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

rr13 <- ggplot(solid_rep_cliff,
               aes(x = R1, y = R3)) +
  geom_abline() +
  geom_point() +
  facet_wrap(.~arm*density, nrow = 1) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = 'spearman') +
  labs(x = "Cliff's Delta in R1",
       y = "Cliff's Delta in R3") +
  coord_cartesian(xlim = c(-1,1),
                  ylim = c(-1,1)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

rr23 <- ggplot(solid_rep_cliff,
               aes(x = R2, y = R3)) +
  geom_abline() +
  geom_point() +
  facet_wrap(.~arm*density, nrow = 1) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = 'spearman') +
  labs(x = "Cliff's Delta in R2",
       y = "Cliff's Delta in R3") +
  coord_cartesian(xlim = c(-1,1),
                  ylim = c(-1,1)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

rr_corr <- ggpubr::ggarrange(rr12, rr13, rr23,
                  ncol = 1)
ggsave(sprintf('%sREPLICATEREPLICATE_CLIFF_CORR.png',out_path), rr_corr,
       height = 200, width = 350, units = 'mm',
       dpi = 300, limitsize = F)
