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
out_path_den <- 'figs/SDPG/DEN/'

##### FIGURE TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### LOAD DATA
load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/solid_results.RData')
load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/liquid_results.RData')

z.test2sam = function(a, b){
  var.a = var(a)
  var.b = var(b)
  n.a = length(a)
  n.b = length(b)
  zeta = (mean(a) - mean(b)) / (sqrt(var.a/n.a + var.b/n.b))
  return(zeta)
}

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


##### COMBINED PHENOTYPE WITH LIQ AUC ONLY
pheno_mat2 <- ggplot(melt(all_phenotype[c(1,2,3,11,14:16)], id.vars = c('orf_name','arm','replicate'),
                         variable.name = 'method', value.name = 'phenotype'),
                    aes(x = method, y = orf_name)) +
  geom_tile(aes(fill = phenotype), col = 'black') +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Deleterious','Neutral','Beneficial'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#212121',
                               'Beneficial'='#FFC107'),
                    na.value = 'grey50') +
  scale_x_discrete(limits = names(all_phenotype[c(11,14:16)]),
                   # labels = str_remove_all(names(all_phenotype[c(11,14:16)]), '_phenotype')) +
                   labels = c('liquid','solid_1536','solid_6144','expectation')) +
  labs(x = 'Method',
       y = 'ORFs') +
  facet_wrap(.~arm*replicate) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sALL_PHENOTYPE_MATRIX2.png',out_path), pheno_mat2,
       height = 270, width = 300, units = 'mm',
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
  scale_x_discrete(labels = c('liquid','solid_1536','solid_6144','expectation')) +
  facet_wrap(.~arm, nrow = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sCOMMON_PHENOTYPE_MATRIX.png',out_path), common_pheno_mat,
       height = 100, width = 300, units = 'mm',
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

##### ALL PAIRWISE COMPARISONS
cor_dat <- NULL
cols <- data.frame()
for (d in unique(solid_stats$density)) {
  for (a in unique(solid_stats$arm[solid_stats$density == d])) {
    for (r in unique(solid_stats$replicate[solid_stats$density == d & 
                                solid_stats$arm == a])) {
      cor_dat <- cbind(cor_dat, t(t(solid_stats$cliff.delta[solid_stats$density == d &
                                                          solid_stats$arm == a &
                                                          solid_stats$replicate == r &
                                                          solid_stats$saturation == 'Saturated'])))
      cols <- rbind(cols,
                    data.frame(title = sprintf('%d_%s_%s',d,a,r),
                               density = d,
                               arm = a,
                               replicate = r))
    }
  }
}
colnames(cor_dat) <- cols$title
head(cor_dat)

ha <- HeatmapAnnotation(
  Density = cols$density,
  Condition = cols$arm,
  Replicate = cols$replicate,
  col = list(Density = c("1536" = "#607D8B", "6144" = "#FFA000"),
             Condition = c("GLU" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F"),
             Replicate = c("R1" = "#E1BEE7", "R2" = "#E040FB", "R3" = "#7B1FA2")),
  gp = gpar(col = "black", lwd = 0.2)
)
ra <- rowAnnotation(
  Density = cols$density,
  Condition = cols$arm,
  Replicate = cols$replicate,
  Saturation = cols$saturation,
  col = list(Density = c("1536" = "#607D8B", "6144" = "#FFA000"),
             Condition = c("GLU" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F"),
             Replicate = c("R1" = "#E1BEE7", "R2" = "#E040FB", "R3" = "#7B1FA2")),
  gp = gpar(col = "black", lwd = 0.2),
  show_annotation_name = FALSE,
  show_legend = FALSE
)
mat <- cor(cor_dat, use = 'complete.obs', method = 'spearman')

ht <- Heatmap(mat,
              name = "Correlation",
              column_title = "The Famous 28 Correlation",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%0.2f", mat[i, j]), x, y, gp = gpar(fontsize = 5))
              },
              column_title_gp = gpar(fontsize = 12, fontface = "bold"),
              show_row_dend = F,
              rect_gp = gpar(col = "black", lwd = 0.2),
              show_row_names = F,
              show_column_names = F,
              border = T,
              top_annotation = ha)
png(file = sprintf('%sALL_CORR.png',out_path),
    width = 180, height = 160, units = "mm", res = 300,
    bg = "transparent")
draw(ra + ht)
dev.off()


#### PHENOTYPE DYNAMICS
pheno_dynamics <- ggplot(solid_stats,
                     aes(x = hours, y = cliff.delta)) +
  geom_hline(yintercept = c(-0.474,0.474), col = '#212121', linetype = 'dashed', lwd = 0.3) +
  geom_hline(yintercept = c(-0.33,0.33), col = '#757575', linetype = 'dashed', lwd = 0.3) +
  geom_hline(yintercept = c(-0.147,0.147), col = '#BDBDBD', linetype = 'dashed', lwd = 0.3) +
  # stat_regline_equation(aes(label =  ..eq.label..)) +
  # stat_cor(method = "pearson") +
  geom_line(aes(group = 1, col = replicate)) +
  geom_point(aes(col = pitt_phenotype)) +
  facet_wrap(.~orf_name*arm*density*replicate,
             ncol = 15,
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
        # axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sPHENOTYPE_DYNAMICS.jpg",out_path),pheno_dynamics,
       height = 1500, width = 1000, units = 'mm', limitsize = F,
       dpi = 300)


##### LIQUID GROWTH CURVES

for (a in unique(all_results2$condition)) {
  for (r in unique(all_results2$exp_rep[all_results2$condition == a])) {
    for (o in unique(all_results2$orf_name[all_results2$condition == a & all_results2$exp_rep == r & all_results2$orf_name != 'BF_control'])) {
      ref_gc <- gc_results[gc_results$arm == a &
                             gc_results$replicate == r &
                             gc_results$orf_name == 'BF_control',]
      orf_gc <- gc_results[gc_results$arm == a &
                             gc_results$replicate == r &
                             gc_results$orf_name == o,]
      ref_res <- all_results2[all_results2$method == 'growthcurver_plc' &
                                all_results2$condition == a &
                                all_results2$exp_rep == r &
                                all_results2$orf_name == 'BF_control',]
      orf_res <- all_results2[all_results2$method == 'growthcurver_plc' &
                                all_results2$condition == a &
                                all_results2$exp_rep == r &
                                all_results2$orf_name == o,]
      ref_auc <- all_results2[all_results2$method == 'growthcurver_plc_auc' &
                                all_results2$condition == a &
                                all_results2$exp_rep == r &
                                all_results2$orf_name == 'BF_control',]
      orf_auc <- all_results2[all_results2$method == 'growthcurver_plc_auc' &
                                all_results2$condition == a &
                                all_results2$exp_rep == r &
                                all_results2$orf_name == o,]
      
      plot_dtime <- ggplot() +
        geom_point(data = ref_res,
                   aes(x = orf_name, y = value, shape = bio_rep, fill = phenotype),
                   size = 3, alpha = 0.5) +
        geom_text_repel(data = ref_res, aes(x = orf_name, y = value, label = sprintf('%0.2f',value)),
                        size = 2) +
        geom_point(data = orf_res,
                   aes(x = orf_name, y = value, shape = bio_rep, fill = phenotype),
                   size = 3) +
        geom_text_repel(data = orf_res, aes(x = orf_name, y = value, label = sprintf('%0.2f',value)),
                        size = 2) +
        # scale_y_continuous(minor_breaks = seq(0,400,20)) +
        # scale_x_discrete(limits = c("aleeza", "alleza_plc", "growthrates", "growthrates_plc", "growthcurver", "growthcurver_plc")) +
        scale_fill_manual(name = 'Phenotype',
                          breaks = c('Beneficial','Neutral','Deleterious'),
                          values = c('Deleterious'='#3F51B5',
                                     'Neutral'='#212121',
                                     'Beneficial'='#FFC107'),
                          na.value = 'grey50', drop = F) +
        scale_shape_manual(name = 'bioRep',
                           breaks = c(1,2),
                           values = c(21,24)) +
        labs(x = 'ORF',
             y = 'Doubling Time (mins)') +
        coord_cartesian(ylim = c(20,360)) +
        theme_linedraw() +
        theme(plot.title = element_text(size = titles + 2,
                                        face = 'bold',
                                        hjust = 0.5),
              axis.title = element_text(size = titles),
              axis.text = element_text(size = txt),
              legend.title = element_text(size = titles),
              legend.text = element_text(size = txt),
              legend.key.size = unit(3, "mm"),
              legend.position = "bottom",
              legend.box.spacing = unit(0.5,"mm"),
              strip.text = element_text(size = txt,
                                        margin = margin(0.1,0,0.1,0, "mm"))) +
        guides(fill = guide_legend(override.aes = list(shape = 22)))
      
      plot_auc <- ggplot() +
        geom_point(data = ref_auc,
                   aes(x = orf_name, y = value, shape = bio_rep, fill = phenotype),
                   size = 3, alpha = 0.5) +
        geom_text_repel(data = ref_auc,
                        aes(x = orf_name, y = value, label = sprintf('%0.2f',value)),
                        size = 2) +
        geom_point(data = orf_auc,
                   aes(x = orf_name, y = value, shape = bio_rep, fill = phenotype),
                   size = 3) +
        geom_text_repel(data = orf_auc,
                        aes(x = orf_name, y = value, label = sprintf('%0.2f',value)),
                        size = 2) +
        # scale_x_discrete(limits = c("growthcurver_auc", "growthcurver_plc_auc", "growthcurver_auc2", "growthcurver_plc_auc2")) +
        # scale_y_continuous(minor_breaks = seq(0,20000,1000)) +
        scale_shape_manual(name = 'bioRep',
                           breaks = c(1,2),
                           values = c(21,24)) +
        scale_fill_manual(name = 'Phenotype',
                          breaks = c('Beneficial','Neutral','Deleterious'),
                          values = c('Deleterious'='#3F51B5',
                                     'Neutral'='#212121',
                                     'Beneficial'='#FFC107'),
                          na.value = 'grey50', drop = F) +
        labs(x = 'ORF',
             y = 'AUC') +
        coord_cartesian(ylim = c(10,17000)) +
        theme_linedraw() +
        theme(plot.title = element_text(size = titles + 2,
                                        face = 'bold',
                                        hjust = 0.5),
              axis.title = element_text(size = titles),
              axis.text = element_text(size = txt),
              legend.title = element_text(size = titles),
              legend.text = element_text(size = txt),
              legend.key.size = unit(3, "mm"),
              legend.position = "bottom",
              legend.box.spacing = unit(0.5,"mm"),
              strip.text = element_text(size = txt,
                                        margin = margin(0.1,0,0.1,0, "mm"))) +
        guides(fill = guide_legend(override.aes = list(shape = 22)))
      
      plot_gc <- ggplot() +
        geom_line(data = ref_gc, aes(x = time, y = od, linetype = bio_rep),
                  alpha = 0.5, lwd = 1) +
        geom_line(data = orf_gc, aes(x = time, y = od, linetype = bio_rep),
                  lwd = 1) +
        # geom_point(data = ref_gc, aes(x = time, y = od, shape = bio_rep),
        #           alpha = 0.5, size = 1) +
        # geom_point(data = orf_gc, aes(x = time, y = od, shape = bio_rep),
        #           size = 1) +
        geom_text(aes(x = 0, y = 5.5, label = 'BF_control'),
                  alpha = 0.5, size = 2.5, hjust = 0) +
        geom_text(aes(x = 0, y = 5.1, label = 'Mutant'),
                  size = 2.5, hjust = 0) +
        # scale_y_log10() +
        labs(x = 'Time (mins)',
             y = 'PLC OD') +
        scale_linetype_discrete(name = 'bioRep') +
        # scale_shape_manual(name = 'bioRep',
        #                    breaks = c(1,2),
        #                    values = c(21,24)) +
        coord_cartesian(ylim = c(0.01,6)) +
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
      
      plot_all <- annotate_figure(ggpubr::ggarrange(plot_gc, ggpubr::ggarrange(plot_dtime, plot_auc,
                                                                               nrow = 1, align = 'hv',
                                                                               widths = c(1,1),
                                                                               common.legend = T,
                                                                               legend = 'right'),
                                                    ncol = 2,
                                                    widths = c(1,1)),
                                  top = text_grob(sprintf('%s | %s | %s',a,r,o), face = "bold", size = 10))
      ggsave(sprintf('%s%s_%s_%s.png',out_path_gc,a,r,o), plot_all,
             height = 100, width = 300, units = 'mm',
             dpi = 300, limitsize = F)
    }
  }
}

##### FITNESS DISTRIBUTION
for (o in unique(solid_fit$orf_name[!(solid_fit$orf_name %in% c('YHR201W-A','BOR','BF_control'))])) {
  for (a in unique(solid_fit$arm[solid_fit$orf_name == o])) {
    for (d in unique(solid_fit$density[solid_fit$orf_name == o & solid_fit$arm == a])) {
      for (r in unique(solid_fit$replicate[solid_fit$orf_name == o & solid_fit$arm == a & solid_fit$density == d & solid_fit$saturation == 'Saturated'])) {
        den_plot <- ggplot(merge(solid_fit[solid_fit$saturation == 'Saturated' &
                                 solid_fit$orf_name %in% c(o, 'BF_control') &
                                 solid_fit$density == d &
                                 solid_fit$arm == a &
                                 solid_fit$replicate == r &
                                 !is.na(solid_fit$orf_name),], solid_phenotype, by = c('arm','replicate','orf_name'), all.x = T),
               aes(x = fitness, y = orf_name, group = orf_name, fill = pitt_phenotype_6144)) +
          geom_density_ridges(quantile_lines = TRUE,
                              scale = 3, alpha = 0.9, size = 0.3,
                              vline_size = 0.2, vline_color = "black",
                              na.rm = T) +
          labs(x = "Fitness",
               y = 'ORF') +
          scale_fill_manual(name = 'Phenotype',
                            breaks = c('Beneficial','Neutral','Deleterious'),
                            values = c('Deleterious'='#3F51B5',
                                       'Neutral'='#212121',
                                       'Beneficial'='#FFC107'),
                            na.value = 'grey50', drop = F) +
          facet_wrap(.~arm*density*replicate,
                     nrow = 3) +
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
                                          margin = margin(0.1,0,0.1,0, "mm"))) +
          coord_cartesian(xlim = c(0.8, 1.2))
        ggsave(sprintf("%s%s_%s_%d_%s.jpg",out_path_den,o,a,d,r), den_plot,
               height = 100, width = 110, units = 'mm',
               dpi = 300)
      }
    }
  }
}


