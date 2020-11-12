##### SDPG LIQUID GROWTH CURVES
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 11/09/2020

##### INITIALIZATION
library(ggplot2)
library(ggExtra)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(stringr)
library(zoo)

path.out <- 'figs/SDPG/'
path.out.gc <- 'figs/SDPG/GC/'
expt_name <- 'SDPG' 

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

# LOADING ALL THE DATA
# load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/growthcurver_results.RData')
# gc_results <- df_final; gc_stats <- growth_crv; gc_models <- models_all
# load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/growthcurver_raw_results.RData')
# gc_raw_results <- df_final; gc_raw_stats <- growth_crv; gc_raw_models <- models_all
# load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/growthrates_results.RData')
# gr_results <- growth_res; gr_data <- growth_dat; gr_exp <- growth_exp; gr_fit <- growth_fit
# load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/growthrates_raw_results.RData')
# gr_raw_results <- growth_res; gr_raw_data <- growth_dat; gr_raw_exp <- growth_exp; gr_raw_fit <- growth_fit
# load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/aleeza_results.RData')
# al_results <- res_plc; al_raw_results <- res_blk

# gc_results$arm <- factor(gc_results$arm, c('GLU','CAS','SDA'))
# gc_results$replicate <- factor(gc_results$replicate, c('R1','R2','R3'))
# gc_results$orf_name <- factor(gc_results$orf_name, sort(unique(gc_results$orf_name)))
# gc_results$bio_rep <- factor(gc_results$bio_rep, c(1,2))
# 
# gc_stats$arm <- factor(gc_stats$arm, c('GLU','CAS','SDA'))
# gc_stats$replicate <- factor(gc_stats$replicate, c('R1','R2','R3'))
# gc_stats$orf_name <- factor(gc_stats$orf_name, sort(unique(gc_stats$orf_name)))
# gc_stats$bio_rep <- factor(gc_stats$bio_rep, c(1,2))
# 
# gr_results$arm <- factor(gr_results$arm, c('GLU','CAS','SDA'))
# gr_results$replicate <- factor(gr_results$replicate, c('R1','R2','R3'))
# gr_results$orf_name <- factor(gr_results$orf_name, sort(unique(gr_results$orf_name)))
# gr_results$bio_rep <- factor(gr_results$bio_rep, c(1,2))
# 
# gr_data$arm <- factor(gr_data$arm, c('GLU','CAS','SDA'))
# gr_data$replicate <- factor(gr_data$replicate, c('R1','R2','R3'))
# gr_data$orf_name <- factor(gr_data$orf_name, sort(unique(gr_data$orf_name)))
# gr_data$bio_rep <- factor(gr_data$bio_rep, c(1,2))
# 
# gr_exp$arm <- factor(gr_exp$arm, c('GLU','CAS','SDA'))
# gr_exp$replicate <- factor(gr_exp$replicate, c('R1','R2','R3'))
# gr_exp$orf_name <- factor(gr_exp$orf_name, sort(unique(gr_exp$orf_name)))
# gr_exp$bio_rep <- factor(gr_exp$bio_rep, c(1,2))
# 
# gr_fit$arm <- factor(gr_fit$arm, c('GLU','CAS','SDA'))
# gr_fit$replicate <- factor(gr_fit$replicate, c('R1','R2','R3'))
# gr_fit$orf_name <- factor(gr_fit$orf_name, sort(unique(gr_fit$orf_name)))
# gr_fit$bio_rep <- factor(gr_fit$bio_rep, c(1,2))
# 
# al_results$Arm <- factor(al_results$Arm, c('GLU','CAS','SDA'))
# al_results$exp_rep <- factor(al_results$exp_rep, c('R1','R2','R3'))
# al_results$orf_name <- factor(al_results$orf_name, sort(unique(al_results$orf_name)))
# al_results$replicate <- factor(al_results$replicate, c(1,2))
# 
# gc_raw_results$arm <- factor(gc_raw_results$arm, c('GLU','CAS','SDA'))
# gc_raw_results$replicate <- factor(gc_raw_results$replicate, c('R1','R2','R3'))
# gc_raw_results$orf_name <- factor(gc_raw_results$orf_name, sort(unique(gc_raw_results$orf_name)))
# gc_raw_results$bio_rep <- factor(gc_raw_results$bio_rep, c(1,2))
# 
# gc_raw_stats$arm <- factor(gc_raw_stats$arm, c('GLU','CAS','SDA'))
# gc_raw_stats$replicate <- factor(gc_raw_stats$replicate, c('R1','R2','R3'))
# gc_raw_stats$orf_name <- factor(gc_raw_stats$orf_name, sort(unique(gc_raw_stats$orf_name)))
# gc_raw_stats$bio_rep <- factor(gc_raw_stats$bio_rep, c(1,2))
# 
# gr_raw_results$arm <- factor(gr_raw_results$arm, c('GLU','CAS','SDA'))
# gr_raw_results$replicate <- factor(gr_raw_results$replicate, c('R1','R2','R3'))
# gr_raw_results$orf_name <- factor(gr_raw_results$orf_name, sort(unique(gr_raw_results$orf_name)))
# gr_raw_results$bio_rep <- factor(gr_raw_results$bio_rep, c(1,2))
# 
# gr_raw_data$arm <- factor(gr_raw_data$arm, c('GLU','CAS','SDA'))
# gr_raw_data$replicate <- factor(gr_raw_data$replicate, c('R1','R2','R3'))
# gr_raw_data$orf_name <- factor(gr_raw_data$orf_name, sort(unique(gr_raw_data$orf_name)))
# gr_raw_data$bio_rep <- factor(gr_raw_data$bio_rep, c(1,2))
# 
# gr_raw_exp$arm <- factor(gr_raw_exp$arm, c('GLU','CAS','SDA'))
# gr_raw_exp$replicate <- factor(gr_raw_exp$replicate, c('R1','R2','R3'))
# gr_raw_exp$orf_name <- factor(gr_raw_exp$orf_name, sort(unique(gr_raw_exp$orf_name)))
# gr_raw_exp$bio_rep <- factor(gr_raw_exp$bio_rep, c(1,2))
# 
# gr_raw_fit$arm <- factor(gr_raw_fit$arm, c('GLU','CAS','SDA'))
# gr_raw_fit$replicate <- factor(gr_raw_fit$replicate, c('R1','R2','R3'))
# gr_raw_fit$orf_name <- factor(gr_raw_fit$orf_name, sort(unique(gr_raw_fit$orf_name)))
# gr_raw_fit$bio_rep <- factor(gr_raw_fit$bio_rep, c(1,2))
# 
# al_raw_results$Arm <- factor(al_raw_results$Arm, c('GLU','CAS','SDA'))
# al_raw_results$exp_rep <- factor(al_raw_results$exp_rep, c('R1','R2','R3'))
# al_raw_results$orf_name <- factor(al_raw_results$orf_name, sort(unique(al_raw_results$orf_name)))
# al_raw_results$replicate <- factor(al_raw_results$replicate, c(1,2))
# 
# gc_results <- gc_results[order(gc_results$arm, gc_results$replicate, gc_results$orf_name, gc_results$bio_rep),]
# gc_stats <- gc_stats[order(gc_stats$arm, gc_stats$replicate, gc_stats$orf_name, gc_stats$bio_rep),]
# gr_results <- gr_results[order(gr_results$arm, gr_results$replicate, gr_results$orf_name, gr_results$bio_rep),]
# gr_data <- gr_data[order(gr_data$arm, gr_data$replicate, gr_data$orf_name, gr_data$bio_rep),]
# gr_exp <- gr_exp[order(gr_exp$arm, gr_exp$replicate, gr_exp$orf_name, gr_exp$bio_rep),]
# gr_fit <- gr_fit[order(gr_fit$arm, gr_fit$replicate, gr_fit$orf_name, gr_fit$bio_rep),]
# al_results <- al_results[order(al_results$Arm, al_results$exp_rep, al_results$orf_name, al_results$replicate),]
# 
# gc_raw_results <- gc_raw_results[order(gc_raw_results$arm, gc_raw_results$replicate, gc_raw_results$orf_name, gc_raw_results$bio_rep),]
# gc_raw_stats <- gc_raw_stats[order(gc_raw_stats$arm, gc_raw_stats$replicate, gc_raw_stats$orf_name, gc_raw_stats$bio_rep),]
# gr_raw_results <- gr_raw_results[order(gr_raw_results$arm, gr_raw_results$replicate, gr_raw_results$orf_name, gr_raw_results$bio_rep),]
# gr_raw_data <- gr_raw_data[order(gr_raw_data$arm, gr_raw_data$replicate, gr_raw_data$orf_name, gr_raw_data$bio_rep),]
# gr_raw_exp <- gr_raw_exp[order(gr_raw_exp$arm, gr_raw_exp$replicate, gr_raw_exp$orf_name, gr_raw_exp$bio_rep),]
# gr_raw_fit <- gr_raw_fit[order(gr_raw_fit$arm, gr_raw_fit$replicate, gr_raw_fit$orf_name, gr_raw_fit$bio_rep),]
# al_raw_results <- al_raw_results[order(al_raw_results$Arm, al_raw_results$exp_rep, al_raw_results$orf_name, al_raw_results$replicate),]
# 
# # ALL RESULTS AND PHENOTYPE
# all_results <- data.frame(condition = al_results$Arm,
#            exp_rep = al_results$exp_rep,
#            orf_name = al_results$orf_name,
#            bio_rep = al_results$replicate,
#            aleeza = al_raw_results$dtime,
#            alleza_plc = al_results$dtime,
#            growthrates = gr_raw_results$dtime,
#            growthrates_plc = gr_results$dtime,
#            growthcurver = gc_raw_stats$t_gen,
#            growthcurver_plc = gc_stats$t_gen,
#            growthcurver_auc = gc_raw_stats$auc_l,
#            growthcurver_plc_auc = gc_stats$auc_l,
#            growthcurver_auc2 = gc_raw_stats$auc_e,
#            growthcurver_plc_auc2 = gc_stats$auc_e)
# 
# for (a in unique(all_results$condition)) {
#   for (r in unique(all_results$exp_rep[all_results$condition == a])) {
#     for (m in unique(names(all_results[all_results$condition == a & all_results$exp_rep == r,5:14]))) {
#       ref_dtime <- all_results[all_results$condition == a & all_results$exp_rep == r & all_results$orf_name == 'BF_control', m]
#       for (o in all_results$orf_name[all_results$condition == a & all_results$exp_rep == r & all_results$orf_name != 'BF_control']) {
#         orf_dtime <- all_results[all_results$condition == a & all_results$exp_rep == r & all_results$orf_name == o, m]
#         stat_res <- t.test(orf_dtime, ref_dtime)
#         all_results[all_results$condition == a & all_results$exp_rep == r & all_results$orf_name == o, sprintf('%s_p',m)] <-
#           stat_res$p.value
#         all_results[all_results$condition == a & all_results$exp_rep == r & all_results$orf_name == o, sprintf('%s_tstat',m)] <-
#           stat_res$statistic[[1]]
#         if (stat_res$p.value <= 0.05) {
#           if (stat_res$statistic[[1]] < 0) {
#             if (!str_detect(m, 'auc')) {
#               phenotype <- 'Beneficial'
#             } else {
#               phenotype <- 'Deleterious'
#             }
#           } else {
#             if (!str_detect(m, 'auc')) {
#               phenotype <- 'Deleterious'
#             } else {
#               phenotype <- 'Beneficial'
#             }
#           }
#         } else {
#           phenotype <- 'Neutral'
#         }
#         all_results[all_results$condition == a & all_results$exp_rep == r & all_results$orf_name == o, sprintf('%s_phenotype',m)] <- phenotype
#       }
#     }
#   }
# }
# 
# melt1 <-melt(all_results[,c(1:14)],
#              id.vars = c('condition','exp_rep','orf_name','bio_rep'),
#              variable.name = 'method')
# melt1$method <- as.character(melt1$method)
# melt2 <- melt(all_results[,c(1:4, c(1:ncol(all_results))[str_detect(names(all_results), 'phenotype')])],
#               id.vars = c('condition','exp_rep','orf_name','bio_rep'),
#               variable.name = 'method',
#               value.name = 'phenotype')
# melt2$method <- as.character(melt2$method)
# melt2$method <- str_remove_all(melt2$method, '_phenotype')
# 
# all_results2 <- merge(melt1, melt2,
#                       by = c('condition','exp_rep','orf_name','bio_rep','method'))
# 
# # SAVING EVERYTHING TOGETHER
# save(gc_results, gc_stats, gc_models,
#      gc_raw_results, gc_raw_stats, gc_raw_models,
#      gr_results, gr_data, gr_exp, gr_fit,
#      gr_raw_results, gr_raw_data, gr_raw_exp, gr_raw_fit,
#      al_results, al_raw_results,
#      all_results, all_results2,
#      file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/all_methods.RData')

# LOADING THE COMBINED RESULTS
load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/all_methods.RData')
head(all_results)

##### GROTHCURVER RESULTS
liq_dtime <- ggplot(all_results,
                    aes(x = orf_name, y = growthcurver_plc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(col = exp_rep), size = 1) +
  scale_color_discrete(name = 'Replicate') +
  # geom_text(aes(label = bio_rep), size = 3, color = 'white') +
  # scale_color_manual(name = 'Phenotype',
  #                      breaks = c('Beneficial','Neutral','Deleterious'),
  #                      values = c('Deleterious'='#3F51B5',
  #                                 'Neutral'='#212121',
  #                                 'Beneficial'='#FFC107'),
  #                      na.value = 'grey50') +
  labs(y = 'Doubling Time (mins)',
       x = 'ORFs') +
  stat_compare_means(method = 't.test', ref.group = "BF_control",
                     label = "p.signif", label.y = 350,
                     size = 1.5, angle = 0, vjust = 0.5) +
  # stat_compare_means(aes(group = orf_name),
  #                    label = "p.signif", label.y = 320,
  #                    size = 1.5, hide.ns = T) +
  coord_flip() +
  facet_wrap(~condition) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sLIQUID_DTIME2.png',path.out), liq_dtime,
       height = 200, width = 400, units = 'mm',
       dpi = 300, limitsize = F)

liq_auc <- ggplot(gc_stats,
       aes(x = orf_name, y = auc_e)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(col = replicate)) +
  stat_compare_means(method = 't.test', ref.group = "BF_control",
                     label = "p.signif",
                     size = 1.5, angle = 0, vjust = 0.5) +
  scale_color_discrete(name = 'Replicate') +
  labs(y = 'Area Under the Curve',
       x = 'ORFs') +
  coord_flip() +
  facet_wrap(~arm) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sLIQUID_AUC.png',path.out), liq_auc,
       height = 200, width = 400, units = 'mm',
       dpi = 300, limitsize = F)


liq_gc <- ggplot(gc_results,
       aes(x = time)) +
  # geom_point(aes(y = od, col = replicate), alpha = 0.2) +
  geom_line(aes(y = pred.od, col = bio_rep)) +
  scale_color_discrete(name = 'BioRep') +
  geom_text(data = gc_stats, aes(x = 2000, y = 1, label = round(t_gen,2), col = bio_rep),
            position = position_dodge(width = 2000), size = 1.5) +
  geom_text(data = gc_stats, aes(x = 2000, y = 0.5, label = round(auc_l,2), col = bio_rep),
            position = position_dodge(width = 2000), size = 1.5) +
  labs(y = 'OD',
       x = 'Time') +
  facet_wrap(.~arm*replicate*orf_name, ncol = 14) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf('%sLIQUID_GC.png',path.out), liq_gc,
       height = 500, width = 500, units = 'mm',
       dpi = 300, limitsize = F)

##### PLOTTING COMBINED RESULTS
head(all_results2)
all_results2$phenotype <- factor(all_results2$phenotype, c('Deleterious','Neutral','Beneficial'))

# a <- 'GLU'
# r <- 'R1'
# o <- 'YBR196C-A'

for (a in unique(all_results2$condition)) {
  for (r in unique(all_results2$exp_rep[all_results2$condition == a])) {
    for (o in unique(all_results2$orf_name[all_results2$condition == a & all_results2$exp_rep == r & all_results2$orf_name != 'BF_control'])) {
      ref_gc <- gc_results[gc_results$arm == a &
                             gc_results$replicate == r &
                             gc_results$orf_name == 'BF_control',]
      orf_gc <- gc_results[gc_results$arm == a &
                             gc_results$replicate == r &
                             gc_results$orf_name == o,]
      ref_res <- all_results2[!str_detect(all_results2$method,'auc') &
                                all_results2$condition == a &
                                all_results2$exp_rep == r &
                                all_results2$orf_name == 'BF_control',]
      orf_res <- all_results2[!str_detect(all_results2$method,'auc') &
                                all_results2$condition == a &
                                all_results2$exp_rep == r &
                                all_results2$orf_name == o,]
      ref_auc <- all_results2[str_detect(all_results2$method,'auc') &
                                all_results2$condition == a &
                                all_results2$exp_rep == r &
                                all_results2$orf_name == 'BF_control',]
      orf_auc <- all_results2[str_detect(all_results2$method,'auc') &
                                all_results2$condition == a &
                                all_results2$exp_rep == r &
                                all_results2$orf_name == o,]
      
      plot_dtime <- ggplot() +
        geom_point(data = ref_res,
                   aes(x = method, y = value, shape = bio_rep, fill = phenotype),
                   size = 3, alpha = 0.5) +
        geom_text_repel(data = ref_res, aes(x = method, y = value, label = sprintf('%0.2f',value)),
                        size = 2) +
        geom_point(data = orf_res,
                   aes(x = method, y = value, shape = bio_rep, fill = phenotype),
                   size = 3) +
        geom_text_repel(data = orf_res, aes(x = method, y = value, label = sprintf('%0.2f',value)),
                        size = 2) +
        # scale_y_continuous(minor_breaks = seq(0,400,20)) +
        scale_x_discrete(limits = c("aleeza", "alleza_plc", "growthrates", "growthrates_plc", "growthcurver", "growthcurver_plc")) +
        scale_fill_manual(name = 'Phenotype',
                          breaks = c('Beneficial','Neutral','Deleterious'),
                          values = c('Deleterious'='#3F51B5',
                                     'Neutral'='#212121',
                                     'Beneficial'='#FFC107'),
                          na.value = 'grey50', drop = F) +
        scale_shape_manual(name = 'bioRep',
                           breaks = c(1,2),
                           values = c(21,24)) +
        labs(x = 'Methods',
             y = 'Doubling Time (mins)') +
        # coord_cartesian(ylim = c(20,360)) +
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
                   aes(x = method, y = value, shape = bio_rep, fill = phenotype),
                   size = 3, alpha = 0.5) +
        geom_text_repel(data = ref_auc,
                        aes(x = method, y = value, label = sprintf('%0.2f',value)),
                        size = 2) +
        geom_point(data = orf_auc,
                   aes(x = method, y = value, shape = bio_rep, fill = phenotype),
                   size = 3) +
        geom_text_repel(data = orf_auc,
                        aes(x = method, y = value, label = sprintf('%0.2f',value)),
                        size = 2) +
        scale_x_discrete(limits = c("growthcurver_auc", "growthcurver_plc_auc", "growthcurver_auc2", "growthcurver_plc_auc2")) +
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
        labs(x = 'Methods',
             y = 'AUC') +
        # coord_cartesian(ylim = c(10,17000)) +
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
                                                                               widths = c(6,4),
                                                                               common.legend = T,
                                                                               legend = 'right'),
                                                    ncol = 1, align = 'h',
                                                    heights = c(1,2)),
                                  top = text_grob(sprintf('%s | %s | %s',a,r,o), face = "bold", size = 10))
      ggsave(sprintf('%s%s_%s_%s.png',path.out.gc,a,r,o), plot_all,
             height = 200, width = 350, units = 'mm',
             dpi = 300, limitsize = F)
    }
  }
}

##### ROLLING DIFFERENCE
df$Diff <- df$Count - dplyr::lag(df$Count, n = 5)
ref_gc$diff <- ref_gc$od - dplyr::lag(ref_gc$od, n = 1)

plot(ref_gc$od - dplyr::lag(ref_gc$od, n = 1))
