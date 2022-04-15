hi <- merge(data.fit[data.fit$attempt == 'init' & data.fit$condition == 'SA' & data.fit$hours == 39, 
                     c('condition','strain_id','orf_name','category','fitness','fitness_','plate_no','plate_col','plate_row')],
            data.fit[data.fit$attempt == 'lin' & data.fit$condition == 'SA' & data.fit$hours == 42, 
                     c('condition','strain_id','orf_name','fitness','fitness_','plate_no','plate_col','plate_row')],
            by = c('condition','strain_id','orf_name','plate_no','plate_col','plate_row'), suffixes = c('_init','_lin'))
head(hi)

hi %>%
  ggplot(aes(x = fitness_init, y = fitness_lin)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  facet_wrap(.~category)

hi %>%
  group_by(category, strain_id) %>%
  summarize(fitness_init = median(fitness_init, na.rm = T),
            fitness_lin = median(fitness_lin, na.rm = T)) %>%
  ggplot(aes(x = fitness_init, y = fitness_lin)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  facet_wrap(.~category)


hi2 <- merge(data.fit.sum[data.fit.sum$attempt == 'init' & data.fit.sum$condition == 'GA' & data.fit.sum$hours == 21 &
                            data.fit.sum$category == 'Transient',
                          c('condition','strain_id','orf_name','fit.median')],
             data.fit.sum[data.fit.sum$attempt == 'lin' & data.fit.sum$condition == 'GA' & data.fit.sum$hours == 24 &
                            data.fit.sum$category == 'Transient',
                          c('condition','strain_id','orf_name','fit.median')],
             by = c('condition','strain_id','orf_name'), suffixes = c('_init','_lin'))
head(hi2)

hi2 %>%
  ggplot(aes(x = fit.median_init, y = fit.median_lin)) +
  geom_abline() +
  geom_point() +
  stat_cor(method = 'spearman') +
  coord_cartesian(xlim = c(0,1.5),
                  ylim = c(0,1.5))


hi3 <- merge(data.fit.sum[data.fit.sum$attempt == 'init' & data.fit.sum$condition == 'GA' & data.fit.sum$hours == 21 &
                            data.fit.sum$category == 'Transient',
                          c('condition','strain_id','orf_name','fitness.median')],
             data.fit.sum[data.fit.sum$attempt == 'valid' & data.fit.sum$condition == 'GA' & data.fit.sum$hours == 26 &
                            data.fit.sum$category == 'Transient',
                          c('condition','strain_id','orf_name','fitness.median')],
             by = c('condition','strain_id','orf_name'), suffixes = c('_init','_valid'))
head(hi3)

hi3 %>%
  ggplot(aes(x = fit.median_init, y = fit.median_valid)) +
  geom_abline() +
  geom_point() +
  stat_cor(method = 'spearman') +
  coord_cartesian(xlim = c(0,1.5),
                  ylim = c(0,1.5))


###### CORR HEATMAP PLOT
#data.fit.sum3 from oe_analysis
mat <- cor(data.fit.sum3[,c(4:45)], use = 'pairwise.complete.obs', method = 'spearman')
mat 

col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

ha <- HeatmapAnnotation(
  Data = str_split(col.names[4:45], '_', simplify = T)[,1],
  Attempt = str_split(col.names[4:45], '_', simplify = T)[,2],
  Condition  = str_split(col.names[4:45], '_', simplify = T)[,3],
  col = list(Data = c("cs" = "#E1BEE7", "auc" = "#E040FB", "growth" = "#7B1FA2"),
             Attempt = c("init" = "Navy", "valid" = "#607D8B", "lin" = "#FFA000"),
             Condition = c('GA' = '#312921',
                           'GA2' = '#312921',
                           'SA' = '#ccaf9b',
                           'HU' = '#a06b39',
                           'HO' = '#873c1e',
                           'CF' = '#c6c3b3',
                           'DM' = '#a63d4c',
                           'FL' = '#c8777b',
                           'TN' = '#dcaab1')),
  gp = gpar(col = "black", lwd = 0.2)
)

ra <- rowAnnotation(
  Data = str_split(col.names[4:45], '_', simplify = T)[,1],
  Attempt = str_split(col.names[4:45], '_', simplify = T)[,2],
  Condition  = str_split(col.names[4:45], '_', simplify = T)[,3],
  col = list(Data = c("cs" = "#E1BEE7", "auc" = "#E040FB", "growth" = "#7B1FA2"),
             Attempt = c("init" = "Navy", "valid" = "#607D8B", "lin" = "#FFA000"),
             Condition = c('GA' = '#312921',
                           'GA2' = '#312921',
                           'SA' = '#ccaf9b',
                           'HU' = '#a06b39',
                           'HO' = '#873c1e',
                           'CF' = '#c6c3b3',
                           'DM' = '#a63d4c',
                           'FL' = '#c8777b',
                           'TN' = '#dcaab1')),
  gp = gpar(col = "black", lwd = 0.2),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

hm <- Heatmap(mat,
              name = "Correlation",
              column_title = "Fitness Estimations",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%0.2f", mat[i, j]), x, y, gp = gpar(fontsize = 4))
              },
              jitter = T,
              column_title_gp = gpar(fontsize = 12, fontface = "bold"),
              show_row_dend = F,
              rect_gp = gpar(col = "black", lwd = 0.2),
              show_row_names = F,
              show_column_names = F,
              border = T,
              top_annotation = ha)
draw(hm+ra)