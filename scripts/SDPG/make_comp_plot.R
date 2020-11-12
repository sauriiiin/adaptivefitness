
make_comp_plot <- function(pair_comp, mat_type,
                           g1_orf_name, g1_density, g1_arm, g1_rep,
                           g2_orf_name, g2_density, g2_arm, g2_rep) {
  
  ##### TEXT SIZE
  titles <- 11
  txt <- 11
  lbls <- 13
  
  g <- ggplot(pair_comp[pair_comp$g1_orf_name %in% g1_orf_name &
                          pair_comp$g1_density %in% g1_density &
                          pair_comp$g1_arm %in% g1_arm &
                          pair_comp$g1_rep %in% g1_rep &
                          pair_comp$g2_orf_name %in% g2_orf_name &
                          pair_comp$g2_density %in% g2_density &
                          pair_comp$g2_arm %in% g2_arm &
                          pair_comp$g2_rep %in% g2_rep,])
  
  if (mat_type == 'magnitude') {
    temp_comp_mat <- g + geom_tile(aes(x = g1, y = g2, fill = magnitude), col = 'black') +
      scale_fill_brewer(name = 'Magnitude',
                        breaks = c('negligible','small','medium','large'),
                        palette = 'Set1', drop = F)
  } else if (mat_type == 'pvalue') {
    temp_comp_mat <- g + geom_tile(aes(x = g1, y = g2, fill = pvalue), col = 'black') +
      geom_text(aes(x = g1, y = g2, label = round(pvalue,3)), col = 'black', size = 5) +
      scale_fill_distiller(name = 'P-Value',
                          palette = 'RdYlGn',
                          limits = c(0,1),
                          breaks = c(0,0.5,1),
                          direction = 1)
  } else if (mat_type == 'cliffdelta') {
    temp_comp_mat <- g + geom_tile(aes(x = g1, y = g2, fill = cliffdelta), col = 'black') +
      geom_text(aes(x = g1, y = g2, label = round(cliffdelta,3)), col = 'black', size = 5) +
      scale_fill_distiller(name = "Cliff's Delta",
                           palette = 'RdYlGn',
                           limits = c(-1,1),
                           breaks = c(-1,0,1),
                           direction = 1)
  } else {
    temp_comp_mat <- g + geom_tile(aes(x = g1, y = g2, fill = phenotype), col = 'black') +
      # geom_text(aes(x = g1, y = g2), col = 'black', size = 5) +
      scale_fill_manual(name = 'Phenotype',
                         breaks = c('Deleterious','Neutral','Beneficial','NA'),
                         values = c('Deleterious'='#3F51B5',
                                    'Neutral'='#212121',
                                    'Beneficial'='#FFC107',
                                    'NA'='transparent'),
                         drop = F) 
  }
  
  comp_mat <- temp_comp_mat +
    theme_linedraw() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.title = element_text(size = titles),
          legend.text = element_text(size = txt),
          legend.key.width = unit(10, "mm"),
          legend.key.height = unit(2, "mm"),
          # legend.box.spacing = unit(1,"mm"),
          legend.position = 'top')
  
  comp_x <- g +
    geom_point(aes(x = g1, y = 'ORF', col = g1_orf_name), size = 3) +
    geom_tile(aes(x = g1, y = 'Replicate', fill = g1_rep), col = 'black') +
    geom_tile(aes(x = g1, y = 'Condition', fill = g1_arm), col = 'black') +
    geom_tile(aes(x = g1, y = 'Density', fill = g1_density), col = 'black') +
    geom_text(aes(x = g1, y = 'Replicate', label = g1_rep), col = 'black', size = 3) +
    geom_text(aes(x = g1, y = 'Condition', label = g1_arm), col = 'black', size = 3) +
    geom_text(aes(x = g1, y = 'Density', label = g1_density), col = 'black', size = 3) +
    labs(x = 'Group 1',
         y = '') +
    scale_y_discrete(breaks = c('Density','Condition','Replicate','ORF'),
                     limits = c('Density','Condition','Replicate','ORF')) +
    scale_fill_manual(name = 'Axis',
                      breaks = c("1536","6144",
                                 "GLU","CAS","SDA",
                                 "R1","R2","R3"),
                      values = c("1536" = "#607D8B", "6144" = "#FFA000",
                                 "GLU" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F",
                                 "R1" = "#E1BEE7", "R2" = "#E040FB", "R3" = "#7B1FA2"),
                      drop = F) +
    scale_color_discrete(name = 'ORFs',
                       breaks = sort(unique(pair_comp$g1_orf_name[pair_comp$g1_orf_name != 'REF'])),
                       drop = F) +
    theme_linedraw() +
    theme(axis.text = element_blank(),
          axis.title = element_text(size = titles),
          # axis.title.y = element_blank(),
          legend.title = element_text(size = titles),
          legend.text = element_text(size = txt),
          legend.key.size = unit(3, "mm"),
          legend.background = element_blank(),
          # legend.box.spacing = unit(1,"mm"),
          legend.position = 'bottom')
  
  comp_y <- g +
    geom_point(aes(y = g2, x = 'ORF', col = g2_orf_name), size = 3) +
    geom_tile(aes(y = g2, x = 'Replicate', fill = g2_rep), col = 'black') +
    geom_tile(aes(y = g2, x = 'Condition', fill = g2_arm), col = 'black') +
    geom_tile(aes(y = g2, x = 'Density', fill = g2_density), col = 'black') +
    geom_text(aes(y = g2, x = 'Replicate', label = g2_rep), col = 'black', size = 3, angle = 90) +
    geom_text(aes(y = g2, x = 'Condition', label = g2_arm), col = 'black', size = 3, angle = 90) +
    geom_text(aes(y = g2, x = 'Density', label = g2_density), col = 'black', size = 3, angle = 90) +
    labs(y = 'Group 2') +
    scale_x_discrete(breaks = c('Density','Condition','Replicate','ORF'),
                     limits = c('Density','Condition','Replicate','ORF')) +
    scale_fill_manual(name = 'Axis',
                      breaks = c("1536","6144",
                                 "GLU","CAS","SDA",
                                 "R1","R2","R3"),
                      values = c("1536" = "#607D8B", "6144" = "#FFA000",
                                 "GLU" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F",
                                 "R1" = "#E1BEE7", "R2" = "#E040FB", "R3" = "#7B1FA2"),
                      guide = F) +
    scale_color_discrete(name = 'ORFs',
                       breaks = sort(unique(pair_comp$g2_orf_name[pair_comp$g2_orf_name != 'REF'])),
                       drop = F,
                       guide = F) +
    theme_linedraw() +
    theme(axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = titles),
          legend.title = element_text(size = titles),
          legend.text = element_text(size = txt),
          legend.key.size = unit(3, "mm")
          # legend.box.spacing = unit(1,"mm")
          )
  
  comp_top <- ggpubr::ggarrange(comp_y, comp_mat,
                                widths = c(1,6),
                                align = 'hv',
                                legend = 'none')
  comp_bot <- ggpubr::ggarrange(NULL, comp_x,
                                widths = c(1,6),
                                align = 'v',
                                legend = 'none')
  comp_leg1 <- get_legend(comp_mat)
  comp_leg2 <- get_legend(comp_x)
  
  comp_main <- ggpubr::ggarrange(comp_top, comp_bot, ncol = 1,
                                 heights = c(5,1))
  comp_leg <- ggpubr::ggarrange(as_ggplot(comp_leg2), as_ggplot(comp_leg1), ncol = 1)
  comp_plot <- ggpubr::ggarrange(comp_main, comp_leg, ncol = 1,
                                 heights = c(6,1))
  return(comp_plot)
}


