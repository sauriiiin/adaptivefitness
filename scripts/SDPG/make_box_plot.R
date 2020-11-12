
make_box_plot <- function(lidres,
                          g1_orf_name, g1_density, g1_arm, g1_rep,
                          g2_orf_name, g2_density, g2_arm, g2_rep) {
  
  ##### TEXT SIZE
  titles <- 11
  txt <- 11
  lbls <- 13
  
  data <- lidres[lidres$time == 'Saturated' &
                   !is.na(lidres$orf_name) &
                   ((lidres$orf_name %in% c(g1_orf_name) &
                       lidres$density %in% c(g1_density) &
                       lidres$rep %in% c(g1_rep) &
                       lidres$arm %in% c(g1_arm)) |
                      (lidres$orf_name %in% c(g2_orf_name) &
                         lidres$density %in% c(g2_density) &
                         lidres$rep %in% c(g2_rep) &
                         lidres$arm %in% c(g2_arm))),]
  
  box_plot <- ggplot(data,
         aes(x = orf_name, y = fitness, fill = orf_name)) +
    geom_hline(yintercept = 1, col = "red", linetype = 'dashed', lwd = 1) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "ORF",
         y = 'Fitness') +
    scale_fill_discrete(name = 'ORF',
                        breaks = sort(unique(lidres$orf_name)),
                        drop = F,
                        guide = F) +
    theme_linedraw() +
    facet_wrap(.~arm*density*rep,
               nrow = 1) +
    theme(plot.title = element_blank(),
          axis.title = element_text(size = titles),
          axis.text = element_text(size = txt),
          axis.text.x = element_text(angle = 90),
          legend.title = element_blank(),
          legend.text = element_text(size = txt),
          legend.key.size = unit(3, "mm"),
          legend.position = "bottom",
          legend.box.spacing = unit(0.5,"mm"),
          strip.text = element_text(size = txt,
                                    margin = margin(0.1,0,0.1,0, "mm"))) +
    coord_cartesian(ylim = c(0.8, 1.2))
  
  return(box_plot)
}

