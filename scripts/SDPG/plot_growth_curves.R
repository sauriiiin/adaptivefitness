plot_growth_curves <- function(a, r, o, gc, n_col) {

  load(file = '/home/sbp29/R/Projects/adaptivefitness/figs/SDPG/all_methods.RData')
  
  titles <- 7
  txt <- 7
  lbls <- 9
  
  melt1 <- melt(all_results[all_results$condition %in% a & all_results$exp_rep %in% r & all_results$orf_name %in% c(o,'BF_control'),1:10],
                id.vars = c('condition', 'exp_rep', 'orf_name', 'bio_rep'),
                variable.name = 'method', value.name = 'dtime')
  melt2 <- melt(all_results[all_results$condition %in% a & all_results$exp_rep %in% r & all_results$orf_name %in% c(o,'BF_control'),c(1:4,13,16,19,22,25,28)],
                id.vars = c('condition', 'exp_rep', 'orf_name', 'bio_rep'),
                variable.name = 'method2', value.name = 'phenotype')
  all_results_mini <- cbind(melt1, melt2[,c(5,6)])
  
  sum_plot <- ggplot(all_results_mini,
                     aes(x = orf_name, y = dtime)) +
    geom_point(aes(col = phenotype), size = 3) +
    geom_text(aes(label = bio_rep), col = 'white', size = 2) +
    scale_color_manual(name = 'Phenotype',
                       breaks = c('Beneficial','Neutral','Deleterious'),
                       values = c('Deleterious'='#3F51B5',
                                  'Neutral'='#212121',
                                  'Beneficial'='#FFC107'),
                       na.value = 'grey50') +
    labs(x = 'ORFs',
         y = 'Doubling Time (mins)') +
    facet_wrap(.~condition*exp_rep*method, ncol = 6) +
    theme_linedraw() +
    theme(plot.title = element_text(size = titles + 1, hjust = 0.5, face = 'bold'),
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
  
  if (gc == 'Y') {
    raw_plot <- ggplot(gr_data[gr_data$arm %in% a & gr_data$replicate %in% r & gr_data$orf_name %in% c(o,'BF_control'),]) +
      geom_point(aes(x = time, y = od, col = as.factor(bio_rep))) +
      scale_color_discrete(name = 'Replicate') +
      coord_cartesian(xlim = c(0,4000),
                      ylim = c(0.01,6)) +
      # scale_y_log10() +
      labs(title = 'Raw Data',
           x = 'Time (min)',
           y = 'PLC OD') +
      facet_wrap(.~arm*replicate*orf_name,
                 ncol = n_col) +
      theme_linedraw() +
      theme(plot.title = element_text(size = titles + 1, hjust = 0.5, face = 'bold'),
            axis.title = element_text(size = titles),
            axis.text = element_text(size = txt),
            legend.title = element_text(size = titles),
            legend.text = element_text(size = txt),
            legend.key.size = unit(3, "mm"),
            legend.position = "bottom",
            legend.box.spacing = unit(0.5,"mm"),
            strip.text = element_text(size = txt,
                                      margin = margin(0.1,0,0.1,0, "mm")))
    
    gr_plot <- ggplot() +
      geom_point(data = gr_data[gr_data$arm %in% a & gr_data$replicate %in% r & gr_data$orf_name %in% c(o,'BF_control'),],
                 aes(x = time, y = od, col = as.factor(bio_rep)), alpha = 0.2) +
      geom_line(data = gr_data[gr_data$arm %in% a & gr_data$replicate %in% r & gr_data$orf_name %in% c(o,'BF_control'),],
                aes(x = time, y = y, col = as.factor(bio_rep))) +
      geom_point(data = gr_exp[gr_exp$arm %in% a & gr_exp$replicate %in% r & gr_exp$orf_name %in% c(o,'BF_control'),],
                 aes(x = time, y = y), col = 'black', shape = 1) +
      # geom_line(data = gr_fit[gr_fit$arm %in% a & gr_fit$replicate %in% r & gr_fit$orf_name %in% c(o,'BF_control'),],
      #           aes(x = time, y = exp(y), col = as.factor(bio_rep))) +
      coord_cartesian(xlim = c(0,4000),
                      ylim = c(0.01,6)) +
      scale_color_discrete(name = 'Replicate') +
      # scale_y_log10() +
      labs(title = 'Growthrates Results',
           x = 'Time (min)',
           y = 'PLC OD') +
      facet_wrap(.~arm*replicate*orf_name, ncol = n_col) +
      theme_linedraw() +
      theme(plot.title = element_text(size = titles + 1, hjust = 0.5, face = 'bold'),
            axis.title = element_text(size = titles),
            axis.text = element_text(size = txt),
            legend.title = element_text(size = titles),
            legend.text = element_text(size = txt),
            legend.key.size = unit(3, "mm"),
            legend.position = "bottom",
            legend.box.spacing = unit(0.5,"mm"),
            strip.text = element_text(size = txt,
                                      margin = margin(0.1,0,0.1,0, "mm")))
    
    gc_plot <- ggplot(gc_results[gc_results$arm %in% a & gc_results$replicate %in% r & gc_results$orf_name %in% c(o,'BF_control'),]) +
      geom_point(aes(x = time, y = od, col = as.factor(bio_rep)), alpha = 0.2) +
      geom_line(aes(x = time, y = pred.od, col = as.factor(bio_rep))) +
      coord_cartesian(xlim = c(0,4000),
                      ylim = c(0.01,6)) +
      scale_color_discrete(name = 'Replicate') +
      # scale_y_log10() +
      labs(title = 'Growthcurver Results',
           x = 'Time (min)',
           y = 'PLC OD') +
      facet_wrap(.~arm*replicate*orf_name, ncol = n_col) +
      theme_linedraw() +
      theme(plot.title = element_text(size = titles + 1, hjust = 0.5, face = 'bold'),
            axis.title = element_text(size = titles),
            axis.text = element_text(size = txt),
            legend.title = element_text(size = titles),
            legend.text = element_text(size = txt),
            legend.key.size = unit(3, "mm"),
            legend.position = "bottom",
            legend.box.spacing = unit(0.5,"mm"),
            strip.text = element_text(size = txt,
                                      margin = margin(0.1,0,0.1,0, "mm")))
    
    all_plot <- ggpubr::ggarrange(raw_plot, gr_plot, gc_plot,
                                      nrow = 1,
                                      common.legend = T,
                                      legend = 'bottom')
    # all_plot <- ggpubr::ggarrange(sum_plot, all_plot_bot,
    #                               ncol = 1,
    #                               heights = c(1.5,2))
  } else {
    all_plot <- sum_plot
  }
  
  return(all_plot)
}