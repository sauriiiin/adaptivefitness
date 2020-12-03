##### EMPIRICAL HYPOTHESIS TESTING
## eg. ~/R/Projects/adaptivefitness/scripts/paper/FIGURES_REV.R

empirical_p <- function(data, ref, value, p, es, avoid, fpr) {
  for (h in sort(unique(data$hours))) {
    ref_fit <- data[data$hours == h & data$orf_name == ref, value]
    for (o in sort(unique(data$orf_name[!(data$orf_name %in% avoid)]))) {
      orf_fit <- data[data$hours == h & data$orf_name == o, value]
      es_val <- (orf_fit - mean(ref_fit, na.rm = T))/mean(ref_fit, na.rm = T) * 100
      
      p_val1 <- sum(orf_fit < ref_fit)/length(ref_fit)
      p_val2 <- sum(orf_fit > ref_fit)/length(ref_fit)
      p_val <- min(c(p_val1, p_val2))*2
      data[data$hours == h & data$orf_name == o, p] <- p_val
      data[data$hours == h & data$orf_name == o, es] <- es_val
    }
    data[data$hours == h, sprintf('%s_thresh',p)] <- quantile(data[data$hours == h, p], probs = fpr, na.rm = T)
  }
  return(data)
}