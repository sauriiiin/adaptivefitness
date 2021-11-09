##### ALL F28 f28.dat
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 04/07/2021

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
library(RMariaDB)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(PoiClaClu)
library(RColorBrewer)

source("R/functions/isoutlier.R")
source("R/functions/empirical_p.R")
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
out_path = 'figs/f28fu/repeat/'

##### FIGURE SIZE
one.c <- 90 #single column 
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

##### LOAD DATA
load("/home/sbp29/R/Projects/adaptivefitness/figs/f28fu/repeat/f28alldata.RData")
head(f28.dat)

f28.dat$rep <- as.numeric(str_trunc(as.character(f28.dat$pos), 4, side = 'left', ellipsis = ''))

##### EMPIRICAL P-VALUES
temp <- NULL
ref <- 'BF_control'
avoid <- c(ref,'BOR','REF','YHR021W-A','YNR015W')
i <- 1

for (p in unique(f28.dat$plasmid[f28.dat$location == "Pittsburgh"])) {
  for (s in unique(f28.dat$stage[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p])) {
    for (c in unique(f28.dat$arm[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p & f28.dat$stage == s])) {
      for (d in unique(f28.dat$density[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p & f28.dat$stage == s & f28.dat$arm == c])) {
        for (h in sort(unique(f28.dat$hours[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p & f28.dat$stage == s & f28.dat$arm == c & f28.dat$density == d]))) {
          temp_fit <- f28.dat$fitness[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p & f28.dat$stage == s &
                                        f28.dat$arm == c & f28.dat$density == d & f28.dat$hours == h & 
                                        f28.dat$orf_name == ref & !is.na(f28.dat$fitness)]
          temp_outliers <- isoutlier(temp_fit, 3)
          temp_fit[temp_outliers] <- NA
          f28.dat$fitness[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p & f28.dat$stage == s &
                            f28.dat$arm == c & f28.dat$density == d & f28.dat$hours == h &
                            f28.dat$orf_name == ref & !is.na(f28.dat$fitness)] <- temp_fit
          
          ref_fit <- f28.dat[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p & f28.dat$stage == s &
                               f28.dat$arm == c & f28.dat$density == d & f28.dat$hours == h & 
                               f28.dat$orf_name == ref,] %>%
            group_by(rep, replicate) %>%
            dplyr::summarise(median = median(fitness, na.rm = T), .groups = "keep") %>%
            data.frame()
          ref_fit <- ref_fit$median[!is.na(ref_fit$median)]
          
          for (o in sort(unique(f28.dat$orf_name[!(f28.dat$orf_name %in% avoid) & !is.na(f28.dat$orf_name)]))) {
            temp_fit <- f28.dat$fitness[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p & f28.dat$stage == s &
                                          f28.dat$arm == c & f28.dat$density == d & f28.dat$hours == h & 
                                          f28.dat$orf_name == o & !is.na(f28.dat$fitness)]
            temp_outliers <- isoutlier(temp_fit, 3)
            temp_fit[temp_outliers] <- NA
            f28.dat$fitness[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p & f28.dat$stage == s &
                              f28.dat$arm == c & f28.dat$density == d & f28.dat$hours == h & 
                              f28.dat$orf_name == o & !is.na(f28.dat$fitness)] <- temp_fit
            
            orf_fit <- median(f28.dat$fitness[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p & f28.dat$stage == s &
                                                f28.dat$arm == c & f28.dat$density == d & f28.dat$hours == h & 
                                                f28.dat$orf_name == o], na.rm = T)
            es_val <- (orf_fit - mean(ref_fit, na.rm = T))/mean(ref_fit, na.rm = T) * 100
            
            p_val1 <- sum(orf_fit < ref_fit, na.rm = T)/length(ref_fit)
            p_val2 <- sum(orf_fit > ref_fit, na.rm = T)/length(ref_fit)
            p_val <- min(c(p_val1, p_val2))*2
            
            temp$density[i] <- d
            temp$orf_name[i] <- o
            temp$hours[i] <- h
            # temp$replicate[i] <- f28.dat$plasmid[f28.dat$location == "Pittsburgh" & f28.dat$stage == s &
            #                                        f28.dat$arm == c & f28.dat$density == d & f28.dat$hours == h][1]
            temp$arm[i] <- c
            temp$stage[i] <- s
            temp$location[i] <- "Pittsburgh"
            temp$plasmid[i] <- f28.dat$plasmid[f28.dat$location == "Pittsburgh" & f28.dat$plasmid == p & f28.dat$stage == s &
                                                 f28.dat$arm == c & f28.dat$density == d & f28.dat$hours == h][1]
            temp$p[i] <- p_val
            temp$es[i] <- es_val
            
            i <- i + 1
          }
        } 
      }
    }
  }
}
temp <- data.frame(temp, stringsAsFactors = F)
head(temp)

temp[temp$p <= 0.05 & temp$plasmid == 'No',]

##### FITNESS SUMMARY
f28.sum <- f28.dat[!(f28.dat$orf_name %in% c('BOR','REF','YHR021W-A','YNR015W')) &
                     f28.dat$hours > 0 & !is.na(f28.dat$orf_name),] %>%
  group_by(location, plasmid, stage, density, arm, hours, orf_name) %>%
  summarise(fitness = mean(fitness, na.rm = T), .groups = "keep") %>%
  data.frame()
head(f28.sum)
f28.sum$id <- apply(f28.sum[,c(1:6)], 1, paste0, collapse="_")

orffitness <- ggplot(f28.sum[f28.sum$location != 'SanDiego',], aes(x = id, y = fitness)) +
  geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,
                          ymin=-Inf,ymax=Inf),
            fill="white") +
  geom_rect(data=NULL,aes(xmin='Pittsburgh_Yes_FS_1536_GLU',xmax=Inf,
                          ymin=-Inf,ymax=Inf),
            fill="grey80") +
  geom_vline(xintercept = seq(0,100), col = '#757575', lwd = 0.3, alpha =0.9) +
  geom_hline(yintercept = seq(0,2,0.1), col = '#757575', lwd = 0.3, alpha =0.9) +
  geom_point(aes(col = stage, shape = as.factor(density))) +
  coord_cartesian(ylim = c(0.8,1.2)) +
  facet_wrap(.~orf_name, ncol = 7) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sF28ALL_ORFFITNESS.jpg",out_path), orffitness,
       height = two.c, width = 400, units = 'mm',
       dpi = 300)  


##### BOXPLOTS
orffitness.box <- f28.dat[f28.dat$hours > 0,] %>%
  ggplot(aes(x = orf_name, y = fitness)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(.~location*plasmid*arm*stage*density*hours) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sF28ALL_ORFFITNESS_BOX.jpg",out_path), orffitness.box,
       height = two.c, width = 400, units = 'mm',
       dpi = 300) 

##### CORRELATIONS
f28.sum %>%
  group_by(location, plasmid, stage, density, arm, hours) %>%
  count() %>%
  data.frame()

f28.mat <- NULL
f28.mat$orf_name <- unique(f28.sum$orf_name)
col.names <- 'orf_name'
for (i in unique(f28.sum$id)) {
  col.names <- c(col.names, i)
  f28.mat <- merge(f28.mat, f28.sum[f28.sum$id == i,c(7,8)], by = 'orf_name')
}
colnames(f28.mat) <- col.names
head(f28.mat)

ids <- c(2,3,4,5,6,7,8,9,11,13,15,17)

mat <- cor(f28.mat[,ids], use = 'complete.obs', method = 'spearman')
mat  

col_fun = colorRamp2(c(0.3, 0.5, 1), c("white", "#FF5252", "#512DA8"))
ha <- HeatmapAnnotation(
  # Location = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,1],
  # Plasmid = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,2],
  Stage = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,3],
  Density = as.numeric(str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,4]),
  Condition  = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,5],
  # Hours = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,6],
  col = list(Location = c("Pittsburgh" = "Black", "SanDiego" = "White"),
             Plasmid = c("Yes" = "Black", "No" = "White"),
             Stage = c("FS" = "#E1BEE7", "PS2" = "#E040FB", "S3" = "#7B1FA2"),
             Density = c("384" = "Navy", "1536" = "#607D8B", "6144" = "#FFA000"),
             Condition = c("GAL" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F", "GLU" = "White", "5FOA" = "Black")),
  gp = gpar(col = "black", lwd = 0.2)
)

ra <- rowAnnotation(
  # Location = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,1],
  # Plasmid = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,2],
  Stage = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,3],
  Density = as.numeric(str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,4]),
  Condition  = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,5],
  # Hours = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,6],
  col = list(Location = c("Pittsburgh" = "Black", "SanDiego" = "White"),
             Plasmid = c("Yes" = "Black", "No" = "White"),
             Stage = c("FS" = "#E1BEE7", "PS2" = "#E040FB", "S3" = "#7B1FA2"),
             Density = c("384" = "Navy", "1536" = "#607D8B", "6144" = "#FFA000"),
             Condition = c("GAL" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F", "GLU" = "White", "5FOA" = "Black")),
  gp = gpar(col = "black", lwd = 0.2),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

hm <- Heatmap(mat,
        name = "Correlation",
        column_title = "The Famous 28 Heatmap",
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



##### SUMMARY PLOT
ids <- c(2,3,4,5,6,7,8,9,11,13,15,17)

plot.sum.violin <- f28.sum[f28.sum$id %in% colnames(f28.mat[,ids]),] %>%
  ggplot(aes(x = id, y = fitness)) +
  geom_violin() +
  # geom_jitter(aes(col = orf_name)) +
  labs(y = 'Fitness') +
  scale_color_discrete(name = 'Strain') +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

plot.sum.density <- f28.sum[f28.sum$id %in% colnames(f28.mat[,ids]),] %>%
  ggplot(aes(x = id, y = 1)) +
  geom_tile(aes(fill = as.factor(density)),col = 'black') +
  geom_text(aes(label = density)) +
  labs(y = 'Density') +
  scale_fill_manual(guide = F,
                    values = c("384" = "Navy", "1536" = "#607D8B", "6144" = "#FFA000")) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

plot.sum.arm <- f28.sum[f28.sum$id %in% colnames(f28.mat[,ids]),] %>%
  ggplot(aes(x = id, y = 1)) +
  geom_tile(aes(fill = arm),col = 'black') +
  geom_text(aes(label = arm)) +
  labs(y = 'Arm') +
  scale_fill_manual(guide = F,
                    values = c("GAL" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F", "GLU" = "White", "5FOA" = "Black")) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

plot.sum.stage <- f28.sum[f28.sum$id %in% colnames(f28.mat[,ids]),] %>%
  ggplot(aes(x = id, y = 1)) +
  geom_tile(aes(fill = stage),col = 'black') +
  geom_text(aes(label = stage)) +
  labs(y = 'Stage') +
  scale_fill_manual(guide = F,
                    values = c("FS" = "#E1BEE7", "PS2" = "#E040FB", "S3" = "#7B1FA2")) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

plot.sum.v <- ggarrange(plot.sum.violin, plot.sum.arm, plot.sum.density, plot.sum.stage,
          ncol = 1,
          heights = c(10,1,1,1))
ggsave(sprintf("%sF28SUM_FITNESS_VIOLIN.jpg",out_path), plot.sum.v,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 300) 



plot.sum.line <- f28.sum[f28.sum$id %in% colnames(f28.mat[,ids]),] %>%
  ggplot(aes(x = id, y = fitness)) +
  geom_violin() +
  geom_line(aes(col = orf_name, group = orf_name)) +
  labs(y = 'Fitness') +
  scale_color_discrete(name = 'Strain',
                       guide = F) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))

plot.sum.l <- ggarrange(plot.sum.line, plot.sum.arm, plot.sum.density, plot.sum.stage,
                        ncol = 1,
                        heights = c(10,1,1,1))
ggsave(sprintf("%sF28SUM_FITNESS_LINE.jpg",out_path), plot.sum.l,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 300)

##### PAIRWISE MATRIX (for seaborn style matrix)
f28.pairs <- NULL
for (g1 in unique(f28.sum$id)) {
  for (g2 in unique(f28.sum$id)) {
    f28.pairs <- rbind(f28.pairs, data.frame(merge(f28.sum[f28.sum$id == g1,c(7,8)],
                                                   f28.sum[f28.sum$id == g2,c(7,8)],
                                                   by = 'orf_name',
                                                   suffixes = c('.g1','.g2')),
                                             group1 = g1, group2 = g2))
                       
  }
}

ids <- c(2,3,4,5,6,7,8,9,11,13,15,17)
# ids <- c(2,3)

pairwise.cor <- f28.pairs[f28.pairs$group1 %in% colnames(f28.mat[,ids]) &
            f28.pairs$group2 %in% colnames(f28.mat[,ids]),] %>%
  ggplot(aes(x = fitness.g1, y = fitness.g2)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  facet_wrap(.~group1*group2,
             ncol = 12) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "right",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sF28SUM_PAIRWISE_CORR.jpg",out_path), pairwise.cor,
       height = 500, width = 500, units = 'mm',
       dpi = 300)
