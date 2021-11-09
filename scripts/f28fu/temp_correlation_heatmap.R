
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

out_path = 'figs/f28fu/repeat2/'

##### FIGURE SIZE
one.c <- 90 #single column 
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

load("/home/sbp29/R/Projects/adaptivefitness/figs/f28fu/repeat/f28fur_fs_data.RData")
data$attempt <- 'original'
f28fur_fs <- data
load("/home/sbp29/R/Projects/adaptivefitness/figs/f28fu/repeat2/f28fur2_fs1_data.RData")
data$attempt <- 'repeat'
f28fur2_fs1 <- data
load("/home/sbp29/R/Projects/adaptivefitness/figs/f28fu/repeat2/f28fur2_fs2_data.RData")
data$attempt <- 'repeat2'
f28fur2_fs2 <- data

##### CORRELATION AT SATURATION
data <- rbind(f28fur_fs[(f28fur_fs$density == 6144 & f28fur_fs$hours == 36 & f28fur_fs$condition != 'GLU') |
                            (f28fur_fs$density == 1536 & f28fur_fs$hours == 60),],
              f28fur2_fs1[(f28fur2_fs1$density == 6144 & f28fur2_fs1$hours == 36) |
                      (f28fur2_fs1$density == 1536 & f28fur2_fs1$hours == 100),],
              f28fur2_fs2[(f28fur2_fs2$density == 6144 & f28fur2_fs2$hours == 36) |
                            (f28fur2_fs2$density == 1536 & f28fur2_fs2$hours == 100),])

data.sum <- data[!(data$orf_name %in% c('BOR','REF','BF_control','YHR021W-A')) & !is.na(data$orf_name),] %>%
  group_by(attempt,condition, density, hours, orf_name) %>%
  summarise(fitness = mean(fitness_clean, na.rm = T),
            emp_p = mean(emp_p, na.rm = T),
            effsize = mean(effsize, na.rm = T),
            p_thresh = mean(emp_p_thresh, na.rm = T),
            .groups = "keep") %>%
  data.frame()
head(data.sum)
data.sum$id <- apply(data.sum[,c(1:4)], 1, paste0, collapse="_")

f28.mat <- NULL
f28.mat$orf_name <- unique(data.sum$orf_name)
col.names <- 'orf_name'
for (i in unique(data.sum$id)) {
  col.names <- c(col.names, i)
  f28.mat <- merge(f28.mat, data.sum[data.sum$id == i,c(5,6)], by = 'orf_name')
}
colnames(f28.mat) <- col.names
head(f28.mat)

ids <- -1
mat <- cor(f28.mat[,ids], use = 'complete.obs', method = 'spearman')
mat  

col_fun = colorRamp2(c(0.3, 0.5, 1), c("white", "#FF5252", "#512DA8"))
ha <- HeatmapAnnotation(
  Attempt = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,1],
  Condition  = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,2],
  Density = as.numeric(str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,3]),
  # Hours = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,6],
  col = list(
    # Location = c("Pittsburgh" = "Black", "SanDiego" = "White"),
             # Plasmid = c("Yes" = "Black", "No" = "White"),
             Attempt = c("original" = "#E1BEE7", "repeat" = "#E040FB", "repeat2" = "#7B1FA2"),
             Density = c("384" = "Navy", "1536" = "#607D8B", "6144" = "#FFA000"),
             Condition = c("GAL" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F", "GLU" = "White", "5FOA" = "Black")),
  gp = gpar(col = "black", lwd = 0.2)
)

ra <- rowAnnotation(
  Attempt = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,1],
  Condition  = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,2],
  Density = as.numeric(str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,3]),
  # Hours = str_split(colnames(f28.mat[,ids]), '_', simplify = T)[,6],
  col = list(
    # Location = c("Pittsburgh" = "Black", "SanDiego" = "White"),
    # Plasmid = c("Yes" = "Black", "No" = "White"),
    Attempt = c("original" = "#E1BEE7", "repeat" = "#E040FB", "repeat2" = "#7B1FA2"),
    Density = c("384" = "Navy", "1536" = "#607D8B", "6144" = "#FFA000"),
    Condition = c("GAL" = "#B3E5FC", "CAS" = "#448AFF", "SDA" = "#303F9F", "GLU" = "White", "5FOA" = "Black")),
  gp = gpar(col = "black", lwd = 0.2),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

hm <- Heatmap(mat,
              name = "Correlation",
              column_title = "The Famous 28 FOA Experiment",
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

##### SIGNIFICANT FITNESS CHANGES
head(data.sum)

sig.o <- data.sum %>%
  filter(emp_p <= 0.05, attempt == 'original') %>%
  group_by(attempt, density, condition) %>%
  data.frame()

sig.r1 <- data.sum %>%
  filter(emp_p <= 0.05, attempt == 'repeat') %>%
  group_by(attempt, density, condition) %>%
  data.frame()

sig.r2 <- data.sum %>%
  filter(emp_p <= 0.05, attempt == 'repeat2') %>%
  group_by(attempt, density, condition) %>%
  data.frame()

sig.all <- merge(merge(sig.o[,c(1,2,3,5)], sig.r1[,c(1,2,3,5)], by = c('condition', 'density', 'orf_name'), all = T),
      sig.r2[,c(1,2,3,5)], by = c('condition', 'density', 'orf_name'), all = T)
colnames(sig.all) <- c('condition','density','orf_name','original','repeat1','repeat2')
is.na(sig.all) <- ''
sig.all$attempts <- str_remove_all(str_remove_all(paste(sig.all$original, sig.all$repeat1, sig.all$repeat2, sep = '\n'), 'NA\n'), '\nNA')
sig.all$attempts <- factor(sig.all$attempts, levels = c('original','repeat2','repeat\nrepeat2','original\nrepeat\nrepeat2')) 

sig.overlap <- sig.all %>%
  ggplot(aes(x = attempts, fill = attempts)) +
  geom_histogram(stat = 'count') +
  geom_text_repel(aes(label = orf_name, y = 1), position = 'identity', size = 2) +
  scale_y_continuous(minor_breaks = NULL) +
  scale_fill_discrete(guide = F) +
  facet_wrap(.~density*condition, nrow = 1) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sSIGNIFICANT_OVERLAP.jpg",out_path), sig.overlap,
       height = 100, width = 300, units = 'mm',
       dpi = 300)


data.sum %>%
  filter(emp_p <= 0.05) %>%
  group_by(orf_name) %>%
  count() %>%
  data.frame()


##### FITNESS BOXPLOTS
data <- rbind(f28fur_fs[(f28fur_fs$density == 6144 & f28fur_fs$hours == 36 & f28fur_fs$condition != 'GLU') |
                          (f28fur_fs$density == 1536 & f28fur_fs$hours == 60),],
              f28fur2_fs1[(f28fur2_fs1$density == 6144 & f28fur2_fs1$hours == 36) |
                            (f28fur2_fs1$density == 1536 & f28fur2_fs1$hours == 100),],
              f28fur2_fs2[(f28fur2_fs2$density == 6144 & f28fur2_fs2$hours == 36) |
                            (f28fur2_fs2$density == 1536 & f28fur2_fs2$hours == 100),])
data$condition <- factor(data$condition, levels = c('GLU','GAL','CAS','SDA'))

all.fit.den <- data[!(data$orf_name %in% c('BOR','REF','BF_control','YHR021W-A')) & !is.na(data$orf_name),] %>%
  # group_by(attempt,condition, density, hours, orf_name,rep) %>%
  # summarise(fitness_clean = mean(fitness_clean, na.rm = T),
  #           .groups = "keep") %>%
  # data.frame() %>%
  ggplot(aes(x = orf_name, y = fitness_clean)) +
  geom_boxplot(aes(fill = attempt), outlier.shape = NA) +
  # geom_density_ridges(quantile_lines = TRUE,
  #                     aes(fill = attempt),
  #                     scale = 2, alpha = 0.7, size = 0.2,
  #                     vline_size = 0.2, vline_color = "black") +
  labs(x = 'ORFs', y = 'Fitness') +
  scale_fill_discrete(name = 'Attempt') +
  facet_wrap(.~density*condition, ncol = 4) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFITNESS_BOX_ALL.jpg",out_path), all.fit.den,
       height = 300, width = 500, units = 'mm',
       dpi = 300)


##### GROWTH CURVES
gc.data <- rbind(f28fur_fs, f28fur2_fs1, f28fur2_fs2)

gc.fitness <- gc.data[!(gc.data$orf_name %in% c('BOR','REF','YHR021W-A')),]%>%
  filter(attempt != 'original') %>%
  ggplot(aes(x = hours, y = fitness_clean, col = orf_name)) +
  geom_smooth(method = 'loess') +
  labs(x = 'Time (hours)', y = 'Fitness') +
  facet_wrap(.~condition*attempt*density, scales = 'free_x') +
  scale_color_discrete(name = 'ORFs') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sFITNESS_GC_ALL.jpg",out_path), gc.fitness,
       height = 300, width = 500, units = 'mm',
       dpi = 300)


gc.cs <- gc.data[!(gc.data$orf_name %in% c('BOR','REF','YHR021W-A')),]%>%
  filter(attempt != 'original') %>%
  ggplot(aes(x = hours, y = average_clean, col = orf_name)) +
  geom_smooth(method = 'loess') +
  labs(x = 'Time (hours)', y = 'Colony Size (pixels)') +
  facet_wrap(.~condition*attempt*density, scales = 'free_x') +
  scale_color_discrete(name = 'ORFs') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = "bottom",
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%sCOLONY_SIZE_GC_ALL.jpg",out_path), gc.cs,
       height = 300, width = 500, units = 'mm',
       dpi = 300)
