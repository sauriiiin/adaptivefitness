##### OVEREXPRESSION SCREEN ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 10/28/2021 

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(ggridges)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(ggforce)
library(plotly)
library(scales)
library(reshape2)
library(locfit)
library(growthcurver)
library(rstatix)
library(gtools)
library(growthrates)
library(RMariaDB)
library(genefilter)
library(apeglm)
library(clusterProfiler)
library(org.Sc.sgd.db)

out_path <- "~/R/Projects/adaptivefitness/output/translatome/"
fig_path <- "~/R/Projects/adaptivefitness/figs/translatome/"

expt.name <- "tr_oe_all"

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

`%notin%` <- Negate(`%in%`)

# load('output/translatome/AllData.RData')

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### INITIALIZE
controls <- read.csv(file = 'rawdata/translatome/condition_controls.csv', stringsAsFactors = F)
pgs <- dbGetQuery(conn, 'select orf_name from PROTOGENES where pg_2012 = 1')

orf_types <- read.csv(file = 'rawdata/translatome/oe_transient')
orf_types$category[orf_types$is_transient + orf_types$translated + orf_types$is_candidate == 3] <- 'Transient'
orf_types$category[orf_types$is_conserved + orf_types$translated + orf_types$is_candidate == 3] <- 'Conserved'
orf_types$category[is.na(orf_types$category)] <- 'Others'

tr.conds <- data.frame(arms = c('ONE','ONE','ONE','ONE','TWO','TWO','TWO'),
                       conds = c('GA','SA','HO','HU','DM','FL','TN'))

borders <- dbGetQuery(conn, 'select * from TR_OE_borderpos')
smudge <- dbGetQuery(conn, 'select * from TR_OE_ONE_FS_GA_smudgebox')


##### GATHER DATA
data.fit.all <- NULL
for (a in unique(tr.conds$arms)) {
  for (c in tr.conds$conds[tr.conds$arms == a]) {
    temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate, b.col, b.row
                                         from TR_OE_ALL_FS_%s_%s_6144_FITNESS a, TR_OE_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate, b.col, b.row', a, c))
    temp.fit$arm <- a
    temp.fit$condition <- c
    temp.fit$saturation <- max(temp.fit$hours)
    data.fit.all <- rbind(data.fit.all, temp.fit)
  }
}
data.fit.all$average[data.fit.all$pos %in% borders$pos] <- NA
data.fit.all$average[data.fit.all$condition == 'HO' & data.fit.all$plate %in% c(2,3,10,18,19,23,26,32)] <- NA
data.fit.all$average[data.fit.all$arm == 'ONE' & data.fit.all$plate == 14] <- NA
data.fit.all$average[data.fit.all$arm == 'ONE' & data.fit.all$pos %in% smudge$pos] <- NA

data.fit.all <- merge(data.fit.all, orf_types, by = 'orf_name', all.x = T)
data.fit.all$category[data.fit.all$orf_name == 'BF_control'] <- 'Reference'
# data.fit.all[is.na(data.fit.all$category),] %>%
#   group_by(orf_name) %>%
#   count() %>% data.frame()
data.fit.all$rep <- as.numeric(str_trunc(as.character(data.fit.all$pos), 5, side = 'left', ellipsis = ''))

data.fit.all <- data.fit.all[!(data.fit.all$condition == 'FL' & data.fit.all$hours == 192) &
                               !(data.fit.all$condition == 'GA' & data.fit.all$hours == 67) &
                               !(data.fit.all$condition == 'SA' & data.fit.all$hours %in% c(51,63,72)),]

data.fit.sum <- data.fit.all %>%
  group_by(arm, condition, hours, strain_id, orf_name, category) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

##### REMOVE OUTLIERS
data.fit.all <- merge(data.fit.all, data.fit.sum,
                      by = c('arm','condition','hours','strain_id','orf_name','category'), all = T)

data.fit.all$average[data.fit.all$average < (data.fit.all$avg.median - 2*data.fit.all$avg.mad) |
                       data.fit.all$average > (data.fit.all$avg.median + 2*data.fit.all$avg.mad)] <- NA
data.fit.all$fitness[data.fit.all$fitness < (data.fit.all$fitness.median - 2*data.fit.all$fitness.mad) |
                       data.fit.all$fitness > (data.fit.all$fitness.median + 2*data.fit.all$fitness.mad)] <- NA

data.fit.sum <- data.fit.all %>%
  group_by(arm, condition, hours, strain_id, orf_name, category) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
head(data.fit.sum)

##### REFERENCE LIMITS
data.fit.lim <- data.fit.all %>%
  filter(orf_name == 'BF_control') %>%
  group_by(arm, condition, hours, orf_name, rep) %>%
  summarize(average = median(average, na.rm = T),
            fitness = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  group_by(arm, condition, hours, orf_name) %>%
  summarize(avg_ll = quantile(average, 0.025, na.rm = T),
            avg_m = median(average, na.rm = T),
            avg_ul = quantile(average, 0.975, na.rm = T),
            fitness_ll = quantile(fitness, 0.025, na.rm = T),
            fitness_m = median(fitness, na.rm = T),
            fitness_ul = quantile(fitness, 0.975, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.fit.sum <- merge(data.fit.sum, data.fit.lim[,-4], by = c('arm','condition','hours'))
data.fit.sum$norm_cs <- data.fit.sum$fitness.median * data.fit.sum$avg_m

##### GROWTH CURVE ANALYSIS
data.fit.cs <- merge(data.fit.all[,c(1:21)], data.fit.lim[,-4], by = c('arm','condition','hours'))

data.temp <- data.fit.cs %>%
  group_by(arm, condition, hours, rep, strain_id, orf_name, category) %>%
  summarize(fitness = median(fitness, na.rm = T),
            avg_m = median(avg_m, na.rm = T), .groups = 'keep') %>%
  data.frame()
data.temp$norm_cs <- data.temp$fitness * data.temp$avg_m

data.temp <- data.temp %>% #[,c('arm', 'condition', 'pos', 'rep', 'strain_id', 'orf_name', 'category', 'hours', 'norm_cs')]
  group_by(arm, condition, rep, strain_id, orf_name, category, hours) %>%
  summarize(norm_cs = norm_cs, .groups = 'keep') %>%
  data.frame()

data.strains <- data.temp %>%
  group_by(arm, condition, rep, strain_id, orf_name, category) %>%
  count() %>% data.frame()

# data.strains %>%
#   group_by(arm, condition, category) %>%
#   count() %>% data.frame()

data.gr.res <- NULL
data.gc.res <- NULL
for (a in unique(data.temp$arm)) {
  for (c in unique(data.temp$condition[data.temp$arm == a])) {
    data.pred <- NULL
    col.names <- NULL
    for (r in unique(data.temp$rep[data.temp$arm == a & data.temp$condition == c])) {
      for (p in unique(data.temp$pos[data.temp$arm == a & data.temp$condition == c & data.temp$rep == r])) {
        for (s in unique(data.temp$strain_id[data.temp$arm == a & data.temp$condition == c & data.temp$rep == r &
                                             data.temp$pos == p])) {
          temp <- data.temp[data.temp$arm == a &
                              data.temp$condition == c &
                              data.temp$rep == r &
                              data.temp$strain_id == s &
                              data.temp$pos == p,]
          if (sum(is.na(temp$norm_cs)) <= 5) {
            lo <- loess.smooth(temp$hours, log(temp$norm_cs),
                               span = 0.6, evaluation = 100, degree = 2,
                               family = 'gaussian')
            data.pred <- cbind(data.pred,exp(lo$y))
            # data.pred <- cbind(data.pred, temp$norm_cs)
            col.names <- cbind(col.names,paste(a,c,r,p,s,sep = ','))

            # temp.plot <- ggplot() +
            #   geom_line(data = data.frame(lo), aes(x = x, y = y)) +
            #   geom_point(data = temp, aes(x = hours, y = log(norm_cs)), shape = 1) +
            #   labs(y = 'Log( Colony Size (pixels) )',
            #        x = 'Time (hours)',
            #        title = paste(a,c,r,s,sep = ' | ')) +
            #   theme_linedraw() +
            #   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
            #         axis.title = element_text(size = titles),
            #         axis.text = element_text(size = txt),
            #         legend.title = element_text(size = titles),
            #         legend.text = element_text(size = txt),
            #         legend.position = 'bottom',
            #         legend.key.size = unit(3, "mm"),
            #         legend.box.spacing = unit(0.5,"mm"),
            #         strip.text = element_text(size = txt,
            #                                   face = 'bold',
            #                                   margin = margin(0.1,0,0.1,0, "mm"))) +
            #   coord_cartesian(ylim = c(3,9))
            # ggsave(sprintf("%s/%s/growthcurves/%s_%s_%s.jpg",fig_path, expt.name, o, c, b), temp.plot,
            #        height = 100, width = 100, units = 'mm',
            #        dpi = 600)
          }
        }
      }
    }
    data.pred <- cbind(lo$x, data.pred)
    # data.pred <- cbind(temp$hours, data.pred)
    data.pred <- data.frame(data.pred)
    colnames(data.pred) <- c('Time',col.names)
    head(data.pred)

    ## GROWTH CURVE ANALYSIS
    temp.gr.res <- NULL
    for (i in 2:dim(data.pred)[2]) {
      fit0 <- fit_easylinear(data.pred$Time[2:dim(data.pred)[1]], data.pred[2:dim(data.pred)[1],i], h = 8, quota = 1);

      temp_res <- data.frame(colnames(data.pred[i]), maxgr = coef(fit0)[[3]],
                             dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
      temp.gr.res <- rbind(temp.gr.res, temp_res)
    }
    temp.gr.res <- data.frame(temp.gr.res)
    colnames(temp.gr.res) <- c('sample','gr','dtime','lag')
    temp <- str_split(temp.gr.res$sample, ',', simplify = T)
    colnames(temp) <- c('arm','condition','rep','pos','strain_id')
    temp.gr.res <- cbind(temp, temp.gr.res)

    temp.gc.res <- SummarizeGrowthByPlate(data.pred)
    temp <- str_split(temp.gc.res$sample, ',', simplify = T)
    colnames(temp) <- c('arm','condition','rep','pos','strain_id')
    temp.gc.res <- cbind(temp, temp.gc.res)

    data.gr.res <- rbind(data.gr.res, temp.gr.res)
    data.gc.res <- rbind(data.gc.res, temp.gc.res)
  }
}
data.gcgr.res <- merge(data.gc.res, data.gr.res, by = c('arm','condition','rep','pos','strain_id','sample'))
data.gcgr.res <- merge(data.strains[,-7], data.gcgr.res, by = c('arm','condition','rep','pos','strain_id'))


##### REFERENCE LIMITS
head(data.gcgr.res)
data.gc.lim <- data.gcgr.res[,c(-3,-4,-6,-7,-16)] %>%
  filter(orf_name == 'BF_control') %>%
  # group_by(arm, condition, orf_name) %>%
  melt(id.vars = c('arm','condition','orf_name')) %>%
  group_by(arm, condition, orf_name, variable) %>%
  summarize(ll = quantile(value, 0.005, na.rm = T),
            m = median(value, na.rm = T),
            ul = quantile(value, 0.995, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.gcgr.sum <- melt(data.gcgr.res[!is.na(data.gcgr.res$category),c(1:14,17:19)] %>%
       filter(category != 'Reference'), 
       id.vars = c('arm','condition','rep','strain_id','orf_name','category','sample'))
data.gcgr.sum <- merge(data.gcgr.sum, data.gc.lim[,-3], by = c('arm','condition','variable'), all.x = T)

data.gcgr.sum$es <- (data.gcgr.sum$value - data.gcgr.sum$m)/data.gcgr.sum$m * 100

data.gcgr.sum$annotation[str_detect(data.gcgr.sum$orf_name, 'smo')] <- 'Unannotated'
data.gcgr.sum$annotation[is.na(data.gcgr.sum$annotation)] <- 'Annotated'

##### SAVING FITNESS AND GROWTH PARAMETERS DATA
data.fit.cs <- merge(data.fit.all[,c(1:21)], data.fit.lim[,-4], by = c('arm','condition','hours'))
data.temp <- data.fit.cs %>%
  group_by(arm, condition, hours, rep, strain_id, orf_name, category) %>%
  summarize(fitness = median(fitness, na.rm = T),
            cs = median(average, na.rm = T),
            control_cs_median = median(avg_m, na.rm = T), .groups = 'keep') %>%
  data.frame()
data.temp$norm_cs <- data.temp$fitness * data.temp$control_cs_median

write.csv(data.gcgr.res, file = 'output/translatome/tr_oe_growth_parameters.csv')
write.csv(data.temp, file = 'output/translatome/tr_oe_fitandcs_highres.csv')

##### FALSE DISCOVERY
data.gcgr.fdr <- merge(merge(data.gcgr.sum %>%
                               filter(variable == 'gr', value > m, category == 'Transient') %>%
                               group_by(arm, condition, category, annotation) %>%
                               count() %>% data.frame(),
                             data.gcgr.sum %>%
                               filter(variable == 'gr', value > ul, category == 'Transient') %>%
                               group_by(arm, condition, category, annotation) %>%
                               count() %>% data.frame(),
                             by = c('arm','condition','category','annotation'), 
                             suffixes = c('_m','_ul')),
                       merge(data.gcgr.sum %>%
                               filter(variable == 'gr', value > ul, category == 'Transient', es > 10) %>%
                               group_by(arm, condition, category, annotation) %>%
                               count() %>% data.frame(),
                             data.gcgr.sum %>%
                               filter(variable == 'gr', value > ul, category == 'Transient', es > 30) %>%
                               group_by(arm, condition, category, annotation) %>%
                               count() %>% data.frame(),
                             by = c('arm','condition','category','annotation'),
                             suffixes = c('_es10','_es30'), all = T),
                       by = c('arm','condition','category','annotation'), all.x = T)
data.gcgr.fdr$fp <- data.gcgr.fdr$n_m * 0.005
data.gcgr.fdr$fdr <- data.gcgr.fdr$fp/data.gcgr.fdr$n_ul * 100
data.gcgr.fdr


# ##### ENRICHMENT ANALYSIS
# goe <- NULL
# kegg <- NULL
# for (a in unique(data.gcgr.sum$arm[data.gcgr.sum$variable == 'gr' &
#                 data.gcgr.sum$annotation == 'Annotated' &
#                 data.gcgr.sum$strain_id < 1000000])) {
#   for (c in unique(data.gcgr.sum$condition[data.gcgr.sum$variable == 'gr' &
#                                      data.gcgr.sum$annotation == 'Annotated' &
#                                      data.gcgr.sum$strain_id < 1000000 &
#                                      data.gcgr.sum$arm == a])) {
#     temp <- data.gcgr.sum[data.gcgr.sum$variable == 'gr' &
#                                      data.gcgr.sum$annotation == 'Annotated' &
#                             data.gcgr.sum$category != 'Transient' &
#                                      data.gcgr.sum$strain_id < 1000000 &
#                                      data.gcgr.sum$arm == a &
#                                      data.gcgr.sum$condition == c,]
#     
#     allgenes <- unique(temp$orf_name)
#     allgenes <- bitr(allgenes, fromType = "ORF",
#                      toType = c("ENTREZID","GENENAME","ENSEMBL"),
#                      OrgDb = org.Sc.sgd.db)
#     
#     temp.deg <- head(temp$orf_name[order(-temp$value)],200)
#     temp.deg <- bitr(temp.deg, fromType = "ORF",
#                      toType = c("ENTREZID","GENENAME","ENSEMBL"),
#                      OrgDb = org.Sc.sgd.db)
#     
#     temp.goe <- enrichGO(gene          = temp.deg$ENSEMBL,
#                          universe      = allgenes$ENSEMBL,
#                          OrgDb         = org.Sc.sgd.db,
#                          keyType       = "ENSEMBL",
#                          ont           = "ALL",
#                          pAdjustMethod = "BH",
#                          pvalueCutoff  = 0.01,
#                          qvalueCutoff  = 0.05)
#     if (length(temp.goe) == 1) {
#       if (dim(temp.goe)[1] == 0) {
#       } else{
#         goe <- rbind(goe, data.frame(arm = a, condition = c, variable = 'gr', phenotype = 'Beneficial', temp.goe))
#       }
#     }
#     
#     temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
#                             universe     = allgenes$ENSEMBL,
#                             organism     = 'sce',
#                             pvalueCutoff = 0.05)
#     if (length(temp.kegg) == 1) {
#       if (dim(temp.kegg)[1] == 0) {
#       } else{
#         kegg <- rbind(kegg, data.frame(arm = a, condition = c, variable = 'gr', phenotype = 'Beneficial', temp.kegg))
#       }
#     }
#   }
# }
# goe$GeneRatio <- as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,1])/
#   as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,2])
# goe$BgRatio <- as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,1])/
#   as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,2])
# goe$GO <- paste0(goe$ONTOLOGY, '_', goe$Description)
# 
# kegg$GeneRatio <- as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,1])/
#   as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,2])
# kegg$BgRatio <- as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,1])/
#   as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,2])
# 
# goe <- goe[order(goe$arm,goe$condition,goe$phenotype,-goe$GeneRatio,-goe$Count,goe$qvalue),]
# kegg <- kegg[order(kegg$arm,kegg$condition,kegg$phenotype,-kegg$GeneRatio,-kegg$Count,kegg$qvalue),]

##### DIFFERENTIAL VALUE
data.gcgr.sum2 <- rbind(merge(data.gcgr.sum[,-c(8,10:13)] %>%
              filter(arm == 'ONE', condition %in% c('HO','HU','SA')),
            data.gcgr.sum[,-c(8,10:13)] %>%
              filter(arm == 'ONE', condition %in% c('GA')),
            by = c('arm','variable','rep','strain_id','orf_name','category','annotation'),
            suffixes = c('_stress','_cont')),
      merge(data.gcgr.sum[,-c(8,10:13)] %>%
              filter(arm == 'TWO', condition %in% c('FL','TN')),
            data.gcgr.sum[,-c(8,10:13)] %>%
              filter(arm == 'TWO', condition %in% c('DM')),
            by = c('arm','variable','rep','strain_id','orf_name','category','annotation'),
            suffixes = c('_stress','_cont')))
data.gcgr.sum2$value <- (data.gcgr.sum2$value_stress - data.gcgr.sum2$value_cont)/data.gcgr.sum2$value_cont * 100
head(data.gcgr.sum2)


data.gcgr.ref2 <- data.gcgr.res[,c(-4,-6,-7,-16)] %>%
  filter(orf_name == 'BF_control') %>%
  melt(id.vars = c('arm','condition','orf_name','rep')) %>%
  data.frame()

data.gcgr.ref2 <- rbind(merge(data.gcgr.ref2 %>%
              filter(arm == 'ONE', condition %in% c('HO','HU','SA')),
              data.gcgr.ref2 %>%
              filter(arm == 'ONE', condition %in% c('GA')),
            by = c('arm','variable','rep','orf_name'),
            suffixes = c('_stress','_cont')),
      merge(data.gcgr.ref2 %>%
              filter(arm == 'TWO', condition %in% c('FL','TN')),
            data.gcgr.ref2 %>%
              filter(arm == 'TWO', condition %in% c('DM')),
            by = c('arm','variable','rep','orf_name'),
            suffixes = c('_stress','_cont')))
head(data.gcgr.ref2)

data.gcgr.lim2 <- NULL
for (a in unique(data.gcgr.ref2$arm)) {
  for (s in unique(data.gcgr.ref2$condition_stress[data.gcgr.ref2$arm == a])) {
    for (c in unique(data.gcgr.ref2$condition_cont[data.gcgr.ref2$arm == a &
                                                   data.gcgr.ref2$condition_stress == s])) {
      for (v in unique(data.gcgr.ref2$variable[data.gcgr.ref2$arm == a &
                                               data.gcgr.ref2$condition_stress == s &
                                               data.gcgr.ref2$condition_cont == c])) {
        temp.vstress <- data.gcgr.ref2$value_stress[data.gcgr.ref2$arm == a &
                                                      data.gcgr.ref2$condition_stress == s &
                                                      data.gcgr.ref2$condition_cont == c &
                                                      data.gcgr.ref2$variable == v]
        temp.vcont <- data.gcgr.ref2$value_cont[data.gcgr.ref2$arm == a &
                                                      data.gcgr.ref2$condition_stress == s &
                                                      data.gcgr.ref2$condition_cont == c &
                                                      data.gcgr.ref2$variable == v]
        temp.q <- NULL
        set.seed(8784334)
        for (i in 1:10) {
          temp1 <- sample(temp.vstress)
          temp2 <- sample(temp.vcont)
          temp.q <- c(temp.q, (temp1 - temp2)/temp2 * 100)
        }
        
        temp.q <- quantile(temp.q, c(0.005,0.5,0.995), na.rm = T)
        
        data.gcgr.lim2 <- rbind(data.gcgr.lim2,
                                data.frame(arm = a, condition_stress = s, condition_cont = c, variable = v, 
                                           ll = temp.q[[1]], m = temp.q[[2]], ul = temp.q[[3]]))
        
      }
    }
  }
}
head(data.gcgr.lim2)

data.gcgr.sum2 <- merge(data.gcgr.sum2[,c(1:12)], data.gcgr.lim2,
                        by = c('arm','condition_stress','condition_cont','variable'))
data.gcgr.sum2$es <- (data.gcgr.sum2$value - data.gcgr.sum2$m)/abs(data.gcgr.sum2$m) * 100

##### FALSE DISCOVERY USING DIFFERENTIAL
data.gcgr.fdr2 <- merge(merge(data.gcgr.sum2 %>%
                               filter(variable == 'gr', value > m, category == 'Transient') %>%
                               group_by(arm, condition_stress, category, annotation) %>%
                               count() %>% data.frame(),
                             data.gcgr.sum2 %>%
                               filter(variable == 'gr', value > ul, category == 'Transient') %>%
                               group_by(arm, condition_stress, category, annotation) %>%
                               count() %>% data.frame(),
                             by = c('arm','condition_stress','category','annotation'), 
                             suffixes = c('_m','_ul')),
                       merge(data.gcgr.sum2 %>%
                               filter(variable == 'gr', value > ul, category == 'Transient', es > 10) %>%
                               group_by(arm, condition_stress, category, annotation) %>%
                               count() %>% data.frame(),
                             data.gcgr.sum2 %>%
                               filter(variable == 'gr', value > ul, category == 'Transient', es > 30) %>%
                               group_by(arm, condition_stress, category, annotation) %>%
                               count() %>% data.frame(),
                             by = c('arm','condition_stress','category','annotation'),
                             suffixes = c('_es10','_es30'), all = T),
                       by = c('arm','condition_stress','category','annotation'), all.x = T)
data.gcgr.fdr2$fp <- data.gcgr.fdr2$n_m * 0.005
data.gcgr.fdr2$fdr <- data.gcgr.fdr2$fp/data.gcgr.fdr2$n_ul * 100
data.gcgr.fdr2


##### DIFFERENTIAL ENRICHMENT ANALYSIS
goe2 <- NULL
kegg2 <- NULL
for (a in unique(data.gcgr.sum2$arm)) {
  for (c in unique(data.gcgr.sum2$condition_stress[data.gcgr.sum2$arm == a])) {
    for (v in c('gr')) {
      temp <- data.gcgr.sum2[data.gcgr.sum2$variable == v &
                               data.gcgr.sum2$annotation == 'Annotated' &
                               data.gcgr.sum2$category != 'Transient' &
                               data.gcgr.sum2$strain_id < 1000000 &
                               data.gcgr.sum2$arm == a &
                               data.gcgr.sum2$condition_stress == c,]
      
      allgenes <- unique(temp$orf_name)
      allgenes <- bitr(allgenes, fromType = "ORF",
                       toType = c("ENTREZID","GENENAME","ENSEMBL"),
                       OrgDb = org.Sc.sgd.db)
      
      temp.deg <- head(temp$orf_name[order(-temp$value)],50)
      temp.deg <- bitr(temp.deg, fromType = "ORF",
                       toType = c("ENTREZID","GENENAME","ENSEMBL"),
                       OrgDb = org.Sc.sgd.db)
      
      temp.goe <- enrichGO(gene          = temp.deg$ENSEMBL,
                           universe      = allgenes$ENSEMBL,
                           OrgDb         = org.Sc.sgd.db,
                           keyType       = "ENSEMBL",
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
      if (length(temp.goe) == 1) {
        if (dim(temp.goe)[1] == 0) {
        } else{
          goe2 <- rbind(goe2, data.frame(arm = a, condition = c, variable = v, phenotype = 'Beneficial', temp.goe))
        }
      }
      
      temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                              universe     = allgenes$ENSEMBL,
                              organism     = 'sce',
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05)
      if (length(temp.kegg) == 1) {
        if (dim(temp.kegg)[1] == 0) {
        } else{
          kegg2 <- rbind(kegg2, data.frame(arm = a, condition = c, variable = v, phenotype = 'Beneficial', temp.kegg))
        }
      }
    }
  }
}
goe2$GeneRatio <- as.numeric(str_split(goe2$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe2$GeneRatio,'/',simplify = T)[,2])
goe2$BgRatio <- as.numeric(str_split(goe2$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe2$BgRatio,'/',simplify = T)[,2])
goe2$GO <- paste0(goe$ONTOLOGY, '_', goe2$Description)

kegg2$GeneRatio <- as.numeric(str_split(kegg2$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg2$GeneRatio,'/',simplify = T)[,2])
kegg2$BgRatio <- as.numeric(str_split(kegg2$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg2$BgRatio,'/',simplify = T)[,2])

goe2 <- goe2[order(goe2$arm,goe2$condition,goe2$variable,goe2$phenotype,-goe2$GeneRatio,-goe2$Count,goe2$qvalue),]
kegg2 <- kegg2[order(kegg2$arm,kegg2$condition,kegg2$variable,kegg2$phenotype,-kegg2$GeneRatio,-kegg2$Count,kegg2$qvalue),]

write.csv(goe2, file = 'output/translatome/go_gc_enrichments.csv')
# write.csv(kegg, file = 'output/translatome/kegg_enrichments.csv')
