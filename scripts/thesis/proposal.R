##### THESIS PROPOSAL
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/11/2020

##### INITIALIZE
library(ggplot2)
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

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
data_path = 'figs/paper/';
out_path = 'figs/proposal/';

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### LID FIGURE
load(sprintf('%sSPEDATA.RData',data_path))
load(sprintf('%sSENDATA.RData',data_path))
load(sprintf('%sREFREPDATA.RData',data_path))
# load(sprintf('%sVP2PIEDATA.RData',data_path))
# load(sprintf('%sVP2RNDDATA.RData',data_path))

t <- dat.cnt2$Beneficial[1]
sen.fdr <- ggplot(dat.cnt2) +
  geom_area(aes(x = (cen-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', lwd = 0.5, linetype = 'dashed') +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100) +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#303F9F',
                                'Neutral'='#9E9E9E',
                                'Beneficial'='#FFA000'),
                     guide = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#303F9F',
                               'Neutral'='#9E9E9E',
                               'Beneficial'='#FFA000')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(5, "mm"),
        legend.position = 'right',
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))

spe.1 <- ggplot(stats.tmp[stats.tmp$hours == stats.tmp$cont_hrs,],
                aes(x = p, y = round((1-fpr),4)*100, col = as.factor(hours))) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), aes(group = 1),
               geom="ribbon", color="#212121", fill = "white",alpha = 0.4) +
  stat_summary(fun=mean, geom="line", color="#E040FB", lwd =1) +
  # stat_summary(fun=mean, geom="point", color="#FFC107", size = 2) +
  labs(x = "p-value",
       y = "Specificity") +
  scale_x_continuous(breaks = seq(-1,1,0.025),
                     minor_breaks = seq(-1,1,0.0125)) +
  scale_y_continuous(breaks = seq(0,200,1),
                     minor_breaks = seq(0,200,0.5),
                     labels = paste(sprintf('%d',seq(0,200,1)),'%', sep = '')) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(90, 100)) +
  scale_color_discrete() +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(5, "mm")) #+
# coord_cartesian(xlim = c(0, 0.1), ylim = c(95, 100))

refrep$power[refrep$hours < refrep$cont_hrs] <- refrep$Deleterious[refrep$hours < refrep$cont_hrs]/910 * 100
refrep$power[refrep$hours > refrep$cont_hrs] <- refrep$Beneficial[refrep$hours > refrep$cont_hrs]/910 * 100
refrep$abs_cen <- abs(1-refrep$cen) * 100
refrep$rep <- as.factor(refrep$rep)
refrep$ref_prop <- as.factor(refrep$ref_prop)

sen.rep <- ggplot(refrep[round(refrep$abs_cen) <= 5,],
                  aes(x = rep, y = power, fill = ref_prop)) +
  geom_boxplot() +
  # geom_smooth(aes(x = rep, y = power), method = 'loess', se = F, lwd = 1.2) +
  labs(x = 'No. of Replicates',
       y = 'Sensitivity') +
  scale_x_discrete(breaks = seq(0,16,2)) +
  scale_y_continuous(breaks = seq(0,100,10),
                     minor_breaks = seq(0,105,5),
                     labels = paste(seq(0,100,10),'%',sep='')) +
  scale_fill_manual(name = 'Reference\nProportion',
                    breaks = c(0.0625,0.125,0.1875,0.25),
                    values = c('#C5CAE9','#448AFF','#3F51B5','#212121'),
                    labels = c(sprintf('%0.2f%%',c(0.0625,0.125,0.1875,0.25)*100))) +
  coord_cartesian(ylim = c(0,100)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

fig.lid.vp1 <- ggarrange(spe.1, sen.fdr, sen.rep,
                         nrow = 1, ncol = 3,
                         labels = c('a','b','c'),
                         label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                           hjust=-1))

ggsave(sprintf("%sFIG_VP1.jpg",out_path), fig.lid.vp1,
       height = 55, width = two.c, units = 'mm',
       dpi = 300)


load(sprintf('%sBEANSPEDATA.RData',data_path))

t <- dat.cnt2$Beneficial[1]
sen.fdr.bean <- ggplot(dat.cnt2) +
  geom_area(aes(x = (cen-1)*100, y = Neutral, fill = 'Neutral'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Beneficial, fill = 'Beneficial'), alpha = 1) +
  geom_area(aes(x = (cen-1)*100, y = Deleterious, fill = 'Deleterious'), alpha = 1) +
  geom_vline(xintercept = seq(-2,2,0.025)*100, col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_hline(yintercept = c(t * seq(0,1,0.05)), col = '#757575', lwd = 0.5, alpha =0.5) +
  geom_vline(xintercept = c(-5,5), col = 'red', lwd = 0.5, linetype = 'dashed') +
  labs(x = 'Fitness Effect',
       y = 'Sensitivity') +
  scale_y_continuous(breaks = c(t * seq(0,1,0.1)),
                     minor_breaks = c(t * seq(0,1,0.05)),
                     labels = c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_x_continuous(breaks = seq(-2,2,0.05)*100,
                     minor_breaks = seq(-2,2,0.025)*100) +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious'),
                     values = c('Deleterious'='#303F9F',
                                'Neutral'='#9E9E9E',
                                'Beneficial'='#FFA000'),
                     guide = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#303F9F',
                               'Neutral'='#9E9E9E',
                               'Beneficial'='#FFA000')) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(5, "mm"),
        legend.position = 'right',
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
  coord_cartesian(xlim = c(-10,10),
                  ylim = c(0,t))

fig.lidbean.vp1 <- ggarrange(spe.1, sen.fdr, sen.rep, sen.fdr.bean,
                         nrow = 2, ncol = 2,
                         labels = c('a','b','c','d'),
                         label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                           hjust=-1))

ggsave(sprintf("%sFIG_withBean_VP1.jpg",out_path), fig.lidbean.vp1,
       height = 160, width = two.c, units = 'mm',
       dpi = 300)

##### VP2 RESULTS
load("output/4C4_rnd_fit_data.RData")
load("output/4C4_rnd_ref_data.RData")
load("output/4C4_rnd_sta_data.RData")

sta.data$res <- sta.data$effect_p
sta.data$res[sta.data$effect_p == 'Beneficial' & sta.data$from == 'more'] <- 'True Beneficial'
sta.data$res[sta.data$effect_p == 'Beneficial' & sta.data$from == 'less'] <- 'False Positive Switch'
sta.data$res[sta.data$effect_p == 'Deleterious' & sta.data$from == 'more'] <- 'False Positive Switch'
sta.data$res[sta.data$effect_p == 'Deleterious' & sta.data$from == 'less'] <- 'True Deleterious'
sta.data$res[sta.data$effect_p == 'Neutral' & sta.data$from == 'more'] <- 'False Negative'
sta.data$res[sta.data$effect_p == 'Neutral' & sta.data$from == 'less'] <- 'False Negative'

sta.data <- sta.data[!sta.data$hours == 0,]
sta.data$method[sta.data$method == "BEAN"] <- "MCAT"

pie.dat <- data.frame(rbind(cbind(hours = sta.data$hours[sta.data$method == 'MCAT'],
                                  value = sta.data$from[sta.data$method == 'MCAT'],
                                  method = 'TRUTH'),
                            cbind(hours = sta.data$hours[sta.data$method == 'MCAT'],
                                  value = sta.data$res[sta.data$method == 'MCAT'],
                                  method = 'MCAT'),
                            cbind(hours = sta.data$hours[sta.data$method == 'LID'],
                                  value = sta.data$res[sta.data$method == 'LID'],
                                  method = 'LID')))

pie.dat$value[pie.dat$value == 'less'] <- 'True Deleterious'
pie.dat$value[pie.dat$value == 'more'] <- 'True Beneficial'

pie.dat <- arrange(transform(pie.dat,
                             method=factor(method,levels=c('LID','TRUTH','MCAT'))),
                   method)

pie.per <- plyr::count(pie.dat, vars = c('method','value'))
pie.per$per <- pie.per$freq/(sum(pie.per$freq)/3) * 100
# rnd_data$hours <- factor(rnd_data$hours, levels = c(sort(unique(rnd_data$hours))))
rnd_data$virtual <- sprintf('Ref. at %0.2f hr', rnd_data$hours)
rnd_data$virtual <- factor(rnd_data$virtual, levels = sprintf('Ref. at %0.2f hr', c(unique(sort(rnd_data$hours)))))

pie.per$value <- factor(pie.per$value, levels = c('True Deleterious','False Negative','False Positive Switch','True Beneficial'))

fig.lid.vp2 <- ggplot(data = pie.per, aes(x = "", y = per, fill = value)) +
  geom_bar(stat = 'identity', col = 'black') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",per)),
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('True Deleterious',
                               'False Negative',
                               'False Positive Switch',
                               'True Beneficial'),
                    values = c('False Negative'='#9E9E9E',
                               'True Deleterious'='#303F9F',
                               'False Positive Switch'='#F5F5F5',
                               'True Beneficial'='#FFA000')) +
  
  facet_wrap(.~method) +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = txt),
        legend.box.spacing = unit(0.5,"mm"))

ggsave(sprintf("%sFIG_VP2.jpg",out_path), fig.lid.vp2,
       height = 60, width = two.c, units = 'mm',
       dpi = 300)

##### OES Pilots
lm_eqn <- function(df){
  m <- lm(bg ~ average, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

fdata <- dbGetQuery(conn, 'select e.*, f.nt_len from
                    (select c.*, d.pg_2012 from
                    (select a.strain_id, a.orf_name,
                    a.cs_mean oesp1_cs_mean, a.cs_median oesp1_cs_median,
                    c.p oesp1_p, c.stat oesp1_stat, c.es oesp1_es,
                    b.cs_mean oesp2_cs_mean, b.cs_median oesp2_cs_median,
                    d.p oesp2_p, d.stat oesp2_stat, d.es oesp2_es
                    from OESP1_FS_6144_FITNESS_STATS a, OESP2_FS_6144_FITNESS_STATS b,
                    OESP1_FS_6144_PVALUE c, OESP2_FS_6144_PVALUE d
                    where a.strain_id = b.strain_id and a.hours = 22 and b.hours = 24
                    and c.hours = a.hours and d.hours = b.hours and a.strain_id = c.strain_id and a.strain_id = d.strain_id) c
                    left join
                    (select * from PROTOGENES) d
                    on c.orf_name = d.orf_name) e
                    left join
                    (select orf_name, length(seq_nt_atg) nt_len from ORFS_SEQUENCES) f
                    on e.orf_name = f.orf_name')

fdata$pg_2012[is.na(fdata$pg_2012)] <- 0
fdata$pg_2012[str_detect(fdata$orf_name, 'smorf')] <- 2
fdata$pg_2012 <- as.factor(fdata$pg_2012)
sum(fdata$pg_2012 == 2)

fdata <- fdata[(abs(fdata$oesp1_cs_mean-fdata$oesp2_cs_mean)/rowMeans(data.frame(fdata$oesp1_cs_mean,fdata$oesp2_cs_mean)))*100 < 50,]

fdata$oesp1_pheno[fdata$oesp1_p <= 0.05 & fdata$oesp1_stat > 0] <- 'Beneficial'
fdata$oesp1_pheno[fdata$oesp1_p <= 0.05 & fdata$oesp1_stat < 0] <- 'Deleterious'
fdata$oesp1_pheno[is.na(fdata$oesp1_pheno)] <- 'Neutral'

fdata$oesp2_pheno[fdata$oesp2_p <= 0.05 & fdata$oesp2_stat > 0] <- 'Beneficial'
fdata$oesp2_pheno[fdata$oesp2_p <= 0.05 & fdata$oesp2_stat < 0] <- 'Deleterious'
fdata$oesp2_pheno[is.na(fdata$oesp2_pheno)] <- 'Neutral'

pcc <- cor(fdata$oesp1_cs_mean, fdata$oesp2_cs_mean, method = 'pearson')

pilot.fit <- ggplot(fdata,
       aes(x = oesp1_cs_mean, y = oesp2_cs_mean)) +
  geom_point(aes(col = pg_2012), alpha = 0.7) +
  geom_abline(col = 'red', linetype = 'dashed') +
  geom_smooth(method = 'lm', se = F, col = 'black') +
  # geom_density_2d(col = 'black', lwd = 0.1) +
  annotate("text", x = 0.3, y = 0.1,
           label = lm_eqn(data.frame(average = fdata$oesp1_cs_mean, bg = fdata$oesp2_cs_mean)),
           size = 2, parse = T,
           hjust = 0) +
  annotate("text",x = 0.3, y = 1,
           label = sprintf('Pearson correlation = %0.2f',pcc),
           size = 2) +
  labs(x = 'Pilot #1 Mean Fitness',
       y = 'Pilot #2 Mean Fitness') +
  scale_color_manual(name = 'ORF Type',
                       breaks = c('0', '1','2'),
                       values = c('0' = 'red', '1' = 'blue', '2' = 'green'),
                       labels = c('0' = 'Gene', '1' = 'Annotated\nProto-gene', '2' = 'Unannotated\nProto-gene')) +
  coord_cartesian(xlim = c(0, 1.1),
                ylim = c(0,1.1)) +
  scale_x_continuous(breaks = seq(-2,2,0.2),
                     minor_breaks = seq(-2,2,0.05)) +
  scale_y_continuous(breaks = seq(-2,2,0.2),
                     minor_breaks = seq(-2,2,0.05)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'right',
        strip.text = element_text(size = txt),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0, unit = 'pt'))
# ggsave(sprintf("%sFIG_P1vsP2_fitness.jpg",out_path), pilot.fit,
#        height = 75, width = one.c, units = 'mm',
#        dpi = 300)


fdata$phenotype <- paste(strtrim(fdata$oesp1_pheno,3), strtrim(fdata$oesp2_pheno,3), sep = '/')
fdata$es <- rowMeans(cbind(fdata$oesp1_es, fdata$oesp2_es), na.rm = T)

pie.per <- plyr::count(fdata, vars = c('phenotype'))
pie.per$per <- pie.per$freq/sum(pie.per$freq) * 100

(805+231)/(805+231+155+129)

pilot.pheno <- ggplot(data = pie.per, aes(x = "", y = per, fill = phenotype)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%\n(%d)",per,freq)),
                   position = position_stack(vjust = 0.5),
                   label.size = 0.15,
                   size = 2,
                   # box.padding = 0.1,
                   seed = 100,
                   show.legend = F) +
  scale_fill_discrete(name = 'Phenotype\n(P1/P2)') +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles, color = 'white'),
        axis.text = element_blank(),
        axis.ticks = element_line(color = 'white'),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'right',
        strip.text = element_text(size = txt),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0, unit = 'pt'))

# ggsave(sprintf("%sOESP_EFFECTS.jpg",out_path),
#        height = one.c, width = one.c, units = 'mm',
#        dpi = 300)

pilot.phenoes <- ggplot(fdata,
       aes(x = oesp1_es, y = oesp2_es)) +
  geom_point(aes(col = phenotype)) +
  # geom_density_2d(col = 'black', lwd = 0.1) +
  labs(x = 'Pilot #1 Effect Size',
       y = 'Pilot #2 Effect Size') +
  scale_color_discrete(name = 'Phenotype\n(P1/P2)') +
  coord_cartesian(xlim = c(-1, 0.1),
                  ylim = c(-1,0.1)) +
  scale_x_continuous(breaks = seq(-2,2,0.2),
                     minor_breaks = seq(-2,2,0.05)) +
  scale_y_continuous(breaks = seq(-2,2,0.2),
                     minor_breaks = seq(-2,2,0.05)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'right',
        strip.text = element_text(size = txt),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0, unit = 'pt'))


fig.pilot <- ggarrange(pilot.fit, pilot.pheno,
                         nrow = 1, ncol = 2,
                         labels = c('a','b'),
                         label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                           hjust=-1))

# ggsave(sprintf("%sFIG_PILOT.jpg",out_path), fig.pilot,
#        height = 80, width = two.c, units = 'mm',
#        dpi = 300)


cnt.data <- fdata[c('pg_2012','oesp1_pheno','oesp2_pheno')]
cnt.data <- melt(cnt.data, id.vars = 'pg_2012')
cnt.data <- plyr::count(cnt.data, vars = c('pg_2012','variable','value'))

colnames(cnt.data) <- c('orf_type','pilot','phenotype','counts')

for (pil in unique(cnt.data$pilot)) {
  for (orf in unique(cnt.data$orf_type)) {
    cnt.data$total[cnt.data$pilot == pil & cnt.data$orf_type == orf] <- sum(cnt.data$counts[cnt.data$pilot == pil & cnt.data$orf_type == orf])
  }
}

cnt.data$minus <- cnt.data$total - cnt.data$counts

or <- NULL
i <- 1
for (pil in unique(cnt.data$pilot)) {
  for (orf in unique(cnt.data$orf_type)) {
    if (orf > 0) {
      for (pheno in unique(cnt.data$phenotype)) {
        or$pilot[i] <- pil
        or$orf[i] <- orf
        or$phenotype[i] <- pheno
        hello <- fisher.test(matrix(t(c(cnt.data$counts[cnt.data$pilot == pil & cnt.data$orf_type == orf & cnt.data$phenotype == pheno],
                                        cnt.data$minus[cnt.data$pilot == pil & cnt.data$orf_type == orf & cnt.data$phenotype == pheno],
                                        cnt.data$counts[cnt.data$pilot == pil & cnt.data$orf_type == 0 & cnt.data$phenotype == pheno],
                                        cnt.data$minus[cnt.data$pilot == pil & cnt.data$orf_type == 0 & cnt.data$phenotype == pheno])),
                                    nrow = 2))
        or$or[i] <- hello$estimate
        or$bottom[i] <- hello$conf.int[1]
        or$top[i] <- hello$conf.int[2]
        or$pvalue[i] <- hello$p.value
        i <- i +1
      }
    }
  }
}

or <- data.frame(or)
or$or <- as.double(as.character(or$or))
or$pilot <- as.character(or$pilot)
or$pilot[or$pilot == 'oesp1_pheno'] <- 'Pilot #1'
or$pilot[or$pilot == 'oesp2_pheno'] <- 'Pilot #2'
or$orf <- as.character(or$orf)
or$orf[or$orf == '1'] <- 'Annotated\nProto-gene'
or$orf[or$orf == '2'] <- 'Unannotated\nProto-gene'

pilot.or <- ggplot(or,aes(x=orf,y=or,ymin=bottom,ymax=top)) +
  geom_point(stat="identity", shape=21, size=2, stroke=1, fill = "white") +
  geom_errorbar(width=0.25) +
  facet_wrap (~pilot*phenotype , nrow=2) +
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  labs(x = 'ORF Type',
       y = 'Odds Ratio (log)') +
  # scale_x_discrete(level = c('Annotated\nProto-gene', 'Unannotated\nProto-gene'))
  theme_linedraw() +
  # theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),
  #       axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),
  #       legend.position="none") + 
  scale_y_continuous(trans="log10", breaks=c(0.2,0.5,1,2,5,10)) +
  annotation_logticks(sides="l") +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.15,0,0.15,0, "mm")),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0, unit = 'pt'))

fig.pilot2 <- ggarrange(pilot.fit, pilot.pheno, pilot.or,
                        nrow = 1, ncol = 3,
                        labels = c('a','b','c'),
                        label.args = list(gp=gpar(font = 2, fontsize = lbls),
                                          hjust=-1))

ggsave(sprintf("%sFIG_PILOT2.jpg",out_path), fig.pilot2,
       height = 60, width = two.c, units = 'mm',
       dpi = 300)

##### NEW PROTO_GENES WITH EACH CONDITION
sd.ben.pg <- dbGetQuery(conn, 'select exp_id, orf_name
                        from brian_031918.DATASET_6
                        where exp_id in (28,31,33,91,93)
                        and effect_cs = "beneficial"
                        and orf_name in
                        (select orf_name from PROTOGENES
                         where pg_2012 = 1)
                        order by exp_id')

cnt.unique <- NULL
temp <- NULL
for (exp in unique(sd.ben.pg$exp_id)) {
  cnt.unique <- rbind(cnt.unique, c(exp,sum(sd.ben.pg$orf_name[sd.ben.pg$exp_id == exp] %in% sd.ben.pg$orf_name[sd.ben.pg$exp_id == exp]),
                                sum(sd.ben.pg$orf_name[sd.ben.pg$exp_id == exp] %in% temp),
                                length(unique(c(temp, sd.ben.pg$orf_name[sd.ben.pg$exp_id == exp])))))
  temp <- unique(c(temp, sd.ben.pg$orf_name[sd.ben.pg$exp_id == exp]))
}
cnt.unique <- data.frame(cnt.unique)
colnames(cnt.unique) <- c('exp_id','exp_ben','overlap','total_ben')
cnt.unique$unique_ben <- cnt.unique$exp_ben - cnt.unique$overlap
sum(cnt.unique$unique_ben)

cnt.unique$condition[cnt.unique$exp_id == 28] <- 'N+,C+'
cnt.unique$condition[cnt.unique$exp_id == 31] <- 'N+,C++'
cnt.unique$condition[cnt.unique$exp_id == 91] <- 'N++,C+'
cnt.unique$condition[cnt.unique$exp_id == 33] <- 'N++,C++'
cnt.unique$condition[cnt.unique$exp_id == 93] <- 'N-,C++'

cum.pgs <- ggplot(cnt.unique,
       aes(x = as.factor(exp_id), y = total_ben)) +
  geom_line(group = 1, col = 'Red') +
  geom_point(size = 3) +
  scale_x_discrete(labels=cnt.unique$condition) +
  scale_y_continuous(breaks = seq(0,30,5)) +
  labs(x = 'Conditions Tested (Vakirlis et al., Nature 2020)',
       y = 'Cummulative Adaptive Proto-genes') +
  theme_linedraw() +
  coord_cartesian(ylim = c(0, 30)) +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.15,0,0.15,0, "mm")),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0, unit = 'pt'))

ggsave(sprintf("%sFIG_CUMULATIVE_PGS.jpg",out_path), cum.pgs,
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)


##### HEATMAP FROM SANDIEGO DATA
cond5 <- dbGetQuery(conn, 'select distinct a.orf_name, hours, n, colony_size normalized_cs, q_cs, effect_cs,exp_id,
                    AGE, chromosome, translation_media, category
                    from brian_031918.DATASET_6 a, ANNE.SUMMARY_2012 b where a.orf_name = b.orf_name and exp_id in (28,31,91,33, 93) and a.orf_name !="BF_control"')
cond5$condition[cond5$exp_id == 28] <- 'N+/C+'
cond5$condition[cond5$exp_id == 31] <- 'N+/C++'
cond5$condition[cond5$exp_id == 91] <- 'N++/C+'
cond5$condition[cond5$exp_id == 33] <- 'N++/C++'
cond5$condition[cond5$exp_id == 93] <- 'N-/C++'

# cond5.hm <- with(cond5, tapply(effect_cs, list(orf_name, exp_id), c))
# cond5.hm[cond5.hm == 'beneficial'] <- 1
# cond5.hm[cond5.hm == 'neutral'] <- 0
# cond5.hm[cond5.hm == 'deleterious'] <- -1
# cond5.hm <- data.frame(cond5.hm, stringsAsFactors = F)

cond5.hm <- with(cond5, tapply(normalized_cs, list(orf_name, exp_id), c))
cond5.hm[is.na(cond5.hm)] <- 0

cond5.hm <- data.frame(cond5.hm)
cond5.hm <- cond5.hm[cond5.hm$X28 > 0.01 &
           cond5.hm$X31 > 0.01 &
           cond5.hm$X33 > 0.01 &
           cond5.hm$X91 > 0.01 &
           cond5.hm$X93 > 0.01 ,]
cond5.hm[cond5.hm == 0,]

orf_type <- NULL
age <- NULL
for (orf in rownames(cond5.hm)) {
  orf_type <- c(orf_type, cond5$category[cond5$orf_name == orf][1])
  age <- c(age, cond5$AGE[cond5$orf_name == orf][1])
}
cond5.hm$orf_type <- orf_type
cond5.hm$age <- age

annot_df <- data.frame(orf_type = cond5.hm$orf_type,
                       age = cond5.hm$age)
col = list(orf_type = c("gene" = "green", "proto-gene" = "black"),
           age = circlize::colorRamp2(c(0, 10),
                                      c("lightblue", "purple")))
ha <- HeatmapAnnotation(df = annot_df, col = col,
                        which = 'row')

Heatmap(cond5.hm[1:5], name = "fitness",
        left_annotation  = ha,
        show_row_names = F)

##### PLATEMAPS
p2c.dat <- dbGetQuery(conn, 'select * from 4C4_pos2coor a, 4C4_pos2orf_name b
                      where a.pos = b.pos')
p2c.dat$colony[is.na(p2c.dat$orf_name)] <- 'Gap'
p2c.dat$colony[p2c.dat$orf_name == 'BF_control'] <- 'Reference'
p2c.dat$colony[p2c.dat$orf_name != 'BF_control'] <- 'Mutant'

stocks <- ggplot(p2c.dat[p2c.dat$density == 384,],
                 # aes(x = col, y = row, col = orf_name)) +
                 aes(x = col, y = row, col = colony)) +
  geom_point(size = 1, shape = 15) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant","Gap"),
                     values = c("Reference" = "#3F51B5",
                                "Mutant" = "#FFC107",
                                "Gap" = "#D32F2F"),
                     guide = F) +
  labs(x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 2) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")))

wctp1 <- ggplot(p2c.dat[p2c.dat$density == 384,],
                # aes(x = col, y = row, col = orf_name)) +
                aes(x = col, y = row, col = colony)) +
  geom_point(size = 1, shape = 16) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant","Gap"),
                     values = c("Reference" = "#3F51B5",
                                "Mutant" = "#FFC107",
                                "Gap" = "#D32F2F"),
                     guide = F) +
  labs(x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 2) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")))

up1tp2 <- ggplot(p2c.dat[p2c.dat$density == 1536,],
                 # aes(x = col, y = row, col = orf_name)) +
                 aes(x = col, y = row, col = colony)) +
  geom_point(size = 0.5, shape = 16) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant","Gap"),
                     values = c("Reference" = "#3F51B5",
                                "Mutant" = "#FFC107",
                                "Gap" = "#D32F2F")) +
  labs(x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 2) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title.x = element_text(size = titles),
        axis.title.y = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm"))) +
  guides(color = guide_legend(override.aes = list(size = 3)))

up2 <- ggplot(p2c.dat[p2c.dat$density == 6144,],
              # aes(x = col, y = row, col = orf_name)) +
              aes(x = col, y = row, col = colony)) +
  geom_point(size = 0.2, shape = 16) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Reference","Mutant","Gap"),
                     values = c("Reference" = "#3F51B5",
                                "Mutant" = "#FFC107",
                                "Gap" = "#D32F2F"),
                     guide = F) +
  labs(x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 2) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title.x = element_text(size = titles),
        axis.title.y = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")))

# plt.maps <- ggarrange(stocks, wctp1, up1tp2, up2,
#                       ncol = 1, nrow = 4)
plt.maps <- ggarrange(wctp1, up1tp2, up2,
                      ncol = 3, nrow = 1)

ggsave(sprintf("%sPLATEMAP.jpg",out_path), plt.maps,
       height = one.c, width = two.c*1.5, units = 'mm',
       dpi = 300)

#### PLATEMAP - REPLICATE
p2c.dat <- dbGetQuery(conn, 'select * from 4C4_pos2coor a, 4C4_pos2orf_name b
                      where a.pos = b.pos')
p2c.dat$colony[is.na(p2c.dat$orf_name)] <- 'Others'
p2c.dat$colony[p2c.dat$orf_name == 'BFC424'] <- 'Example'
p2c.dat$colony[p2c.dat$orf_name != 'BFC424'] <- 'Others'

wctp1 <- ggplot(p2c.dat[p2c.dat$density == 384,],
                # aes(x = col, y = row, col = orf_name)) +
                aes(x = col, y = row, col = colony)) +
  geom_point(size = 1, shape = 16) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Example","Others"),
                     values = c("Example" = "#E91E63",
                                "Others" = "#BDBDBD"),
                     guide = F) +
  labs(x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 2) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")))

up1tp2 <- ggplot(p2c.dat[p2c.dat$density == 1536,],
                 # aes(x = col, y = row, col = orf_name)) +
                 aes(x = col, y = row, col = colony)) +
  geom_point(size = 0.5, shape = 16) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Example","Others"),
                     values = c("Example" = "#E91E63",
                                "Others" = "#BDBDBD")) +
  labs(x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 2) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title.x = element_text(size = titles),
        axis.title.y = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm"))) +
  guides(color = guide_legend(override.aes = list(size = 3)))

up2 <- ggplot(p2c.dat[p2c.dat$density == 6144,],
              # aes(x = col, y = row, col = orf_name)) +
              aes(x = col, y = row, col = colony)) +
  geom_point(size = 0.2, shape = 16) +
  # scale_color_discrete(guide = F) +
  scale_y_reverse() +
  scale_color_manual(name = 'Colony Type',
                     breaks = c("Example","Others"),
                     values = c("Example" = "#E91E63",
                                "Others" = "#BDBDBD"),
                     guide = F) +
  labs(x = 'Column', y = "Row") +
  facet_wrap(.~density*plate,
             scales = 'free',
             nrow = 2) +
  theme_linedraw() +
  theme(title = element_text(size = titles),
        axis.title.x = element_text(size = titles),
        axis.title.y = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        strip.text = element_text(size = txt,
                                  margin = margin(0.5,0,0.5,0, "mm")))

plt.map.rep <- ggarrange(wctp1, up1tp2, up2,
                      ncol = 3, nrow = 1)

ggsave(sprintf("%sPLATEMAP_REP.jpg",out_path), plt.map.rep,
       height = one.c, width = two.c*1.5, units = 'mm',
       dpi = 300)



