# LIBRARIES
library(ggplot2)
library(plyr)
library(ggpubr)
library(effsize)
library(nnet)
library(reshape2)
library(dplyr)
library(tidyr)

out_path = 'figs/proposal/';

###
orf_table<-read.csv(paste(out_path, "Data1.csv",sep=""),header=TRUE)
orf_table $protogene <-factor(orf_table $emergence_status, levels=c("Emerging ORFs","Established ORFs"))
orf_table $overexpression_relative_fitness <-factor(orf_table $overexpression_relative_fitness, levels=c("increased","decreased","unchanged"))
fitness_table<-read.csv(paste(out_path, "Data3.csv",sep=""),header=TRUE)

###
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$emergence_status =="Emerging ORFs",c("orf_name","no_helices_TMHMM")]
df$is_tm<-ifelse(df$no_helices_TMHMM>0,1,0)
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","emergence_status")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df$fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate proportions with at least 1 TM domain
stat_df<-df[,c("fitness_category","is_tm")]
domains<-NULL
domains $total<- as.vector(table(stat_df $fitness_category))
domains $yes<-as.vector(sapply(levels(stat_df $fitness_category), function(x) sum(stat_df[stat_df $fitness_category ==x,]$is_tm)))
domains $fraction<-domains $yes/domains $total
domains $sder<-sqrt(domains $fraction * (1-domains $fraction) / domains $total)
domains <- data.frame(domains)
domains $fitness<-c("Adaptive","Neutral","Deleterious")
domains $fitness<-factor(domains $fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(domains,aes(x= fitness,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge(), alpha=0.5)+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="With TM domain (TMHMMM)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=8), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
image

#
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","emergence_status","no_helices_TMHMM")], all.x=TRUE)
adapt_pg_tm <- fitness_df[fitness_df$emergence_status == 'Emerging ORFs' &
             fitness_df$effect_cs=="increased" &
             fitness_df$no_helices_TMHMM > 0, c("orf_name", "exp_environment", "exp_id")]
count(adapt_pg_tm, vars = 'exp_environment')

adapt_pg_tm$orf_name <- as.character(adapt_pg_tm$orf_name)
adapt_pg_tm$exp_environment <- as.character(adapt_pg_tm$exp_environment)
cnt.unique2 <- NULL
temp <- NULL
for (exp in unique(sort(adapt_pg_tm$exp_id))) {
  cnt.unique2 <- rbind(cnt.unique2, c(exp,unique(adapt_pg_tm$exp_environment[adapt_pg_tm$exp_id == exp]),
                                    sum(adapt_pg_tm$orf_name[adapt_pg_tm$exp_id == exp] %in% adapt_pg_tm$orf_name[adapt_pg_tm$exp_id == exp]),
                                    sum(adapt_pg_tm$orf_name[adapt_pg_tm$exp_id == exp] %in% temp),
                                    length(unique(c(temp, adapt_pg_tm$orf_name[adapt_pg_tm$exp_id == exp])))))
  temp <- unique(c(temp, adapt_pg_tm$orf_name[adapt_pg_tm$exp_id == exp]))
}
cnt.unique2 <- data.frame(cnt.unique2)
colnames(cnt.unique2) <- c('exp_id','exp_environment','exp_ben','overlap','total_ben')
cnt.unique2$exp_id <- as.numeric(as.character(cnt.unique2$exp_id))
cnt.unique2$total_ben <- as.numeric(as.character(cnt.unique2$total_ben))

cum.tm.pgs <- ggplot(cnt.unique2,
                  aes(x = as.factor(exp_id), y = total_ben,
                      col = 'with TM')) +
  geom_line(group = 1) +
  geom_point(size = 3, col = 'black') +
  scale_x_discrete(labels=cnt.unique2$exp_environment) +
  scale_y_continuous(breaks = seq(0,30,5)) +
  labs(x = 'Conditions Tested (Vakirlis et al., Nature 2020)',
       y = 'Cummulative Adaptive Proto-genes') +
  theme_linedraw() +
  coord_cartesian(ylim = c(0, 30)) +
  scale_color_discrete(name = '') +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.15,0,0.15,0, "mm")),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0, unit = 'pt'))

cum.adapt <- cum.tm.pgs + geom_line(data = cnt.unique, aes(x = as.factor(exp_id), y = total_ben, col = 'all'),
                       group = 1) + geom_point(data = cnt.unique, size = 3, col = 'black')

ggsave(sprintf("%sFIG_CUMULATIVE_PGS_ADAPT.jpg",out_path), cum.adapt,
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)

###

pg_adapt_cnt <- count(adapt_pg_tm, vars = 'orf_name')
sum(pg_adapt_cnt$freq[pg_adapt_cnt$freq > 3])
pg_adapt_cnt[pg_adapt_cnt$freq >= 3,]
paired_events <- NULL
i <- 1
for (freq in sort(unique(pg_adapt_cnt$freq))) {
  paired_events$freq_cutoff[i] <- freq
  paired_events$pg_cnt[i] <- length(pg_adapt_cnt$orf_name[pg_adapt_cnt$freq >= freq])
  paired_events$pairs[i] <- sum(pg_adapt_cnt$freq[pg_adapt_cnt$freq >= freq])
  i <- i + 1
}
paired_events <- data.frame(paired_events)

adapt.pgs.pairs <- ggplot(paired_events,
       aes(x = freq_cutoff, y = pg_cnt)) +
  geom_line(aes(col = 'AdPTs'), group = 1) +
  geom_point() +
  geom_line(data = paired_events, 
            aes(x = freq_cutoff, y = pairs, col = 'AdPTs-Cond Pairs'), group = 1) +
  geom_point(data = paired_events, 
             aes(x = freq_cutoff, y = pairs), group = 1) +
  theme_linedraw() +
  labs(x = 'Times Adaptive',
       y = 'Frequency') +
  scale_color_discrete(name = '') +
  scale_x_continuous(labels = paste('>=',seq(1,5),sep='')) +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.15,0,0.15,0, "mm")),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0, unit = 'pt'))
  
adapt.conds <- ggplot(adapt_pg_tm) +
  geom_histogram(aes(x = exp_environment), stat = 'count') +
  labs(x = 'Environment',
       y = 'No. of AdPTs') +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'right',
        strip.text = element_text(size = txt,
                                  margin = margin(0.15,0,0.15,0, "mm")),
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0, unit = 'pt'))

adapt.counts <- ggarrange(adapt.conds, adapt.pgs.pairs,
          nrow = 1)
ggsave(sprintf("%sFIG_ADAPT_COUNTS.jpg",out_path), adapt.counts,
       height = one.c, width = two.c, units = 'mm',
       dpi = 300)

ggplot(paired_events) +
  geom_point(aes(x = pg_cnt, y = pairs))
  
table(adapt_pg_tm)


###
fitness_df
ben.prop.all <- count(fitness_df, vars = c('exp_id','effect_cs'))
ben.prop.all$prop <- ben.prop.all$freq/4647 * 100
ben.prop.all
