hello <- read.csv(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/costanzo2021_singlemutant_differentialfitness.csv')
pgs <- read.csv(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/protogenes.csv')
head(hello)


hello$orf_type[hello$SystematicName %in% pgs$orf_name] <- 'Proto-gene'
hello$orf_type[is.na(hello$orf_type)] <- 'Gene'


hello %>%
  melt(id.vars = c('SystematicName','DiagnosticArray','orf_type'),
       variable.name = 'Condition', value.name = 'DifferentialFitness') %>%
  ggplot(aes(x = DifferentialFitness, y = orf_type)) +
  ggridges::geom_density_ridges(aes(col = orf_type)) +
  scale_color_discrete(guide = F) +
  labs(x = 'Differential Fitness',
       y = '') +
  facet_wrap(.~Condition) +
  theme_linedraw()

hi <- hello %>%
  melt(id.vars = c('SystematicName','DiagnosticArray','orf_type'),
       variable.name = 'Condition', value.name = 'DifferentialFitness') %>%
  data.frame()

for (c in unique(hi$Condition)) {
  temp <- hi[hi$Condition == c,]
  hi$sig[hi$Condition == c] <- hi$DifferentialFitness[hi$Condition == c] <= 
    quantile(hi$DifferentialFitness[hi$Condition == c], c(0.025, 0.975), na.rm = T)[[1]]
  hi$sig[hi$Condition == c] <- hi$DifferentialFitness[hi$Condition == c] >= 
    quantile(hi$DifferentialFitness[hi$Condition == c], c(0.025, 0.975), na.rm = T)[[2]]
}

hi[hi$sig == TRUE & hi$Condition == 'ActinomycinD' & hi$orf_type == 'Proto-gene',]
hi[hi$sig == TRUE & hi$Condition == 'ActinomycinD' & hi$orf_type == 'Proto-gene',]

quantile(hello$ActinomycinD, c(0.025, 0.975), na.rm = T)
dim(hello[hello$orf_type == 'Proto-gene' & hello$ActinomycinD <= -0.2095 & !is.na(hello$ActinomycinD),])


dim(hello[hello$orf_type == 'Proto-gene',])
hello$SystematicName[hello$orf_type == 'Proto-gene' & hello$DiagnosticArray == 'yes']

write.csv(hi, file = '/home/sbp29/R/Projects/adaptivefitness/output/costanzo2021/significant_fitness_effects.csv')


####
hello2 <- read.csv(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/costanzo2021_singlemutant_conditionresponsive.csv')
hello2[hello2$Systematic.Name %in% pgs$orf_name,]

hi %>%
  filter(orf_type == 'Gene')%>%
  group_by(SystematicName, orf_type, sig) %>%
  count() %>%
  data.frame()
