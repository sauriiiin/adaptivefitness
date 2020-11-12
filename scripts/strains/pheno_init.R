
strain_pheno <- read.table(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/strains/yeasts_natural.tsv',
                           header = T, sep = "\t")

head(strain_pheno)

ggplot(strain_pheno[strain_pheno$qvalue <= 0.1,]) +
  geom_point(aes(x = strain, y = score, col = condition))

cool_strains <- plyr::count(strain_pheno[strain_pheno$qvalue <= 0.1,], vars = c('strain'))
cool_strains[cool_strains$freq > 10,]

