##### INITIALIZE
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 03/17/2022

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
library(readxl)
library(RMariaDB)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(PoiClaClu)
library(RColorBrewer)
library(effsize)

out_path <- "~/R/Projects/adaptivefitness/output/translatome/"
fig_path <- "~/R/Projects/adaptivefitness/figs/translatome/"

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

`%notin%` <- Negate(`%in%`)

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### ORF CATEGORY
orf_types <- read.csv(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/oe_transient')
orf_types$category[orf_types$is_transient + orf_types$translated + orf_types$is_candidate == 3] <- 'Transient'
orf_types$category[orf_types$is_conserved + orf_types$translated + orf_types$is_candidate == 3] <- 'Conserved'
orf_types$category[is.na(orf_types$category)] <- 'Others'


