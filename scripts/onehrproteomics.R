library(RMariaDB)
library(readxl)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(grid)
library(tidyverse)

hello <- data.frame(read_xlsx('rawdata/mcp.M113.034769-3.xlsx', sheet = 1))

ggplot(hello) +
  geom_histogram(aes(x = X..AAs), binwidth = 25) +
  geom_vline(xintercept = 100, col = 'red', linetype = 'dashed') +
  labs(x = 'AA Length',
       y = 'Frequency') +
  theme_linedraw()
