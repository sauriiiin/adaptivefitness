##### 4C_ALLDAT
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/23/2019

##### ALL Data from the virtual plates and normal plates taken together

library(readxl)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(egg)

##### INITIALIZE
expt_name = '4C3_GA1'
expt = 'FS1-1'
out_path = 'figs/'
density = 6144

##### LOAD DATA
