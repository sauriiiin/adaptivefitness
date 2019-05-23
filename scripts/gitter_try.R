##### GITTER
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 05/23/2019
##### TESTING OUT GITTER OUTPUT

##### INITIALIZE
# source("https://bioconductor.org/biocLite.R")
# biocLite("EBImage")
# install.packages("tiff")
library(gitter)

jpeg.dir <- "/home/sbp29/RAW_Data/4Control/4C3/Ga/S1/"
jpeg.files <- list.files(path = jpeg.dir,
                         pattern = ".JPG", recursive = TRUE)
junk <- list.files(path = jpeg.dir,
           pattern = ".JPG.", recursive = TRUE)
jpeg.files <- jpeg.files[!jpeg.files %in% junk]

##### RUNNING GITTER ON JPGs
gitter.batch(sprintf("%s%s",jpeg.dir,jpeg.files), plate.format = c(64,96),
             remove.noise = F,
             dat.save = "/home/sbp29/RAW_Data/4Control/4C3/Ga/S1_gitter/")

gitter(sprintf("%s%s",jpeg.dir,jpeg.files)[20], plate.format = c(64,96),
             remove.noise = F,
             dat.save = "/home/sbp29/RAW_Data/4Control/4C3/Ga/S1_gitter/",
      .is.ref = T, .params =  attr(31,'window'))

##### NO INBUILT NORMALIZATION IN GITTER
##### SGA Tools utilizes gitter for colony measurements
#####   normalization --> SGA score relies on an SGA type experimental design
