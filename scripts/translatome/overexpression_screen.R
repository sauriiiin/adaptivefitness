
library(stringr)
library(RMariaDB)

source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")


p2c <- dbGetQuery(conn, 'select a.*, b.orf_name
                  from TR_OE_pos2coor a, TR_OE_pos2orf_name b
                  where a.pos = b.pos')


sb <- dbGetQuery(conn, 'select * from lin.TR_OE_ARM1_smudge')

p2c$rep <- as.numeric(str_trunc(as.character(p2c$pos), 8, side = 'left', ellipsis = ''))
sb$pos <- as.numeric(sb$pos)

dbWriteTable(conn,'TR_OE_ONE_THUMB',p2c[p2c$rep %in% sb$pos[sb$density == 1536] & p2c$density == 6144,])

# temp <- dbGetQuery(conn, 'select * from TR_OE_ONE_FS_HO_smudgebox')
# temp1 <- NULL
# temp1$pos <- unique(temp$pos)
# dbWriteTable(conn,'TR_OE_ONE_FS_HO_smudgebox', data.frame(temp1), overwrite = T)
