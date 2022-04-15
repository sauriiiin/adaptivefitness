library(RMariaDB)
library(dplyr)

`%notin%` <- Negate(`%in%`)

source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

orf_types <- read.csv(file = '/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/oe_transient')
orf_types$category[orf_types$is_transient + orf_types$translated + orf_types$is_candidate == 3] <- 'Transient'
orf_types$category[orf_types$is_conserved + orf_types$translated + orf_types$is_candidate == 3] <- 'Conserved'
orf_types$category[is.na(orf_types$category)] <- 'Others'

orf_types %>%
  group_by(category) %>%
  count()

orfs <- dbGetQuery(conn, 'select * from BARFLEX_SPACE_AGAR_180313 union select * from PROTOGENE_COLLECTION')
orfs$`384plate` <- as.numeric(orfs$`384plate`)
orfs <- merge(orfs, orf_types[,c(1,7)], by = 'orf_name', all.x = T)

orfs %>%
  group_by(category) %>%
  count()

orfs %>%
  filter(category == 'Transient') %>%
  group_by( `384plate`) %>%
  count() %>% data.frame()

orfs %>%
  filter(category == 'Transient',
         `384plate` %in% c(2,3,17,23:28)) %>%
  count()

orfs %>%
  filter(category == 'Transient',
         `384plate` %notin% c(2,3,17,23:28)) %>%
  count()

empty384 <- dbGetQuery(conn, 'select * from PLATE384')
empty384$strain_id[empty384$`384row` %in% c(1,16)] <- -2
empty384$strain_id[empty384$`384col` %in% c(1,24)] <- -2


transient_collection <- NULL
for (i in 1:3) {
  temp <- empty384
  temp$strain_id[temp$pos %in%
                       sample(temp$pos[is.na(temp$strain_id)],5)] <- 0
  transient_collection <- rbind(transient_collection, 
                                cbind(`384plate` = i, temp[,-1]))
  if (i < 3) {
    transient_collection$strain_id[transient_collection$`384plate` == i &
                                     is.na(transient_collection$strain_id)] <- 
      sample(orfs$strain_id[orfs$`384plate` %in% c(2,3,17,23:28) &
                              !is.na(orfs$category) &
                              orfs$category == 'Transient' &
                              orfs$strain_id %notin% transient_collection$strain_id], 303)
    
  } else {
    transient_collection$strain_id[transient_collection$`384plate` == i &
                                     is.na(transient_collection$strain_id)][1:57] <- 
      sample(orfs$strain_id[orfs$`384plate` %notin% c(2,3,17,23:28) &
                              !is.na(orfs$category) &
                              orfs$category == 'Transient' &
                              # orfs$strain_id %notin% transient_collection$strain_id &
                              !is.na(orfs$strain_id)])
    
    transient_collection$strain_id[transient_collection$`384plate` == i &
                                     is.na(transient_collection$strain_id)] <- 
      sample(orfs$strain_id[orfs$`384plate` %in% c(23:28) &
                              !is.na(orfs$category) &
                              orfs$category == 'Transient' &
                              # orfs$strain_id %notin% transient_collection$strain_id &
                              !is.na(orfs$strain_id)], 303 - 57)
  }
}
tecan_file <- merge(transient_collection %>% filter(strain_id %notin% c(0,-2)), orfs, by = 'strain_id', all.x = T,
                              suffixes = c('_target','_source'))
tecan_file <- tecan_file[order(tecan_file$`384plate_target`, tecan_file$`384col_target`,tecan_file$`384row_target`),]

transient_collection$new_strain_id <- NULL
transient_collection$new_strain_id[transient_collection$strain_id %notin% c(0,-2)] <- 
  (max(orfs$strain_id)+1):(max(orfs$strain_id)+
     length(transient_collection$strain_id[transient_collection$strain_id %notin% c(0,-2)]))

# write.csv(transient_collection, file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/transient_collection.csv')
# write.csv(tecan_file, file = '/home/sbp29/R/Projects/adaptivefitness/output/translatome/tecan_file.csv')

length(unique(transient_collection$strain_id))

transient_collection %>%
  ggplot(aes(x = `384col`, y = `384row`)) +
  geom_tile(aes(fill = new_strain_id), col = 'black') +
  facet_wrap(.~`384plate`)

head(transient_collection)
head(orfs)

db.file <- merge(transient_collection %>% filter(strain_id > 0), orfs[,c('strain_id','orf_name')], by = 'strain_id')
db.oldsid2newsid <- db.file[,c('strain_id','new_strain_id')]
db.file <- db.file[,c('orf_name','new_strain_id','384plate','384row','384col')]
colnames(db.file) <- c('orf_name','strain_id','384plate','384row','384col')

dbWriteTable(conn, 'TRANSIENT_COLLECTION', db.file, overwrite = T)
dbWriteTable(conn, 'TRANSIENT_COLLECTION_STRAIN_IDS', db.oldsid2newsid, overwrite = T)
dbWriteTable(conn, 'TRANSIENT_COLLECTION_SOURCE2TARGET', tecan_file[,c('strain_id','orf_name',
                                                                       '384plate_source','384row_source','384col_source',
                                                                       '384plate_target','384row_target','384col_target')],
             overwrite = T)
