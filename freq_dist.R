##### FREQ_DIST
##### proto-gene frequency distribution plots
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 10/11/2018

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("ggplot2")
library(RMariaDB)
library(ggplot2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### FETCH DATA

query = dbSendQuery(conn, ("select orf_name, fitness from PT2_PGLU_FS_6144_FITNESS
where hours = 12 and orf_name = 'BF_control' and fitness > 0"))
data.control = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

query = dbSendQuery(conn, ("select orf_name, fitness from PT2_PGLU_FS_6144_FITNESS
where hours = 12 and fitness > 0
and orf_name in
(select orf_name from PT2_PGLU_FS_6144_RES_eFDR
where hours = 12 and effect_cs = 1)
order by fitness desc"))
data.orfs = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

query = dbSendQuery(conn, ("select orf_name, fitness from PT2_PGLU_FS_6144_FITNESS
where hours = 12 and fitness > 0
and orf_name in
(select orf_name from PT2_PGLU_FS_6144_RES_eFDR
where hours = 12 and effect_cs = -1)
order by fitness desc"))
data.dels = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

##### PLOTS

#ggplot(data.orfs, aes(x=fitness, fill=orf_name)) +
#  geom_density(alpha = 0.3) +
#  geom_density(data = data.control, alpha = 0.3)

ggplot(data.orfs, aes(x=fitness)) +
  geom_density(fill="#F44336", color="#757575",alpha = 0.3) +
  geom_density(data = data.control, fill="#3F51B5", color="#757575", alpha = 0.3) +
  geom_density(data = data.dels, fill="#BDBDBD", color="#757575", alpha = 0.3) +
  scale_fill_manual(name = "ORF Type", labels = c("beneficial","controls","deleterious")) +
  scale_x_continuous(breaks = round(seq(0, 2, by = 0.1),1)) +
  labs(title = "Fitness Distribution", x = "Fitness", y = "Density") +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=20),
        plot.title=element_text(size=20,hjust =.5))
  

##### END
dbDisconnect(conn)

