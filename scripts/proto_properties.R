##### PROTO_PROPERTIES
##### proto-gene property evaluations
##### Author  : Saurin Parikh (dr.saurin.parikh@gmail.com)
##### Date    : 25/11/2018

##### INITIALIZE
#install.packages("RMariaDB")
#install.packages("ggplot2")
library(RMariaDB)
library(ggplot2)
source("R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")


##### DATA GRAB
query = dbSendQuery(conn, ("select distinct a.orf_name, b.aromaticity, b.aromaticity, b.pI, b.protogene
from PT2_pos2orf_name a, NT_PROPERTIES b
where a.orf_name = b.orf_name"))
nt_prop = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

query = dbSendQuery(conn, ("select b.*
from PT2R_PGAL_FS_6144_RES_eFDR a, NT_PROPERTIES b
where a.hours = 21 and a.effect_cs = 1
and a.orf_name = b.orf_name"))
ben_orfs = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

##### CLEAN UP

for (i in 1:length(nt_prop$protogene)){
  if (nt_prop$protogene[i] == 1) {
    nt_prop$proto[i] <- 'proto-gene'
  } else {
    nt_prop$proto[i] <- 'gene'
  }
}

for (i in 1:length(ben_orfs$protogene)){
  if (ben_orfs$protogene[i] == 1) {
    ben_orfs$proto[i] <- 'proto-gene'
  } else {
    ben_orfs$proto[i] <- 'gene'
  }
}

##### STATS

w_all <- wilcox.test(nt_prop$aromaticity[nt_prop$proto == 'proto-gene'],
                     nt_prop$aromaticity[nt_prop$proto == 'gene'])

w_ben <- wilcox.test(ben_orfs$aromaticity[ben_orfs$proto == 'proto-gene'],
                     ben_orfs$aromaticity[ben_orfs$proto == 'gene'])

##### BOX PLOTS

ggplot(nt_prop, aes(x = proto, y = aromaticity, fill=proto)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=4) +
  labs(title = sprintf("Overall \n p = %.5f", w_all$p.value), x = "ORF Type", y = "Aromaticity", fill = "Proto") +
  ylim(0,0.35) +
  theme_gray()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20),
        plot.title=element_text(size=20,hjust =.5),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.position="none")
ggsave("figs/overall_aromaticity_pt2.png",bg="transparent",height = 9, width = 12)

ggplot(ben_orfs, aes(x = proto, y = aromaticity, fill=proto)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=4) +
  labs(title = sprintf("pgalR \n p = %.5f", w_ben$p.value), x = "ORF Type", y = "Aromaticity", fill = "Proto") +
  ylim(0,.35) +
  theme_gray()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20),
        plot.title=element_text(size=20,hjust =.5),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.position="none")
ggsave("figs/pgalR_aromaticity.png",bg="transparent",height = 9, width = 12)

##### END
dbDisconnect(conn)


