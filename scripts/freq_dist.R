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

query = dbSendQuery(conn, ("select orf_name, fitness from YBR_RAF_6144_FITNESS
where hours = 30 and orf_name = 'BF_control' and fitness > 0"))
data.control = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

query = dbSendQuery(conn, ("select orf_name, fitness from YBR_RAF_6144_FITNESS
where hours = 30 and orf_name = 'YBR196C-A' and fitness > 0"))
data.orfs = dbFetch(query, n=-1)

if (dbHasCompleted(query)) {
  dbClearResult(query)
} else {
  print("Error Running Query")
}

#query = dbSendQuery(conn, ("select orf_name, fitness from PT_SA_CN_6144_FITNESS
#where hours = 12 and fitness between 0 and 0.90
#and orf_name in
#(select orf_name from PT_SA_CN_6144_RES_eFDR
#where hours = 12 and effect_cs = -1)
#order by fitness desc"))
#data.dels = dbFetch(query, n=-1)

#if (dbHasCompleted(query)) {
#  dbClearResult(query)
#} else {
#  print("Error Running Query")
#}

#data.dels$fitness[data.dels$fitness < 0]

##### PLOTS

#ggplot(data.orfs, aes(x=fitness, fill=orf_name)) +
#  geom_density(alpha = 0.3) +
#  geom_density(data = data.control, alpha = 0.3)

#ggplot() +
#  geom_density(data = data.dels, aes(x=fitness) ,fill="#D32F2F", color="#757575",alpha = 0.5) +
#  #scale_y_continuous(trans="log10") +
#  geom_density(data = data.orfs, aes(x=fitness), fill="#4CAF50", color="#757575", alpha = 0.5) +
#  geom_density(data = data.control, aes(x=fitness), fill="#303F9F", color="#757575", alpha = 0.5) +
#  scale_x_continuous(breaks = round(seq(0, 2, by = 0.1),1)) +
#  labs(title = "Probability Distribution", x = "Fitness", y = "Density", fill = "Strain Type", col = "Strain Type") +
#  #scale_fill_manual(labels = c("Deleterious", "Beneficial", "Reference"), values = c("#F44336", "#3F51B5", "#BDBDBD")) +
#  theme_gray()+
#  theme(axis.text=element_text(size=14),
#        axis.title=element_text(size=20),
#        plot.title=element_text(size=20,hjust =.5),
#        #panel.grid.major = element_blank(),
#        #panel.grid.minor = element_blank(),
#        panel.background = element_rect(fill = "transparent", colour = NA),
#        plot.background = element_rect(fill = "transparent", colour = NA))
#ggsave("figs/overall_distribution_PT_SA_noback.png",bg="transparent",height = 9, width = 12)

##### YBR DATA

ggplot() +
  geom_density(data = data.orfs, aes(x=fitness), fill="#4CAF50", color="#757575", alpha = 0.5) +
  geom_density(data = data.control, aes(x=fitness), fill="#303F9F", color="#757575", alpha = 0.5) +
  #scale_x_continuous(breaks = round(seq(0, 2, by = 0.1),1)) +
  labs(title = "Probability Distribution", x = "Fitness", y = "Density", fill = "Strain Type", col = "Strain Type") +
  #scale_fill_manual(labels = c("Deleterious", "Beneficial", "Reference"), values = c("#F44336", "#3F51B5", "#BDBDBD")) +
  xlim(.5, 1.5) +
  theme_gray()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20),
        plot.title=element_text(size=20,hjust =.5),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))
#ggsave("figs/overall_distribution_PT_SA_noback.png",bg="transparent",height = 9, width = 12)

##### END
dbDisconnect(conn)

