

goe.es <- data.frame()
kegg.es <- data.frame()
for (id_ref in c('ONE_GA','TWO_GA','TWO_DM')) {
  allgenes <- unique(data.diff2$orf_name[data.diff2$id_ref == id_ref])
  allgenes <- bitr(allgenes, fromType = "ORF",
                   toType = c("ENTREZID","GENENAME","ENSEMBL"),
                   OrgDb = org.Sc.sgd.db)
  allgenes <- allgenes[!is.na(allgenes$ENSEMBL),]
  for (id_cond in unique(data.diff2$id_cond[data.diff2$id_ref == id_ref])) {
    for (es in seq(0,2,0.2)) {
      for (d in c(-1,1)) {
        temp.deg <- data.diff2$orf_name[data.diff2$id_ref == id_ref & data.diff2$id_cond == id_cond &
                                          data.diff2$fitness_diff*d >= es & !is.na(data.diff2$fitness_diff) &
                                          data.diff2$orf_name != 'BF_control' & !is.na(data.diff2$orf_name)]
        if (length(temp.deg) != 0) {
          temp.deg <- bitr(temp.deg, fromType = "ORF",
                           toType = c("ENTREZID","GENENAME","ENSEMBL"),
                           OrgDb = org.Sc.sgd.db)
          temp.deg <- temp.deg[!is.na(temp.deg$ENSEMBL),]
          
          temp.goe <- enrichGO(gene          = temp.deg$ENSEMBL,
                               universe      = allgenes$ENSEMBL,
                               OrgDb         = org.Sc.sgd.db,
                               keyType       = "ENSEMBL",
                               ont           = "ALL",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05)
          
          if (length(temp.goe) == 1) {
            if (dim(temp.goe)[1] == 0) {
            } else{
              goe.es <- rbind(goe.es, data.frame(temp.goe, id_ref = id_ref, id_cond = id_cond, effect_size = es, direction = d))
            }
          }
          
          temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                                  universe     = allgenes$ENSEMBL,
                                  organism     = 'sce',
                                  pvalueCutoff = 0.05)
          if (length(temp.kegg) == 1) {
            if (dim(temp.kegg)[1] == 0) {
            } else{
              kegg.es <- rbind(kegg.es, data.frame(temp.kegg, id_ref = id_ref, id_cond = id_cond, effect_size = es, direction = d))
            }
          }
        }
      }
    }
  }
}
goe.es$GeneRatio <- as.numeric(str_split(goe.es$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe.es$GeneRatio,'/',simplify = T)[,2])
goe.es$BgRatio <- as.numeric(str_split(goe.es$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe.es$BgRatio,'/',simplify = T)[,2])
goe.es$GO <- paste0(goe.es$ONTOLOGY, '_', goe.es$Description)

kegg.es$GeneRatio <- as.numeric(str_split(kegg.es$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg.es$GeneRatio,'/',simplify = T)[,2])
kegg.es$BgRatio <- as.numeric(str_split(kegg.es$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg.es$BgRatio,'/',simplify = T)[,2])

goe.es <- goe.es[order(goe$id_ref,goe$id_cond,goe.es$direction,goe.es$effect_size,-goe.es$GeneRatio,-goe.es$Count,goe.es$qvalue),]
kegg.es <- kegg.es[order(goe$id_ref,goe$id_cond,kegg.es$direction,kegg.es$effect_size,-kegg.es$GeneRatio,-kegg.es$Count,kegg.es$qvalue),]

write.csv(goe.es, file = 'output/translatome/go_es_enrichments.csv')
write.csv(kegg.es, file = 'output/translatome/kegg_es_enrichments.csv')
