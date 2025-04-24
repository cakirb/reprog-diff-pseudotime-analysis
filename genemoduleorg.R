# gene module organization

# briggs-standard 
library(ggplot2)
library(dplyr)
agmat_brigst <- read.csv('./final_files/aggmat_final//briggs_std_aggmat.csv', row.names = 1)
genmod_brigst <- read.csv('./final_files/genemodules_final/briggs_std_geneModules.csv', row.names = 1)
  sort(apply(agmat_brigst, 1, min))
  sort(apply(agmat_brigst, 1, sd))

p1 <- pheatmap::pheatmap(agmat_brigst, clustering_method="ward.D2", cluster_cols = F)

plot(x$tree_row)
abline(h=1.45, col="red", lty=2, lwd=2)

groups <- as.factor(cutree(p1$tree_row, h=1.45))

agmat_brigst <- cbind(agmat_brigst, groups)
agmat_brigst <- tibble::rownames_to_column(agmat_brigst, 'module')            

agmat_brigst <- tidyr::pivot_longer(agmat_brigst, cols = starts_with('d'))
agmat_brigst$name <- factor(agmat_brigst$name, levels = c('d0', 'd5', 'd12'))

ggplot(agmat_brigst, aes(x = name, y = value, color = groups, group = module)) + geom_line() + 
#  gghighlight(X %in% paste('Module',brig_std_g[[i]])) + 
  geom_hline(aes(yintercept = 0), color = 'gray') + theme_minimal() + facet_wrap(~groups) + 
  xlab('Timepoint') + ylab('Aggregated Gene Counts')
ggsave('briggs_std_trends_7feb22.png')





#example go analysis (tran_a2s)

modules <- names(groups[groups == 1])
modules <- unlist(strsplit(modules, '[ ]'))
modules <- as.integer(modules[modules != 'Module'])
genes <- rownames(gen_mod[gen_mod$module %in% modules,])


library(gprofiler2)
x <- gprofiler2::gost(genes, "mmusculus", ordered_query = FALSE, 
                      custom_bg = rownames(cds), sources = "GO:BP") 
x <- x$result
addWorksheet(wb, paste0("cluster",i,"neg"))
writeData(wb, sheet = paste0("cluster",i,"neg"), x = x)
print(paste0("cluster",i," completed!!"))



# coloring module groups
cds <- readRDS("./final_files/traj_rds_final/tran_a2s_noirr_monocle3v1_5nov21.rds")
gene_module_df <- readr::read_csv("./final_files/genemodules_final/tran_a2s_geneModules_8feb22.csv")
agg_mat <- read.csv("./final_files/aggmat_final/tran_a2s_aggmat_8feb22.csv", row.names = 1)
hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
groups <- as.factor(cutree(hm$tree_row, k=4))
groups <- as.data.frame(groups)
groups$module <- 1:nrow(groups)

for (j in 1:nrow(groups)){
  gene_module_df$module[gene_module_df$module == groups[j,2]] <- groups[j,1]
}


rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(1,2,3,4)),
           label_cell_groups=FALSE, label_branch_points = F,
           label_roots = F, label_leaves = F)
ggsave('tran_a2s_modulesPlotted_8feb22.png',
       height = 19.05, width = 16.52, units = 'cm')



# extended go analysis function (tran_a2s)
library(gprofiler2)
library(openxlsx)


wb <- createWorkbook()

for (k in unique(gene_module_df$module)){
  genes <- gene_module_df[gene_module_df$module == k,]$id
  x <- gprofiler2::gost(genes, "mmusculus", ordered_query = FALSE, 
                        custom_bg = rownames(cds), sources = "GO:BP") 
  #x <- x$result
  p <-gostplot(x, capped = F, interactive = F)
  publish_gostplot(p,  highlight_terms = x$result[1:10,])
  
  addWorksheet(wb, paste0("cluster",i,"neg"))
  writeData(wb, sheet = paste0("cluster",i,"neg"), x = x)
  print(paste0("cluster",i," completed!!"))
}

library(clusterProfiler)
library(org.Mm.eg.db)
for (i in unique(gene_module_df$module)){
  print(paste(i,'started!'))
  genes <- gene_module_df[gene_module_df$module == i,]$id
  ego <- enrichGO(gene          = genes,
                  universe      = rownames(cds),
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "SYMBOL", 
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  if (dim(summary(ego))[1] == 0){next}
  assign(paste0('x',i),barplot(ego, showCategory=10))
  print(paste(i,'done!'))
  #ggsave(paste0('tran_a2s_dotplot.png'), width = 10)
  }

for (i in unique(gene_module_df$module)){
  print(paste(i,'started!'))
  genes <- gene_module_df[gene_module_df$module == i,]$id
  ego <- enrichGO(gene          = genes,
                  universe      = rownames(cds),
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "SYMBOL", 
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  if (dim(summary(ego))[1] == 0){next}
  assign(paste0('x',i),dotplot(ego, showCategory=10))
  print(paste(i,'done!'))
  #ggsave(paste0('tran_a2s_dotplot.png'), width = 10)
}

for (i in list.files(path="./final_files/traj_rds_final")){
  print(paste(i,'started!'))
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  agg_mat <- read.csv(paste0('./final_files/aggmat_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_aggmat_8feb22.csv'), row.names = 1)
  gene_module <- read.csv(paste0('./final_files/genemodules_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_geneModules_8feb22.csv'), row.names = 1)
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module$module[gene_module$module == groups[j,2]] <- groups[j,1]
  }
  if (startsWith(i,'briggs_direct')) { 
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="d5"] <- "d4"
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd4', 'd11'))
  } else if (startsWith(i,'briggs_std')) {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd5', 'd12'))
  } else if (startsWith(i,'tran_a2s')) {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd2', 'd4', 'd6', 'mESC'))
  } else if (startsWith(i,'tran_fbs')){
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC'))
  } else if (startsWith(i,'treutlein')){
    cds@colData$Timepoint <- factor(cds@colData$time_point, levels = c('0', '2', '5', '20', '22'))
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="0"] <- "d0"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="2"] <- "d2"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="5"] <- "d5"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="20"] <- "d20"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="22"] <- "d22"
  }
  rowData(cds)$gene_name <- rownames(cds)
  rowData(cds)$gene_short_name <- rowData(cds)$gene_name
  gene_module <- tibble::rownames_to_column(gene_module, 'genes') 
  plot_cells(cds, color_cells_by = 'Timepoint', label_cell_groups = FALSE)
  ggsave(paste0('./final_files/plots/',i,'_umap_13feb22.png'), 
         width = 17.52, height = 18.92, units = 'cm', scale = 0.8)
  plot_cells(cds,
             genes=gene_module %>% filter(module %in% c(1,2,3,4)) %>% select(1,2),
             label_cell_groups=FALSE, label_branch_points = F,
             label_roots = F, label_leaves = F)
  ggsave(paste0('./final_files/plots/',i,'_umap_modulegroups_13feb22.png'), 
         width = 17.52, height = 18.92, units = 'cm', scale = 0.8)
  print(paste(i,'done!'))
}
if (startsWith(i,'briggs_direct')){colnames(agg_mat) <- c('d0', 'd4', 'd11')}


library(patchwork)
library(ggplot2)
library(stringr)
library(org.Mm.eg.db)

# GO enrichment
for (i in list.files(path="./final_files/traj_rds_final")){
  print(paste(i,'started!'))
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  agg_mat <- read.csv(paste0('./final_files/aggmat_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_aggmat_8feb22.csv'), row.names = 1)
  gene_module <- read.csv(paste0('./final_files/genemodules_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_geneModules_8feb22.csv'), row.names = 1)
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module$module[gene_module$module == groups[j,2]] <- groups[j,1]
  }
  rm(agg_mat)
  rm(hm)
  rm(groups)
  for (k in unique(gene_module$module)){
    print(paste('Module Group',k,'started!'))
    genes <- rownames(gene_module[gene_module$module == k,])
    ego <- enrichGO(gene          = genes,
                    #universe      = rownames(cds),
                    OrgDb         = org.Mm.eg.db,
                    keyType       = "SYMBOL", 
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
    if (dim(summary(ego))[1] == 0){assign(paste0('p',k),grid::textGrob(paste('No significant\nGO term found for\nmodule group',k),))}
    else{assign(paste0('p',k),dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                  scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
  }
  pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
  ggsave(paste0('./final_files/GOplots/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_dotplotss.png'), scale = 1.5)
  print(paste(i,'finished!'))
}

#KEGG enrichment
library(AnnotationDbi)
library(org.Mm.eg.db)
for (i in list.files(path="./final_files/traj_rds_final")){
  print(paste(i,'started!'))
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  agg_mat <- read.csv(paste0('./final_files/aggmat_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_aggmat_8feb22.csv'), row.names = 1)
  gene_module <- read.csv(paste0('./final_files/genemodules_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_geneModules_8feb22.csv'), row.names = 1)
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module$module[gene_module$module == groups[j,2]] <- groups[j,1]
  }
  rm(agg_mat)
  rm(hm)
  rm(groups)
  univ <- mapIds(org.Mm.eg.db, keys = rownames(cds), column = "ENTREZID", keytype = "SYMBOL")
  for (k in unique(gene_module$module)){
    print(paste('Module Group',k,'started!'))
    genes <- rownames(gene_module[gene_module$module == k,])
    genes <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
    ego <- enrichKEGG(gene          = genes,
                      universe      = univ,
                      organism      = "mmu",
                      keyType       = "kegg",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
    if (dim(summary(ego))[1] == 0){assign(paste0('p',k),grid::textGrob(paste('No significant\nKEGG term found for\nmodule group',k),))}
    else{assign(paste0('p',k),dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                  scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
  }
  pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
  ggsave(paste0('./final_files/GOplots/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_dotplotss_kegg.png'), scale = 1.5)
  print(paste(i,'finished!'))
}

#WikiPathways enrichment
library(AnnotationDbi)
library(org.Mm.eg.db)
for (i in list.files(path="./final_files/traj_rds_final")){
  print(paste(i,'started!'))
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  agg_mat <- read.csv(paste0('./final_files/aggmat_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_aggmat_8feb22.csv'), row.names = 1)
  gene_module <- read.csv(paste0('./final_files/genemodules_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_geneModules_8feb22.csv'), row.names = 1)
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module$module[gene_module$module == groups[j,2]] <- groups[j,1]
  }
  rm(agg_mat)
  rm(hm)
  rm(groups)
  univ <- mapIds(org.Mm.eg.db, keys = rownames(cds), column = "ENTREZID", keytype = "SYMBOL")
  for (k in unique(gene_module$module)){
    print(paste('Module Group',k,'started!'))
    genes <- rownames(gene_module[gene_module$module == k,])
    genes <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
    ego <- enrichWP(genes, organism = "Mus musculus")
    if (dim(summary(ego))[1] == 0){assign(paste0('p',k),grid::textGrob(paste('No significant\nWikiPathways term found for\nmodule group',k),))}
    else{assign(paste0('p',k),dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                  scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
  }
  pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
  ggsave(paste0('./final_files/GOplots/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_dotplotss_wp.png'), scale = 1.5)
  print(paste(i,'finished!'))
}


#gprofiler2 !! (august2022) -----------------
library(gprofiler2)
for (i in list.files(path="./final_files/traj_rds_final")){
  print(paste(i,'started!'))
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  agg_mat <- read.csv(paste0('./final_files/aggmat_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_aggmat_8feb22.csv'), row.names = 1)
  gene_module <- read.csv(paste0('./final_files/genemodules_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_geneModules_8feb22.csv'), row.names = 1)
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module$module[gene_module$module == groups[j,2]] <- groups[j,1]
  }
  rm(agg_mat)
  rm(hm)
  rm(groups)
  for (k in unique(gene_module$module)){
    print(paste('Module Group',k,'started!'))
    genes <- rownames(gene_module[gene_module$module == k,])
    ego <- gprofiler2::gost(genes, "mmusculus", ordered_query = FALSE, 
                            custom_bg = rownames(cds), sources = "GO:BP")
    if (is.null(ego)){assign(paste0('p',k),grid::textGrob(paste('No significant\nGO term found for\nmodule group',k),))}else{
      assign(paste0('p',k),publish_gosttable(ego, show_columns = c("term_name"), highlight_terms = ego$result$term_id[1:10]))}
  }
  pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
  ggsave(paste0('./final_files/GOplots/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_dotplotss.png'), scale = 1.5)
  print(paste(i,'finished!'))
}


#KEGG enrichment
library(AnnotationDbi)
library(org.Mm.eg.db)
for (i in list.files(path="./final_files/traj_rds_final")){
  print(paste(i,'started!'))
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  agg_mat <- read.csv(paste0('./final_files/aggmat_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_aggmat_8feb22.csv'), row.names = 1)
  gene_module <- read.csv(paste0('./final_files/genemodules_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_geneModules_8feb22.csv'), row.names = 1)
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module$module[gene_module$module == groups[j,2]] <- groups[j,1]
  }
  rm(agg_mat)
  rm(hm)
  rm(groups)
  univ <- mapIds(org.Mm.eg.db, keys = rownames(cds), column = "ENTREZID", keytype = "SYMBOL")
  for (k in unique(gene_module$module)){
    print(paste('Module Group',k,'started!'))
    genes <- rownames(gene_module[gene_module$module == k,])
    ego <- gprofiler2::gost(genes, "mmusculus", ordered_query = FALSE, 
                            custom_bg = rownames(cds), sources = "KEGG")
    if (is.null(ego)){assign(paste0('p',k),grid::textGrob(paste('No significant\nKEGG term found for\nmodule group',k),))}else{
      assign(paste0('p',k),publish_gosttable(ego, show_columns = c("term_name"), highlight_terms = ego$result$term_id[1:10]))}
  }
  pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
  ggsave(paste0('./final_files/GOplots/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_dotplotss_kegg.png'), scale = 1.5)
  print(paste(i,'finished!'))
}

#WikiPathways enrichment
library(AnnotationDbi)
library(org.Mm.eg.db)
for (i in list.files(path="./final_files/traj_rds_final")){
  print(paste(i,'started!'))
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  agg_mat <- read.csv(paste0('./final_files/aggmat_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_aggmat_8feb22.csv'), row.names = 1)
  gene_module <- read.csv(paste0('./final_files/genemodules_final/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_geneModules_8feb22.csv'), row.names = 1)
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module$module[gene_module$module == groups[j,2]] <- groups[j,1]
  }
  rm(agg_mat)
  rm(hm)
  rm(groups)
  univ <- mapIds(org.Mm.eg.db, keys = rownames(cds), column = "ENTREZID", keytype = "SYMBOL")
  for (k in unique(gene_module$module)){
    print(paste('Module Group',k,'started!'))
    genes <- rownames(gene_module[gene_module$module == k,])
    ego <- gprofiler2::gost(genes, "mmusculus", ordered_query = FALSE, 
                            custom_bg = rownames(cds), sources = "WP")
    if (is.null(ego)){assign(paste0('p',k),grid::textGrob(paste('No significant\nWikiPathways term found for\nmodule group',k),))}else{
      assign(paste0('p',k),publish_gosttable(ego, show_columns = c("term_name"), highlight_terms = ego$result$term_id[1:10]))}
  }
  pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
  ggsave(paste0('./final_files/GOplots/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_dotplotss_wp.png'), scale = 1.5)
  print(paste(i,'finished!'))
}


# GO enrichment (new tran-a2s only)
  i <- "tran_a2s_noirr_monocle3v1_5nov21.rds"     
  print(paste(i,'started!'))
  agg_mat <- read.csv('tran_a2s_aggmat_8feb22.csv', row.names = 1)
  gene_module <- read.csv('tran_a2s_geneModules_8feb22.csv', row.names = 1)
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module$module[gene_module$module == groups[j,2]] <- groups[j,1]
  }
  rm(agg_mat)
  rm(hm)
  rm(groups)
  for (k in unique(gene_module$module)){
    print(paste('Module Group',k,'started!'))
    genes <- rownames(gene_module[gene_module$module == k,])
    ego <- enrichGO(gene          = genes,
                    universe      = rownames(cds),
                    OrgDb         = org.Mm.eg.db,
                    keyType       = "SYMBOL", 
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
    if (dim(summary(ego))[1] == 0){assign(paste0('p',k),grid::textGrob(paste('No significant\nGO term found for\nmodule group',k),))}
    else{assign(paste0('p',k),dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                  scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
  }
  pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
  ggsave(paste0('./final_files/GOplots/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_dotplotss.png'), scale = 1.5)
  print(paste(i,'finished!'))
  
  
# to produce gene profile plots
    
for (i in list.files(path=".", pattern = "graphtest_2aug22.csv")){
  print(paste(i,'started!'))
  
  gt_result <- read.csv(i, row.names = 1)
  setname <- strsplit(i, "_")[[1]][1]
  if (startsWith(i, "tran")){setname <- paste0(setname,"_noirr")}
  cds <- readRDS(paste0("./final_files/traj_rds_final/",setname,"_monocle3v1_5nov21.rds"))
  cds <- preprocess_cds(cds, num_dim = 50)
  
  if (startsWith(i,'briggs_direct')) { 
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="d5"] <- "d4"
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd4', 'd11'))
  } else if (startsWith(i,'briggs_std')) {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd5', 'd12'))
  } else if (startsWith(i,'tran_a2s')) {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd2', 'd4', 'd6', 'mESC'))
  } else if (startsWith(i,'tran_fbs')){
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC'))
  } else if (startsWith(i,'treutlein')){
    cds@colData$Timepoint <- factor(cds@colData$time_point, levels = c('0', '2', '5', '20', '22'))
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="0"] <- "d0"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="2"] <- "d2"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="5"] <- "d5"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="20"] <- "d20"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="22"] <- "d22"
  }
  
  traj_genes <- row.names(subset(gt_result, q_value < 0.05))
  
  gene_module_df <- find_gene_modules(cds[traj_genes,])#, resolution=c(10^seq(-6,-1)))
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=colData(cds)$Timepoint)
  agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F,
                           file=paste0(setname,"_genemodules_2aug22.png"))
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module_df$module[gene_module_df$module == groups[j,2]] <- groups[j,1]
  }
  
  agmat2 <- cbind(agg_mat, groups)
  agmat2 <- tibble::rownames_to_column(agmat2, 'module')
  agmat2 <- tidyr::pivot_longer(agg_mat, cols = starts_with('d'))
  agmat2$name <- factor(agmat2$name, levels = c('d0', 'd5', 'd12'))
  
  ggplot(agmat2, aes(x = name, y = value, color = groups, group = module)) + geom_line() + 
    #  gghighlight(X %in% paste('Module',brig_std_g[[i]])) + 
    geom_hline(aes(yintercept = 0), color = 'gray') + theme_minimal() + facet_wrap(~groups) + 
    xlab('Timepoint') + ylab('Aggregated Gene Counts')
  ggsave(paste0(setname, "_expProfile_4aug22.png"))
  
  print(paste(i,'done!'))
  
}
  
for (i in list.files(path=".", pattern = "graphtest_2aug22.csv")){
  print(paste(i,'started!'))
  
  gt_result <- read.csv(i, row.names = 1)
  setname <- strsplit(i, "_")[[1]][1]
  if (startsWith(i, "tran")){setname <- paste0(setname,"_noirr")}
  cds <- readRDS(paste0("./final_files/traj_rds_final/",setname,"_monocle3v1_5nov21.rds"))
  cds <- preprocess_cds(cds, num_dim = 50)
  
  if (startsWith(i,'briggs_direct')) { 
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="d5"] <- "d4"
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd4', 'd11'))
  } else if (startsWith(i,'briggs_std')) {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd5', 'd12'))
  } else if (startsWith(i,'tran_a2s')) {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd2', 'd4', 'd6', 'mESC'))
  } else if (startsWith(i,'tran_fbs')){
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC'))
  } else if (startsWith(i,'treutlein')){
    cds@colData$Timepoint <- factor(cds@colData$time_point, levels = c('0', '2', '5', '20', '22'))
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="0"] <- "d0"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="2"] <- "d2"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="5"] <- "d5"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="20"] <- "d20"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="22"] <- "d22"
  }
  
  traj_genes <- row.names(subset(gt_result, q_value < 0.05))
  
  gene_module_df <- find_gene_modules(cds[traj_genes,])#, resolution=c(10^seq(-6,-1)))
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=colData(cds)$Timepoint)
  agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F,
                           file=paste0(setname,"_genemodules_2aug22.png"))
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module_df$module[gene_module_df$module == groups[j,2]] <- groups[j,1]
  }
  ggplot(agg_mat, aes(x = name, y = value, color = groups, group = module)) + geom_line() + 
    #  gghighlight(X %in% paste('Module',brig_std_g[[i]])) + 
    geom_hline(aes(yintercept = 0), color = 'gray') + theme_minimal() + facet_wrap(~groups) + 
    xlab('Timepoint') + ylab('Aggregated Gene Counts')
  ggsave(paste0(setname, "_expProfile_4aug22.png"))
  for (k in unique(gene_module_df$module)){
    print(paste('Module Group',k,'started!'))
    genes <- gene_module_df[gene_module_df$module == k,]$id
    ego <- enrichGO(gene          = genes,
                    universe      = rownames(cds),
                    OrgDb         = org.Mm.eg.db,
                    keyType       = "SYMBOL", 
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
    if (dim(summary(ego))[1] == 0){assign(paste0('p',k),grid::textGrob(paste('No significant\nGO term found for\nmodule group',k),))}
    else{assign(paste0('p',k),dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                  scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
  }
  pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
  ggsave(paste0('./final_files/GOplots/',paste(strsplit(i, '_')[[1]][1:2], collapse = '_'),'_dotplotss_2aug22.png'), scale = 1.5)
  print(paste(i,'finished!'))
  print(paste(i,'done!'))
  
  # gene_module <- tibble::rownames_to_column(gene_module, 'genes') 
  # plot_cells(cds, color_cells_by = 'Timepoint', label_cell_groups = FALSE)
  # ggsave(paste0('./final_files/plots/',i,'_umap_13feb22.png'), 
  #        width = 17.52, height = 18.92, units = 'cm', scale = 0.8)
  # plot_cells(cds,
  #            genes=gene_module %>% filter(module %in% c(1,2,3,4)) %>% select(1,2),
  #            label_cell_groups=FALSE, label_branch_points = F,
  #            label_roots = F, label_leaves = F)
  # ggsave(paste0('./final_files/plots/',i,'_umap_modulegroups_13feb22.png'), 
  #        width = 17.52, height = 18.92, units = 'cm', scale = 0.8)
  
}
  