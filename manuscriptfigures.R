## manuscript figures

# GO enrichment
i <- "tran_a2s_noirr_monocle3v1_5nov21.rds" 

library(org.Mm.eg.db)
library(clusterProfiler)
library(patchwork)
library(monocle3)
library(tidyverse)

for (set in list.files('./graph_test_outputs/')){
  setname <- strsplit(set, "_")[[1]][1]
  print(paste0('Reading graph test file for [[',setname,']]...'))
  gt_result <- read.csv(paste0('./graph_test_outputs/',set), row.names = 1)
  
  if (startsWith(set, "tran")){setname <- paste0(setname,"_noirr")}
  
  print('Reading Monocle3 rds file...')
  cds <- readRDS(paste0("./final_files/traj_rds_final/",setname,"_monocle3v1_5nov21.rds"))
  cds <- preprocess_cds(cds, num_dim = 50)
  
  traj_genes <- row.names(subset(gt_result, q_value < 0.05))
  
  print('Calculating gene modules...')
  gene_module_df <- find_gene_modules(cds[traj_genes,], resolution=c(10^seq(-6,-1)))
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=colData(cds)$time_point)
  agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  write.csv(agg_mat,paste0("./final_files_nov22/aggmat/",setname,"_aggmat_7nov22.csv"))

  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module_df$module[gene_module_df$module == groups[j,2]] <- groups[j,1]
  }
  rm(agg_mat)
  rm(hm)
  rm(groups)
  write.csv(gene_module_df,paste0("./final_files_nov22/genemodule/",setname,"_genemodule_7nov22.csv"))
  for (a in c('GO','KEGG','WP')){
    print(paste0('...',a, '... '))
    for (k in unique(gene_module_df$module)){
      print(paste('Module Group',k,'started!'))
      genes <- gene_module_df[gene_module_df$module == k,]$id
      if (a == 'GO'){
        ego <- enrichGO(gene        = genes,
                      universe      = rownames(cds),
                      OrgDb         = org.Mm.eg.db,
                      keyType       = "SYMBOL", 
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
        if (dim(summary(ego))[1] == 0){
          assign(paste0('p',k),ggplot() + annotate("text", x = 1, y = 1, size = 5, label = paste('No significant\nGO term found for\nmodule group',k)) + theme_void())
          } else { assign(paste0('p',k),
                         dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                    scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
        }
      else if (a == 'KEGG'){
        genes <- gene_module_df[gene_module_df$module == k,]$id
        genes <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
        univ <- mapIds(org.Mm.eg.db, keys = rownames(cds), column = "ENTREZID", keytype = "SYMBOL")
        ego <- enrichKEGG(gene          = genes,
                          universe      = univ,
                          organism      = "mmu",
                          keyType       = "kegg",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05)
        if (dim(summary(ego))[1] == 0){
          assign(paste0('p',k),ggplot() + annotate("text", x = 1, y = 1, size = 5, label = paste('No significant\nKEGG term found for\nmodule group',k)) + theme_void())
          } else { assign(paste0('p',k),
                          dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                      scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
        }
      else if (a == 'WP'){
        genes <- gene_module_df[gene_module_df$module == k,]$id
        genes <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
        ego <- enrichWP(genes, organism = "Mus musculus")
        if (dim(summary(ego))[1] == 0){
          assign(paste0('p',k),ggplot() + annotate("text", x = 1, y = 1, size = 5, label = paste('No significant\nWikiPathways term found for\nmodule group',k)) + theme_void())
          } else { assign(paste0('p',k),
                          dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                      scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
        }
    }
    pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
    ggsave(paste0('./final_files_nov22/clusterProfiler_plots/',paste(setname, a, '22nov22', sep = '_'),'_dotplots.png'), scale = 1.5)
    print(paste(setname,'finished!'))
  }
}
  


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







## thesis figures (last updated: 7 nov)

# running everything once again to confirm

# GO enrichment (new tran-a2s only)

setwd("~/Documents/thesispaper")

library(org.Mm.eg.db)
library(clusterProfiler)
library(patchwork)
library(monocle3)
library(tidyverse)

for (set in list.files('./graph_test_outputs/')){
  setname <- strsplit(set, "_")[[1]][1]
  
  if (startsWith(set, "tran")){setname <- paste0(setname,"_noirr")}
  
  # load gene modules
  gene_module_df <- read.csv(paste0('./genemodule/',setname,'_genemodule_7nov22.csv'))
  # prepare cell metadata (timepoints) 
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=colData(cds)$time_point)

  agg_mat <- read.csv(paste0("./final_files_nov22/aggmat/",setname,"_aggmat_7nov22.csv"))
  
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  groups$module <- 1:nrow(groups)
  for (j in 1:nrow(groups)){
    gene_module_df$module[gene_module_df$module == groups[j,2]] <- groups[j,1]
  }
  rm(agg_mat)
  rm(hm)
  rm(groups)
  
  
  write.csv(gene_module_df,paste0("./final_files_nov22/genemodule/",setname,"_genemodule_7nov22.csv"))
  for (a in c('GO','KEGG','WP')){
    print(paste0('...',a, '... '))
    for (k in unique(gene_module_df$module)){
      print(paste('Module Group',k,'started!'))
      genes <- gene_module_df[gene_module_df$module == k,]$id
      if (a == 'GO'){
        ego <- enrichGO(gene        = genes,
                        universe      = rownames(cds),
                        OrgDb         = org.Mm.eg.db,
                        keyType       = "SYMBOL", 
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
        if (dim(summary(ego))[1] == 0){
          assign(paste0('p',k),ggplot() + annotate("text", x = 1, y = 1, size = 5, label = paste('No significant\nGO term found for\nmodule group',k)) + theme_void())
        } else { assign(paste0('p',k),
                        dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                          scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
      }
      else if (a == 'KEGG'){
        genes <- gene_module_df[gene_module_df$module == k,]$id
        genes <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
        univ <- mapIds(org.Mm.eg.db, keys = rownames(cds), column = "ENTREZID", keytype = "SYMBOL")
        ego <- enrichKEGG(gene          = genes,
                          universe      = univ,
                          organism      = "mmu",
                          keyType       = "kegg",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05)
        if (dim(summary(ego))[1] == 0){
          assign(paste0('p',k),ggplot() + annotate("text", x = 1, y = 1, size = 5, label = paste('No significant\nKEGG term found for\nmodule group',k)) + theme_void())
        } else { assign(paste0('p',k),
                        dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                          scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
      }
      else if (a == 'WP'){
        genes <- gene_module_df[gene_module_df$module == k,]$id
        genes <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
        ego <- enrichWP(genes, organism = "Mus musculus")
        if (dim(summary(ego))[1] == 0){
          assign(paste0('p',k),ggplot() + annotate("text", x = 1, y = 1, size = 5, label = paste('No significant\nWikiPathways term found for\nmodule group',k)) + theme_void())
        } else { assign(paste0('p',k),
                        dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                          scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
      }
    }
    pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
    ggsave(paste0('./final_files_nov22/clusterProfiler_plots/',paste(setname, a, '22nov22', sep = '_'),'_dotplots.png'), scale = 1.5)
    print(paste(setname,'finished!'))
  }
}



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

for (i in list.files(path=".", pattern = "graphtest_2aug22.csv")){
  print(paste(i,'started!'))
  
  gt_result <- read.csv(i, row.names = 1)
  setname <- strsplit(i, "_")[[1]][1]
  if (startsWith(i, "tran")){setname <- paste0(setname,"_noirr")}
  cds <- readRDS(paste0("./final_files/traj_rds_final/",setname,"_monocle3v1_5nov21.rds"))
  cds <- preprocess_cds(cds, num_dim = 50)
  
  if (startsWith(i,'briggs_direct')) { 
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="d5"] <- "d4"
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd5', 'd11'))
  } else if (startsWith(i,'briggs_std')) {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd6 ', 'd12'))
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


library(monocle3)

# for (file in list.files('./final_files_nov22/aggmat/')){
#   agg_mat <- read.csv(paste0('./final_files_nov22/aggmat/',file), row.names = 1)
#   if (file == 'treutlein_aggmat_7nov22.csv'){
#     for (c("d0", "d2", "d5", "d20", "d22")
#   }
#   hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
# }

for (file in list.files('./final_files_nov22/aggmat/')){
  agg_mat <- read.csv(paste0('./final_files_nov22/aggmat/',file), row.names = 1)
  groups <- as.factor(cutree(hm$tree_row, k=4))
  groups <- as.data.frame(groups)
  pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F, annotation_row = groups)
}

gene_module_df <- find_gene_modules(cds[traj_genes,])#, resolution=c(10^seq(-6,-1)))


#---


setwd("~/Documents/thesispaper")

library(ggplot2)
library(monocle3)

for (i in list.files('./final_files/traj_rds_final/')){
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  if (setname %in% c('treutlein', 'brigdir', 'brigstd')){
    next
  }
  
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  
  gt_result <- read.csv(paste0('./graph_test_outputs/', setname, '_graphtest_4aug22.csv'), row.names = 1)
  cds <- preprocess_cds(cds, num_dim = 50)
  
  if (setname == 'brigdir') { 
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd5', 'd11'))
    levelss <- c('d0', 'd5', 'd11')
  } else if (setname == 'brigstd') {
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="d5"] <- "d6"
    levelss <- c('d0', 'd6', 'd12')
  } else if (setname == 'trana2s') {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint)
    levelss <- c('d0', 'd2', 'd4', 'd6', 'mESC')
  } else if (setname == 'tranfbs') {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC'))
    levelss <- c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC')
  } else if (setname == 'treutlein') {
    cds@colData$Timepoint <- factor(cds@colData$time_point, levels = c('0', '2', '5', '20', '22'))
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="0"] <- "d0"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="2"] <- "d2"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="5"] <- "d5"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="20"] <- "d20"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="22"] <- "d22"
    levelss <- c('d0', 'd2', 'd5', 'd20', 'd22')
  }
  
  traj_genes <- row.names(subset(gt_result, q_value < 0.05))
  
  gene_module_df <- find_gene_modules(cds[traj_genes,], resolution=c(10^seq(-6,-1)))
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=colData(cds)$Timepoint)
  
  agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  
  write.csv(agg_mat,paste0("./final_files_nov22/aggmat_jan23/",setname,"_aggmat_16jan23.csv"))
  
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  moduleGroup <- as.factor(cutree(hm$tree_row, k=4))
  gene_module_df$module4 <- moduleGroup[as.character(gene_module_df$module)]
  moduleGroup <- as.data.frame(moduleGroup)
  
  pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F, annotation_row = moduleGroup, fontsize = 7,
                     file = paste0('./final_files_nov22/heatmaps_final_jan23/',setname,'_heatmap.png'))
  
  write.csv(gene_module_df,paste0("./final_files_nov22/genemodule_jan23/",setname,"_genemodule_16jan23.csv"))
  
  agmat2 <- cbind(agg_mat, moduleGroup)
  agmat2 <- tibble::rownames_to_column(agmat2, 'module')
  agmat2 <- tidyr::pivot_longer(agmat2, cols = levelss)
  agmat2$name <- factor(agmat2$name, levels = levelss)
  
  ggplot(agmat2, aes(x = name, y = value, color = moduleGroup, group = module)) + geom_line() + 
    #  gghighlight(X %in% paste('Module',brig_std_g[[i]])) + 
    geom_hline(aes(yintercept = 0), color = 'gray') + theme_minimal() + facet_wrap(~moduleGroup) + 
    xlab('Timepoint') + ylab('Aggregated Gene Counts')
  ggsave(paste0('./final_files_nov22/genemodule_profile_jan23/', setname, '_profile_16jan23.png'))
}
  

setwd("~/Documents/thesispaper")

library(ggplot2)
for (i in list.files('./final_files_nov22/aggmat_jan23/')){
  agg_mat <- read.csv(paste0('./final_files_nov22/aggmat_jan23/',i), row.names = 1)
  setname <- strsplit(i, "_")[[1]][1]
  
  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  moduleGroup <- as.factor(cutree(hm$tree_row, k=4))
  moduleGroup <- as.data.frame(moduleGroup)
  
  agmat2 <- cbind(agg_mat, moduleGroup)
  agmat2 <- tibble::rownames_to_column(agmat2, 'module')
  agmat2 <- tidyr::pivot_longer(agmat2, cols = starts_with('d'))
  
  agmat2$name <- factor(agmat2$name, levels = c('d0', 'd5', 'd12'))
  
  
  if (startsWith(i,'brigdir')) { 
    levelss <- c('d0', 'd5', 'd11')
  } else if (startsWith(i,'brigstd')) {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd6 ', 'd12'))
    levelss <- c('d0', 'd6 ', 'd12')
  } else if (startsWith(i,'trana2s')) {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd2', 'd4', 'd6', 'mESC'))
    levelss <- c('d0', 'd2', 'd4', 'd6', 'mESC')
  } else if (startsWith(i,'tranfbs')){
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC'))
    levels <- c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC')
  } else if (startsWith(i,'treutlein')){
    cds@colData$Timepoint <- factor(cds@colData$time_point, levels = c('0', '2', '5', '20', '22'))
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="0"] <- "d0"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="2"] <- "d2"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="5"] <- "d5"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="20"] <- "d20"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="22"] <- "d22"
    levelss <- c('d0', 'd2', 'd5', 'd20', 'd22')
  }
  
  
  ggplot(agmat2, aes(x = name, y = value, color = moduleGroup, group = moduleGroup)) + geom_line() + 
    #  gghighlight(X %in% paste('Module',brig_std_g[[i]])) + 
    geom_hline(aes(yintercept = 0), color = 'gray') + theme_minimal() + facet_wrap(~moduleGroup) + 
    xlab('Timepoint') + ylab('Aggregated Gene Counts')
}


setwd("~/Documents/thesispaper")

library(ggplot2)
library(monocle3)

### PLOT GENE MODULE HEATMAPS
for (i in list.files('./final_files/traj_rds_final/')){
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  
  gt_result <- read.csv(paste0('./graph_test_outputs/', setname, '_graphtest_4aug22.csv'), row.names = 1)
  cds <- preprocess_cds(cds, num_dim = 50)
  
  # reorganise timepoints
  if (setname == 'brigdir') { 
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd5', 'd11'))
    levelss <- c('d0', 'd5', 'd11')
  } else if (setname == 'brigstd') {
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="d5"] <- "d6"
    levelss <- c('d0', 'd6', 'd12')
  } else if (setname == 'trana2s') {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint)
    levelss <- c('d0', 'd2', 'd4', 'd6', 'mESC')
  } else if (setname == 'tranfbs') {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC'))
    levelss <- c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC')
  } else if (setname == 'treutlein') {
    cds@colData$Timepoint <- factor(cds@colData$time_point, levels = c('0', '2', '5', '20', '22'))
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="0"] <- "d0"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="2"] <- "d2"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="5"] <- "d5"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="20"] <- "d20"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="22"] <- "d22"
    levelss <- c('d0', 'd2', 'd5', 'd20', 'd22')
  }
  
  traj_genes <- row.names(subset(gt_result, q_value < 0.05))
  
  gene_module_df <- find_gene_modules(cds[traj_genes,], resolution=c(10^seq(-6,-1)))
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=colData(cds)$Timepoint)
  
  agg_mat <- read.csv(paste0('./aggmat/', setname, '_aggmat_16jan23.csv'), row.names = 1)

  hm <- pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F)
  moduleGroup <- as.factor(cutree(hm$tree_row, k=4))
  gene_module_df$module4 <- moduleGroup[as.character(gene_module_df$module)]
  moduleGroup <- as.data.frame(moduleGroup)
  
  pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", cluster_cols = F, annotation_row = moduleGroup, fontsize = 7,
                     file = paste0('./heatmaps/',setname,'_heatmap.png'))
  
  write.csv(gene_module_df,paste0("./genemodule/",setname,"_genemodule_30jan23.csv"))
  
  agmat2 <- cbind(agg_mat, moduleGroup)
  agmat2 <- tibble::rownames_to_column(agmat2, 'module')
  agmat2 <- tidyr::pivot_longer(agmat2, cols = levelss)
  agmat2$name <- factor(agmat2$name, levels = levelss)
  
  ggplot(agmat2, aes(x = name, y = value, color = moduleGroup, group = module)) + geom_line() + 
    #  gghighlight(X %in% paste('Module',brig_std_g[[i]])) + 
    geom_hline(aes(yintercept = 0), color = 'gray') + theme_minimal() + facet_wrap(~moduleGroup) + 
    xlab('Timepoint') + ylab('Aggregated Gene Counts') + 
    scale_color_manual(values=c("#96CA00", "#FF9289", "#00DAE0", "#E199FF"))
  ggsave(paste0('./genemodule/', setname, '_profile_30jan23.png'))
}


setwd("~/Documents/thesispaper")

library(ggplot2)
library(patchwork)

### SAVE UMAPS
for (i in list.files('./final_files/traj_rds_final/')){
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))

  # reorganise timepoints
  if (setname == 'brigdir') { 
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd5', 'd11'))
    levelss <- c('d0', 'd5', 'd11')
  } else if (setname == 'brigstd') {
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="d5"] <- "d6"
    levelss <- c('d0', 'd6', 'd12')
  } else if (setname == 'trana2s') {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint)
    levelss <- c('d0', 'd2', 'd4', 'd6', 'mESC')
  } else if (setname == 'tranfbs') {
    cds@colData$Timepoint <- factor(cds@colData$Timepoint, levels = c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC'))
    levelss <- c('d0', 'd3', 'd6', 'd9', 'd12', 'mESC')
  } else if (setname == 'treutlein') {
    cds@colData$Timepoint <- factor(cds@colData$time_point, levels = c('0', '2', '5', '20', '22'))
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="0"] <- "d0"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="2"] <- "d2"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="5"] <- "d5"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="20"] <- "d20"
    levels(cds@colData$Timepoint)[levels(cds@colData$Timepoint)=="22"] <- "d22"
    levelss <- c('d0', 'd2', 'd5', 'd20', 'd22')
  }
  
  # load gene modules
  gene_module_df <- read.csv(paste0("./genemodule/",setname,"_genemodule_30jan23.csv"), row.names = 1)
  
  rowData(cds)$gene_name <- rownames(cds)
  rowData(cds)$gene_short_name <- rowData(cds)$gene_name
  
  plot_cells(cds, color_cells_by = 'Timepoint', label_cell_groups = FALSE)
  ggsave(paste0('./plots/',setname,'_umap_2feb23.png'), 
         width = 17.52, height = 18.92, units = 'cm', scale = 0.8)
  plot_cells(cds,
             genes=gene_module_df %>% filter(module4 %in% c(1,2,3,4)) %>% select(id, module4),
             label_cell_groups=FALSE, label_branch_points = F,
             label_roots = F, label_leaves = F) +
    scale_colour_gradient(low = "yellow", high = "red", na.value = NA)
  ggsave(paste0('./plots/',setname,'_umap_modulegroups_2feb23.png'), 
         width = 17.52, height = 18.92, units = 'cm', scale = 0.8)
  
}

library(figpatch)
library(ggplot2)
library(patchwork)
library(grid)

### COMBINE PLOTS INTO A SINGLE FIGURE
for (i in list.files('./final_files/traj_rds_final/')){
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  p1 <- fig(paste0("./plots/",setname,"_umap_2feb23.png"))
  p2 <- fig(paste0("./plots/",setname,"_umap_modulegroups_2feb23.png"))
  p3 <- fig(paste0("./heatmaps/",setname,"_heatmap.png"))
  p4 <- fig(paste0("./genemodule/",setname,"_profile_30jan23.png"))
  
  wrap_plots(p1, p2, p3, p4) + plot_annotation(tag_levels = "A")
  grid.rect(gp = gpar(lwd = 3, col = "black", fill = NA))
  
  ggsave(paste0("./combined_plots/",setname,"_comb.png"))
}



setwd("~/Documents/thesispaper")

library(org.Mm.eg.db)
library(clusterProfiler)
library(patchwork)
library(monocle3)
library(tidyverse)

### RUN GO/KEGG/WP ENRICHMENT ANALYSES 
for (i in list.files('./final_files/traj_rds_final/')){
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  
  gene_module_df <- read.csv(paste0("./genemodule/",setname,"_genemodule_30jan23.csv"), row.names = 1)
  
  for (a in c('GO','KEGG','WP')){
    print(paste0('...',a, '... '))
    for (k in unique(gene_module_df$module4)){
      print(paste('Module Group',k,'started!'))
      genes <- gene_module_df[gene_module_df$module4 == k,]$id
      if (a == 'GO'){
        ego <- enrichGO(gene        = genes,
                        universe      = rownames(cds),
                        OrgDb         = org.Mm.eg.db,
                        keyType       = "SYMBOL", 
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
        if (dim(summary(ego))[1] == 0){
          assign(paste0('p',k),ggplot() + annotate("text", x = 1, y = 1, size = 5, label = paste('No significant\nGO term found for\nmodule group',k)) + theme_void())
        } else { assign(paste0('p',k),
                        dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                          scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
      }
      else if (a == 'KEGG'){
        genes <- gene_module_df[gene_module_df$module4 == k,]$id
        genes <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
        univ <- mapIds(org.Mm.eg.db, keys = rownames(cds), column = "ENTREZID", keytype = "SYMBOL")
        ego <- enrichKEGG(gene          = genes,
                          universe      = univ,
                          organism      = "mmu",
                          keyType       = "kegg",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05)
        if (dim(summary(ego))[1] == 0){
          assign(paste0('p',k),ggplot() + annotate("text", x = 1, y = 1, size = 5, label = paste('No significant\nKEGG term found for\nmodule group',k)) + theme_void())
        } else { assign(paste0('p',k),
                        dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                          scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
      }
      else if (a == 'WP'){
        genes <- gene_module_df[gene_module_df$module4 == k,]$id
        genes <- mapIds(org.Mm.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL")
        ego <- enrichWP(genes, organism = "Mus musculus")
        if (dim(summary(ego))[1] == 0){
          assign(paste0('p',k),ggplot() + annotate("text", x = 1, y = 1, size = 5, label = paste('No significant\nWikiPathways term found for\nmodule group',k)) + theme_void())
        } else { assign(paste0('p',k),
                        dotplot(ego, showCategory=10, font.size = 10, title = paste('Module Group',k)) + 
                          scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
      }
    }
    pp <- p1 + p2 + p3 + p4 + plot_layout(widths = c(2, 2))
    ggsave(paste0('./clusterProfiler_plots/',paste(setname, a, '2feb23', sep = '_'),'_dotplots.png'), scale = 1.5)
    print(paste(setname,'finished!'))
  }
}



#### CALCULATE SPEARMAN CORRELATION
for (i in list.files('./final_files/traj_rds_final/')){
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  
  cds2 <- choose_graph_segments(cds, clear_cds = FALSE)
  cells <- colnames(cds2)
  pseud <- sort(pseudotime(cds)[cells])
  
  counts_cells <- normalized_counts(cds)[,names(pseud)]
  
  cors <- c()
  counter <- 0
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(rownames(cds)), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = " s=")   # Character used to create the bar
  
  for (j in rownames(cds)){
    cors <- append(cors, cor(counts_cells[j,], pseudotime(cds)[cells], method = 'spearman'))
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
  }
  
  cors_genes <- data.frame(rownames(cds), cors)
  colnames(cors_genes) <- c('genes', 'spearman_corr')
  
  write.csv(cors_genes,paste0("./correlation/",setname,"_spearman_10feb23.csv"))
  
  close(pb) # Close the connection
  
}



library(dorothea)
library(dplyr)

## SUBSET 10 MIN/MAX TFS ACCORDING TO SPEARMAN CORRELATIONS

# find the TFs from CellDataSet objects
tfs <- dorothea_mm %>% filter(confidence %in% c('A','B','C')) %>%
  filter(tf %in% rownames(cds))


for (i in list.files('./final_files/traj_rds_final/')){
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  cors_genes <- read.csv(paste0("./correlation/",setname,"_spearman_10feb23.csv"), 
                         row.names = 1)
  
  cors <- bind_rows(cors_genes %>% slice_min(spearman_corr, n = 10),
                    cors_genes %>% slice_max(spearman_corr, n = 10)) %>%
    mutate(if_tf = if_else(genes %in% tfs$tf, TRUE, FALSE)) %>%
    mutate(if_target = if_else(genes %in% tfs$target, TRUE, FALSE))
  
   
  write.csv(cors,paste0("./correlation/",setname,"_tfs_10feb23.csv"))

}

 

## FIND THE TARGETS OF TFS AND CREATE TF-TARGET PAIR LISTS
for (i in list.files('./correlation/', 'tfs')){
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  cors <- read.csv(paste0("./correlation/",setname,"_tfs_10feb23.csv"), row.names = 1)
  
  tfgene <- cors %>% filter(if_tf == TRUE)
  
  colnames(tfgene)[1] <- 'tf'
  
  
  if (nrow(tfgene) >= 1){
    tfgene <- tfs %>% filter(tf %in% tfgene$tf) %>%
      left_join(tfgene, by = 'tf')
  }
  
  write.csv(tfgene,paste0("./final_files_nov22/correlation/",setname,"_tf-gene_12feb23.csv"))
  
}

## FIND THE TFS OF TARGETS AND CREATE TARGET-TF PAIR LISTS
for (i in list.files('./final_files_nov22/correlation/', 'tfs')){
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  cors <- read.csv(paste0("./final_files_nov22/correlation/",setname,"_tfs_10feb23.csv"), row.names = 1)
  
  tfgene <- cors %>% filter(if_target == TRUE)
  
  colnames(tfgene)[1] <- 'target'
  
  if (nrow(tfgene) >= 1){
    tfgene <- tfs %>% filter(target %in% tfgene$target) %>%
      left_join(tfgene, by = 'target')
  }
  
  write.csv(tfgene,paste0("./final_files_nov22/correlation/",setname,"_target-tf_12feb23.csv"))
  
}


##### CALCULATE SPEARMAN CORRELATION FOR SIDE BRANCHES OF TRAJECTORIES
for (i in list.files('./final_files/traj_rds_final/')){
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  cds <- readRDS(paste0("./final_files/traj_rds_final/",i))
  
  # manual selection of graph segments for each dataset
  cds2 <- choose_graph_segments(cds, clear_cds = FALSE)
  if (ncol(cds2) == 0){next}
  cells <- colnames(cds2)
  pseud <- sort(pseudotime(cds)[cells])
  
  counts_cells <- normalized_counts(cds)[,names(pseud)]
  
  cors <- c()
  counter <- 0
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(rownames(cds)), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = " s=")   # Character used to create the bar
  
  for (j in rownames(cds)){
    cors <- append(cors, cor(counts_cells[j,], pseudotime(cds)[cells], method = 'spearman'))
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
  }
  
  cors_genes <- data.frame(rownames(cds), cors)
  colnames(cors_genes) <- c('genes', 'spearman_corr')
  
  write.csv(cors_genes,paste0("./final_files_nov22/correlation/",setname,"_spearman_2ndd_10feb23.csv"))
  
  close(pb) # Close the connection
  
}


library(dorothea)
library(dplyr)

# find the TFs from CellDataSet objects
tfs <- dorothea_mm %>% filter(confidence %in% c('A','B','C')) %>%
  filter(tf %in% rownames(cds))


for (i in list.files('./final_files/traj_rds_final/')){
  ## SUBSET 10 MIN/MAX TFS ACCORDING TO SPEARMAN CORRELATIONS
  setname <- strsplit(i, "_")[[1]][1]
  
  print(setname)
  
  if (setname == 'brigstd'){next}
  cors_genes <- read.csv(paste0("./final_files_nov22/correlation/",setname,"_spearman_2ndd_10feb23.csv"), 
                         row.names = 1)
  
  cors <- bind_rows(cors_genes %>% slice_min(spearman_corr, n = 10),
                    cors_genes %>% slice_max(spearman_corr, n = 10)) %>%
    mutate(if_tf = if_else(genes %in% tfs$tf, TRUE, FALSE)) %>%
    mutate(if_target = if_else(genes %in% tfs$target, TRUE, FALSE))
  
  
  write.csv(cors,paste0("./final_files_nov22/correlation/",setname,"_tfs_2ndd_10feb23.csv"))
  
  ## FIND THE TARGETS OF TFS AND CREATE TF-TARGET PAIR LISTS
  tfgene <- cors %>% filter(if_tf == TRUE)
  
  colnames(tfgene)[1] <- 'tf'
  
  
  if (nrow(tfgene) >= 1){
    tfgene <- tfs %>% filter(tf %in% tfgene$tf) %>%
      left_join(tfgene, by = 'tf')
  }
  
  write.csv(tfgene,paste0("./final_files_nov22/correlation/",setname,"_tf-gene_2ndd_12feb23.csv"))
  
  ## FIND THE TFS OF TARGETS AND CREATE TARGET-TF PAIR LISTS
  tfgene <- cors %>% filter(if_target == TRUE)
  
  colnames(tfgene)[1] <- 'target'
  
  if (nrow(tfgene) >= 1){
    tfgene <- tfs %>% filter(target %in% tfgene$target) %>%
      left_join(tfgene, by = 'target')
  }
  
  write.csv(tfgene,paste0("./final_files_nov22/correlation/",setname,"_target-tf_2ndd_12feb23.csv"))
  
}




library(openxlsx)
wb <- createWorkbook()

tfs <- c()
for (i in list.files('final_files_nov22/correlation/', pattern = 'onlyTFs')){
  setname <- strsplit(i, "_")[[1]][1]
  if (grepl('2ndd', i)){setname <- paste0(setname,'_2ndd')}
  print(setname)
  
  cors_genes <- read.csv(paste0("./final_files_nov22/correlation/",i), row.names = 1)
  
  addWorksheet(wb, setname)
  writeData(wb, sheet = setname, cors_genes)
}

saveWorkbook(wb, "onlyTfsALL.xlsx", overwrite = TRUE)



tfs <- c()
for (i in list.files('final_files_nov22/correlation/', pattern = 'onlyTFs')){
  setname <- strsplit(i, "_")[[1]][1]
  if (grepl('2ndd', i)){next}
  print(setname)
  
  cors_genes <- read.csv(paste0("./final_files_nov22/correlation/",i), row.names = 1)
  
  tfs <- c(tfs,cors_genes$genes)
}


library(pheatmap)
j <- 1
for (i in list.files('final_files_nov22/correlation/', pattern = 'onlyTFs')){
  setname <- strsplit(i, "_")[[1]][1]
  if (grepl('2ndd', i)){next}
  print(setname)
  
  if (j>1){
    x <- read.csv(paste0("./final_files_nov22/correlation/",i), row.names = 1)
    x$dataset <- setname
    cors_genes <- rbind(cors_genes, x)
  } else {
    cors_genes <- read.csv(paste0("./final_files_nov22/correlation/",i), row.names = 1)
    cors_genes$dataset <- setname
  }
  j <- j+1
}

cors_genes2 <- data.frame(pivot_wider(cors_genes, names_from = dataset, values_from = spearman_corr))
rownames(cors_genes2) <- cors_genes2$genes
cors_genes2 <- cors_genes2[,2:ncol(cors_genes2)]

colnames(cors_genes2) <- gsub("..", "", colnames(cors_genes2))

pheatmap::pheatmap(cors_genes2, color=colorRampPalette(c("navy", "white", "red"))(50),
                   clustering_method="ward.D2", cluster_cols = F, cluster_rows = F,
                   file="heatmap_correlation_na0.png", 
                   width=8.27, height=14)

cors_genes2[is.na(cors_genes2)] <- 0
pheatmap::pheatmap(cors_genes2, color=colorRampPalette(c("navy", "white", "red"))(50),
                   clustering_method="ward.D2", cluster_cols = F, 
                   file="heatmap_correlation.png", 
                   width=8.27, height=11.67)



# -----

library(rentrez)
res <- entrez_search(db = "gene", term = "Pou5f1")
esums <- entrez_summary(db = "gene", id = res$ids)
for (i in esums){
  if ((i$name == gene) & (i$organism$scientificname == 'Mus musculus')){
    print(i$name)
    print(i$summary)
  }
}


gene_name <- c()
summary <- c()
for (g in rownames(cors_genes2)){
  res <- entrez_search(db = "gene", term = paste0(g, "[GENE]", " AND Mus musculus[ORGN]"))
  esums <- entrez_summary(db = "gene", id = res$ids)
  if ((length(esums) > 1) & (length(esums) < 21)){
    id <- names(which(extract_from_esummary(esums, 'name') == g))[1]
    gene_res <- get(id, esums)$name
    summary_to_added <- extract_from_esummary(get(id, esums), 'summary')
  } else{
    gene_res <- extract_from_esummary(esums, 'name')
    summary_to_added <- extract_from_esummary(esums, 'summary')
  }
  gene_name <- c(gene_name, gene_res)
  a <- length(summary)
  summary <- c(summary, str_split(summary_to_added, "\\. ")[[1]][1])
  b <- length(summary)
  if (b-a > 1){
    print(paste0('!!!! -> ',gene_name))
  }
  print(paste(gene_res, 'done!!'))
}

write_csv(data.frame(summary), 'summary.txt')

for (g in 1:nrow(cors_genes2)){
  rownames(cors_genes2)[g] <- paste(rownames(cors_genes2)[g], summary[g], sep=": ")
}

rownames(cors_genes2) <- genes


# -----

for (setname in colnames(cors_genes2)){
  genes_neg <- rownames(cors_genes2[which(cors_genes2[[setname]] < 0),])
  genes_pos <- rownames(cors_genes2[which(cors_genes2[[setname]] > 0),])
  genes <- list('Negatively Correlated TFs' = genes_neg,
                'Positively Correlated TFs' = genes_pos)
  for (k in names(genes)){
  ego <- enrichGO(gene        = genes[[k]],
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "SYMBOL", 
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  if (dim(as.data.frame(ego))[1] == 0){
    assign(paste0('p',which(names(genes) == k)),ggplot() + annotate("text", x = 1, y = 1, size = 5, label = paste('No significant\nGO term found for\nmodule group',k)) + theme_void())
  } else { assign(paste0('p',which(names(genes) == k)),
                  dotplot(ego, showCategory=10, font.size = 10, title = k) + 
                    scale_y_discrete(labels=function(x) str_wrap(x, width=45)))}
  }
  pp <- p1 / p2
  ggsave(paste0('./final_files_nov22/tf_go/',paste(setname, '30apr23','dotplots.png', sep = '_')), scale = 1.5)
  print(paste(setname,'finished!'))
}



setwd("~/Documents/thesispaper")
library(ggplot2)
treu <- readRDS("~/Documents/thesispaper/final_files/traj_rds_final/treutlein_monocle3v1_5nov21.rds")

cds <- treu[which(rownames(treu) %in% c('Pou5f1', 'Sox2', 'Klf4', 'Myc', 'Nanog')),]
cds
plot_genes_in_pseudotime(cds, label_by_short_name = FALSE, color_cells_by="time_point")
ggplot2::ggsave('final_files_nov22/treutlein_trends/whole_trends.png')
cds1 <- choose_graph_segments(treu)
cells1 <- colnames(cds1)
cds1 <- choose_graph_segments(treu)
cells2 <- colnames(cds1)
cds1 <- treu[which(rownames(treu) %in% c('Pou5f1', 'Sox2', 'Klf4', 'Myc', 'Nanog')),cells1]
cds2 <- treu[which(rownames(treu) %in% c('Pou5f1', 'Sox2', 'Klf4', 'Myc', 'Nanog')),cells2]
plot_genes_in_pseudotime(cds1, label_by_short_name = FALSE, color_cells_by="time_point")
ggplot2::ggsave('final_files_nov22/treutlein_trends/upper_branch_trends.png')
plot_genes_in_pseudotime(cds2, label_by_short_name = FALSE, color_cells_by="time_point")
ggplot2::ggsave('final_files_nov22/treutlein_trends/lower_branch_trends.png')

