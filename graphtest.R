# "Graph Tests and Module Discovery for Thesis Datasets"
#
# This R script applies the steps below for each dataset:
# * loads RDS files including the trajectory
# * applies graph test
# * finds gene modules of DEGs
# * does GO analysis to the modules

# These steps are applied to the datasets below:
# * Tran-A2S
# * Tran-FBS
# * Treutlein
# * Briggs-standard
# * Briggs-direct



# Tran-A2S
# load the dataset
library(monocle3)
cds <- readRDS("./final_files/traj_rds_final/tran_a2s_noirr_monocle3v1_5nov21.rds")
unique(cds@colData$Timepoint)
# preprocesing to prevent find_gene_modules error
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, resolution=1e-5)
# run the graph test
trana_traj <- graph_test(cds, neighbor_graph="principal_graph", cores = 4, alternative = 'less')
trana_traj_genes <- row.names(subset(trana_traj, q_value < 0.05))
# find the gene modules
gene_module_df <- find_gene_modules(cds[trana_traj_genes,])
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$Timepoint)
# calculate aggregated gene expression and organise the data frames
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
agg_mat <- agg_mat[,c('d0','d2','d4','d6','mESC')]
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", 
                   cluster_cols = F, filename="tran_a2s_genemodules_8feb22.png")
# save all the outputs
write.csv(gene_module_df, "tran_a2s_geneModules_8feb22.csv", row.names = F)
write.csv(as.data.frame(as.matrix(agg_mat)), "tran_a2s_aggmat_8feb22.csv")




# Tran-FBS
# load the dataset
cds <- readRDS("./final_files/traj_rds_final/tran_fbs_noirr_monocle3v1_5nov21.rds")
# preprocesing to prevent find_gene_modules error
cds <- preprocess_cds(cds, num_dim = 50)
#cds <- reduce_dimension(cds)
#cds <- cluster_cells(cds, resolution=1e-5)
# run the graph test
tranf_traj <- graph_test(cds, neighbor_graph="principal_graph", cores = 4, alternative = 'two.sided')
tranf_traj_genes <- row.names(subset(tranf_traj, q_value < 0.05))
# find the gene modules
gene_module_df <- find_gene_modules(cds[tranf_traj_genes,])
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$Timepoint)
# calculate aggregated gene expression and organise the data frames
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
agg_mat <- agg_mat[,c('d0','d3','d6','d9','d12','mESC')]
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", 
                   cluster_cols = F, filename="tran_fbs_genemodules_8feb22.png")
# save all the outputs
write.csv(gene_module_df, "tran_fbs_geneModules_8feb22.csv", row.names = F)
write.csv(as.data.frame(as.matrix(agg_mat)), "tran_fbs_aggmat_8feb22.csv")




# Treutlein:
# load the dataset
cds <- readRDS("./final_files/traj_rds_final/treutlein_monocle3-final.rds")
# preprocesing to prevent find_gene_modules error
#cds <- preprocess_cds(cds, num_dim = 50)
#cds <- reduce_dimension(cds)
#cds <- cluster_cells(cds, resolution=1e-5)
# run the graph test
treut_traj <- graph_test(cds, neighbor_graph="principal_graph", alternative = 'two.sided')
treut_traj_genes <- row.names(subset(treut_traj, q_value < 0.05))
# find the gene modules
gene_module_df <- find_gene_modules(cds[treut_traj_genes,])
cds@colData$time_point <- paste0('d',cds@colData$time_point)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$time_point)
# calculate aggregated gene expression and organise the data frames                        
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
agg_mat <- agg_mat[,c('d0','d2','d5','d20','d22')]
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", 
                   cluster_cols = F, filename="treutlein_genemodules_8feb22.png")
# save all the outputs
write.csv(gene_module_df, "treutlein_geneModules_8feb22.csv", row.names = F)
write.csv(as.data.frame(as.matrix(agg_mat)), "treutlein_aggmat_8feb22.csv")




# Briggs-standard
# load the dataset
cds <- readRDS("./final_files/traj_rds_final/briggs_std_monocle3v1_4nov21.rds")
# preprocesing to prevent find_gene_modules error
# cds <- preprocess_cds(cds, num_dim = 50)
# cds <- reduce_dimension(cds)
# cds <- cluster_cells(cds, resolution=1e-5)
# run the graph test
brist_traj <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
brist_traj_genes <- row.names(subset(brist_traj, q_value < 0.05))
# find the gene modules
gene_module_df <- find_gene_modules(cds[brist_traj_genes,])
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$Timepoint)
# calculate aggregated gene expression and organise the data frames                        
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
agg_mat <- agg_mat[,c('d0','d5','d12')]
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", 
                   cluster_cols = F, filename="briggs_standard_genemodules_8feb22.png")
# save all the outputs
write.csv(gene_module_df, "briggs_std_geneModules_8feb22.csv", row.names = F)
write.csv(as.data.frame(as.matrix(agg_mat)), "briggs_std_aggmat_8feb22.csv")


# Briggs-direct
# load the dataset
cds <- readRDS("./final_files/traj_rds_final/briggs_direct_monocle3v1_allCells_org_4nov21.rds")
# preprocesing to prevent find_gene_modules error
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, resolution=1e-5)
levels(cds@colData$Timepoint) <- c("d0", "d4", "d11")
# run the graph test
bridi_traj <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
bridi_traj_genes <- row.names(subset(bridi_traj, q_value < 0.05))
# find the gene modules
gene_module_df <- find_gene_modules(cds[bridi_traj_genes,])#, resolution=c(10^seq(-6,-1)))
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$Timepoint)
# calculate aggregated gene expression and organise the data frames                        
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
agg_mat <- agg_mat[,c('d0','d4','d11')]
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, clustering_method="ward.D2", 
                   cluster_cols = F, filename="briggs_direct_genemodules_8feb22.png")
# save all the outputs
write.csv(gene_module_df, "briggs_direct_geneModules_8feb22.csv", row.names = F)
write.csv(as.data.frame(as.matrix(agg_mat)), "briggs_direct_aggmat_8feb22.csv")