# Reprogramming and Differentiation Pseudotime Analysis

This repository contains scripts to find the gene modules, create dotplots and heatmaps, calculate correlation with respect to pseudotime and create heatmap, and generate other figures for the preprint titled "Comparative Pseudotime Analysis of Single-Cell Trajectories in Reprogramming and Neural Differentiation"

1. graphtest.R  
for each dataset, run Monocle 3’s graph_test → select significant genes → find gene modules → plot and save heatmap and tables.

2. genemoduleorg.R  
exploratory organization of modules, trend plots, GO analysis examples.

3. manuscriptfigures.R  
end‑to‑end generation of all manuscript figures: GO/KEGG/WP enrichment, heatmaps, profiles, UMAPs, correlations, TF analyses, Entrez summaries, combined figures.
