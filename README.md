# Immune cell type signature discovery and random forest classification for analysis of single cell gene expression datasets 
Aybey et al.

This git repository contains all the scripts to generate the figures and tables in our preprint.

1. preprocess_scripts
This folder includes preprocessing scripts to generate RDS files for single cell experiments.

2. gene_lists_papers
This file contains published signatures we used in the publication.

3. immune_signature_discovery_workflow
- S1_discovery_dataset_integration.R: Integration of 4 discovery datasets
- S2_immune_cell_type_signature_discovery_workflow.R: Discovery of immune cell type signatures
  - Outputs:
    - dbscan_cluster_genes_param_0.txt: intermediate result of all signature genes 
    - **dbscan_cluster_genes_param_0_GO_annot.txt: final annotated clusters**
    - Fig2e_mean_score_our_signatures_heatmap.svg
    - Fig2c_filtered_genes_umap.svg
    - Fig2d_umap_refined_gene_clusters.svg
    - Tab2_final_gene_list.txt: Tab. 2; summarized version of our gene signatures
    - SuppFig1_db_GO_analysis.pdf
- S3_literature_comparison_excel.R: Comparison of our gene signatures in the original datasets
  - Output: dbscan_cluster_genes_lit_annot_updated.tsv
- S4_signature_comparison_literature_jaccard_scores.R: Calculating Jaccard index scores and overlap scores
  - Output: SuppFig2_overlap_and_jaccard.svg

4. benchmarking
  - Benchmarking results for two independent PBMC datasets (Zheng sorted and Kotliarov).
  - For each dataset we ran random forest approach (RF) either using our selected immune cell type signatures or other published signatures --> saved under RF_gene or RF_hao
  - Additionally we included singleR and Seurat predictions. For Kotliarov we tested three reference datasets. For Zheng sorted dataset we have used only Hao reference.
  - Each R file follows similar rule:"(method e.g.rf, singler, seurat, rf_ourgenes, rf_nieto ...)\_(reference dataset e.g. hao)\_(benchmarking dataset e.g. kotliarov).R"
  - Prediction scores are saved as txt files under the corresponding files
  - celltype\_annot\_stats.R: calculates cell type averaged prediction scores
  - kotliarov\_harmonize\_prediction\_scores\_func.R: harmonizes prediction results and Kotliarov cell type annotations and calculates cell type averaged prediction scores
  
**Contact**



