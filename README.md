# Immune Cell Type Signature Discovery and Random Forest Classification for Analysis of Single-Cell Gene Expression Datasets

[doi.org/10.1101/2023.03.24.534078](https://doi.org/10.1101/2023.03.24.534078)

### Authors: 
Bogac Aybey, Sheng Zhao, Benedikt Brors, and Eike Staub

This Git repository contains all the scripts to generate the figures and tables in our preprint. The main analysis is divided into two main parts/folders.

## Data
Discovery TME datasets (/TME_datasets), benchmarking datasets (Kotliarov and Zheng), reference (Hao) and Kartha IFN PBMC datasets are stored. Preprocessing scripts are in their corresponding folders. 

- Benchmarking_dataprep_h5da_celltypist.R: Hao reference, Kotliarov, Zheng and Kartha datasets are converted to h5ad format for CellTypist workflow.

Kartha IFN PBMC dataset 

- Integration_Kartha_GSE178429.R: Kartha dataset has been integrated using Seurat based on each sample (IFN-g treated and control samples across different time points). The integrated matrix was provided for cell type classification in scType and Seurat.

## Discovery

- Discovery_integration.R: Integration of our 7 TME discovery datasets
- Discovery_discovery.R: Discovery workflow
- Discovery_validation_TMEatlasnieto.R: Validating our immune cell type signatures in TME immune cell atlas from Nieto et al.
- Discovery_signature_comparison.R: Comparison of our immune cell type signatures with other published signatures

## RF classifier

You can run RF classifier using appropriate reference dataset, cell type signature genes and query dataset using Example_RF_classifier.R. 
The output of this example is rf_Example_Hao_Zheng_ourgenes_163.csv. RF model is also saved into the directory designated as save_dir.

## Benchmarking

Scripts stored to run different cell type prediction methods (CellTypist, random forest RF, scType, sinleR, CHETAH, and Seurat) in benchmarking datasets.

- Run RF models: Benchmarking_RF_commands_kotliarov.R, Benchmarking_RF_commands_zheng.R
- Run CellTypist: Benchmarking_celltypist_commands.R
- Run other models: Benchmarking_classifiers_commands_kotliarov.R, Benchmarking_classifiers_commands_zheng.R
- Benchmarking_plots.R: Create benchmarking plots comparing RF with other methods (run in default methods) and also RF using our immune cell type genes against other published repertoires.

#### Models

Models from CellTypist and random forest are saved.

#### Prediction Results

All cell type prediction results for CellTypist, random forest, scType, sinleR, CHETAH, and Seurat are saved as CSV files in the corresponding folders for each method. 

- Benchmarking datasets: Kotliarov and Zheng scRNA-seq
- IFN dataset: Kartha
- Unified file name: method_reference_queery_genechoice_numberofgenes.csv (e.g., rf_Hao_Kartha_ourgenes_164.csv)

#### Diffexp_kartha - Downstream analysis bias in IFN PBMC Kartha dataset

Compare DCs labeled by random forest vs scType, CellTypist and Seurat and show downstream analysis bias based on IFN-g effect on DC 

-  kartha_celltyping.R : Cell type labeling by random forest vs scType and Seurat.
-  *CellTypist labels were called using Benchmarking_celltypist_commands.R
-  kartha_DC_mono_classification_analysis.R : compare DCs labeled by RF (DC_RF) vs misclassified monocytes (Mono_RF-classified as monocytes in RF  but DCs in other methods)
    -  Compare expression of DC or monocyte gene markers
    -  Compare IFNg Hallmark gene expression scores 
-  kartha_diff_exp_analysis_DC.R
    -  Calculate p-values for control vs IFNg treatment at each time point focusing only on DCs - based on the cell type classifications by each method
    -  Plot p-value-to-p-value scatter plots for RF vs other methods


## Functions

- discovery_functions.R: Functions used in discovery workflow
- benchmarking_RF_training.R: Functions used in training and testing random forest classifiers
- benchmarking_classifiers.R: Functions used in singleR, Seurat, CHETAH, and scType cell type classification
- benchmarking_celltypist.py: Functions to classify in CellTypist
- benchmarking_prediction_metrics.R: Functions to calculate cell type prediction metrics
- RF_celltype_classifier.R: RF classifier function

  
## Gene lists from published signatures

Gene lists used in this study are stored in /gene_lists_papers.

## Contact

Bogac Aybey: [bogac.aybey@stud.uni-heidelberg.de](mailto:bogac.aybey@stud.uni-heidelberg.de)
