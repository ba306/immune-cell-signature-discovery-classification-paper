import sys
import scanpy as sc
import celltypist
from celltypist import models

#Running cell typist using different parameters
# Either ImmuneAllLow as model
# Or train model using input genes e.g. HVGs with a given reference
# Or train model using default with a given reference
# Saves the model to ./benchmarking/models/
# Saves the predictions to ./benchmarking/prediction_results/celltypist

ref_name=sys.argv[1] # how to name the model "hao"
reference=sys.argv[2] # name of the reference data with full directory './data/Hao/pbmc_hao_ref_up.h5ad'
ref_col=sys.argv[3] # which columns are cell type metadata "harmonized_celltype"

query_name=sys.argv[4] # how to name the predicted query data "zheng"
query=sc.read(sys.argv[5]) # # name of the query data with full directory "./data/Zheng/pbmc_zheng_sorted_2k.h5ad"

model_name=sys.argv[6] # if Immune_All_Low use immune_all_low model from the package if not it will print this name as model name 
genes = sys.argv[7] # list of genes to use c("X,Y,Z")
genes = genes.split(',') 
geneset=sys.argv[8] # which genesets are used "HVGs" or gene signatures e.g. "Our_genes"

# Prediction by cell typist model
# using the default module immune_all_low from the immune reference
if model_name=="Immune_All_Low":
  model = models.Model.load(model='Immune_All_Low.pkl')
  print("Immune_all_low classifying")
  predictions = celltypist.annotate(query, model = model, majority_voting = True)
  predictions.predicted_labels.rename(columns={'majority_voting': 'Predictions'}, inplace=True)
  predictions.predicted_labels.to_csv("./benchmarking/prediction_results/celltypist/celltypist_"+ref_name+"_"+query_name+"_"+"ImmuneAllLow"+"_"+str(len(model.features))+".csv")
# training using other gene sets or HVGs
elif model_name!="Immune_All_Low" and len(genes)>1:
  reference=sc.read(reference) #'./data/pbmc_hao_ref_up.h5ad'
  ref_selected=reference[:,genes]
  no_genes=len(genes)
  model = celltypist.train(ref_selected, labels = ref_col, n_jobs = 10, feature_selection = False,check_expression=False)
  print("Training done")
  model.write("./benchmarking/models/model_from_"+ref_name+"_"+geneset+"_"+str(no_genes)+".pkl")
  print("Model saved")
  predictions = celltypist.annotate(query, model = model, majority_voting = True)
  predictions.predicted_labels.rename(columns={'majority_voting': 'Predictions'}, inplace=True)
  predictions.predicted_labels.to_csv("./benchmarking/prediction_results/celltypist/celltypist_"+ref_name+"_"+query_name+"_"+geneset+"_"+str(no_genes)+".csv")

# training using default top_gnees=300
else:
  model = celltypist.train(reference, labels = ref_col, n_jobs = 10, feature_selection = True)
  print("Training done")
  model.write("./benchmarking/models/model_from_"+ref_name+"_"+geneset+"_"+str(len(model.features))+".pkl")
  print("Model saved")
  predictions = celltypist.annotate(query, model = model, majority_voting = True)
  predictions.predicted_labels.rename(columns={'majority_voting': 'Predictions'}, inplace=True)
  predictions.predicted_labels.to_csv("./benchmarking/prediction_results/celltypist/celltypist_"+ref_name+"_"+query_name+"_"+geneset+"_"+str(len(model.features))+".csv")
