########### Recall and precision ###########################
#cell type names from the original(ground truth) 
Cell_type_stats=function(DF_compare) {
cell_names= sort(unique(DF_compare$original))
tab_conf= table(pred=DF_compare$predicted,orig=DF_compare$original)
cell_names=intersect(cell_names,rownames(tab_conf))

#Classification statistics
TP=diag(tab_conf[cell_names,cell_names])
FP=rowSums(tab_conf)[cell_names]-TP
FN=colSums(tab_conf)[cell_names]-TP
TN=sum(tab_conf)-FP-FN-TP

recall=100*TP/(TP+FN)
precision=100*TP/(TP+FP)
npv=100*TN/(TN+FN)
spec=100*TN/(TN+FP)
accuracy=100*(TP+TN)/(TP+TN+FP+FN)

weight_factor= colSums(tab_conf)[cell_names]
# Weighted average recall(sensitivity)
df=data.frame(Sensitivity=round(sum(recall*weight_factor)/sum(weight_factor),digits = 2),
# Weighted average precision
PPV=round(sum(precision*weight_factor)/sum(weight_factor),digits = 2),
# NPV 
NPV=round(sum(npv*weight_factor)/sum(weight_factor),digits = 2),
# specificity 
Specificity=round(sum(spec*weight_factor)/sum(weight_factor),digits = 2),
Accuracy= round(sum(accuracy*weight_factor)/sum(weight_factor),digits = 2)
)
df$f1_score=round(2*(df$PPV*df$Sensitivity)/(df$PPV+df$Sensitivity),digits=2)
df

}

