# Load the WLS functions
source("Rcode/wls_functions.R")

dan_predTF <- function(scores, labels, data_fit) {
  res_varRep <- varRep(data_fit[,labels.tf1_x], data_fit[,labels.tf1_y], plot = F)
  return(class_predTF(data_fit[,labels.tf1_x],data_fit[,labels.tf2_y],scores.x, scores.y,t=res_varRep$t,sigma=res_varRep$sigma,level=0.99,plot=T))
}
# Set some scores
# For ELK1 vs ETS1

scores.x = c(0.292, 0.297, 0.161, 0.527)
scores.y = c(0.313, 0.272, 0.063, 0.661)

data_fit = read.table("PBMdata/Combined_Ets1_100nM_Elk1_100nM_50nM_Gabpa_100nM_Bound_filtered_normalized.txt",col.names=c("seq","ets1_100nM","elk1_100nM","elk1_50nM","gabpa_100nM"))
# Labels
labels.tf1_x <- "elk1_100nM"  ## Used in Replicate estimate (VarRep) and TF pair estimate (ClassTF)
labels.tf1_y <- "elk1_50nM"   ## Used in Replicate estimate
labels.tf2_y <- "ets1_100nM"  ## Used in TF pair estimate
result = dan_predTF(scores, labels, data_fit)
print(result)
