# Load the WLS functions
source("Rcode/wls_functions.R")

dan_predTF <- function(scores1, scores2, lab, data_fit) {
  # scores1/scores2 must be same length
  # This is so that class_predTF can slice off the first 3 columns
  M1 = matrix(0, nrow = length(scores1), ncol = length(scores1) + 3)
  M2 = matrix(0, nrow = length(scores2), ncol = length(scores2) + 3)
  for(i in 1:length(scores1)) {
    M1[i,i + 3] = scores1[i]
    M2[i,i + 3] = scores2[i]
  }
  res_varRep <- varRep(data_fit[,lab.tf1_x], data_fit[,lab.tf1_y], plot = F)
  return(class_predTF(data_fit[,lab.tf1_x],data_fit[,lab.tf2_y],M1,M2,t=res_varRep$t,sigma=res_varRep$sigma,level=0.99,plot=T))
}
# Set some scores
# For ELK1 vs ETS1

Elk1Scores = c(0.292, 0.297, 0.161, 0.527)
Ets1Scores = c(0.313, 0.272, 0.063, 0.661)

data_fit = read.table("PBMdata/Combined_Ets1_100nM_Elk1_100nM_50nM_Gabpa_100nM_Bound_filtered_normalized.txt",col.names=c("seq","ets1_100nM","elk1_100nM","elk1_50nM","gabpa_100nM"))
# Labels
lab.tf1_x <- "elk1_100nM"  ## Used in Replicate estimate (VarRep) and TF pair estimate (ClassTF)
lab.tf1_y <- "elk1_50nM"   ## Used in Replicate estimate
lab.tf2_y <- "ets1_100nM"  ## Used in TF pair estimate
result = dan_predTF(Elk1Scores, Ets1Scores, lab, data_fit)
print(result)
