library(nlme)

input.path <- "../PBMdata"
setwd(input.path)

################################################################
# Part 1. Identify differentially preferred sites for PBM data
################################################################

# ! Orders of tf1_x, tf1_y, tf2_y are very important, do not change!

###begin hand typing input###
# ETS family, between Elk1 and Ets1
data_ets = read.table("Combined_Ets1_100nM_Elk1_100nM_50nM_Gabpa_100nM_Bound_filtered_normalized.txt",col.names=c("seq","ets1_100nM","elk1_100nM","elk1_50nM","gabpa_100nM"))
data_fit = data_ets
lab.tf1_x <- "elk1_100nM"  ## Used in Replicate estimate (VarRep) and TF pair estimate (ClassTF)
lab.tf1_y <- "elk1_50nM"   ## Used in Replicate estimate
lab.tf2_y <- "ets1_100nM"  ## Used in TF pair estimate

# E2F family, between E2f1 and E2f4
data_e2f = read.table("Combined_E2f1_250nM_200nM_E2f3_250nM_E2f4_500nM_800nM_Bound_filtered_normalized.txt",col.names = c("seq","e2f1_250nM","e2f1_200nM","e2f3_250nM","e2f4_500nM","e2f4_800nM"))
data_fit = data_e2f
lab.tf1_x = "e2f1_250nM"  ## Used in Replicate estimate (VarRep) and TF pair estimate (ClassTF)
lab.tf1_y = "e2f1_200nM"  ## Used in Replicate estimate
lab.tf2_y = "e2f4_500nM"  ## Used in TF pair estimate

# bHLH family, between Mad and Myc
data_bHLH = read.table("Combined_Max_Myc_Mad_MadL_Bound_filtered_normalized.txt",col.names = c("seq","max","myc","mad","madL"))
data_fit = data_bHLH
lab.tf1_x = "mad"
lab.tf1_y = "madL"
lab.tf2_y = "myc"

# Runx family, between Runx1 and Runx2
data_runx = read.table("Combined_Runx1_10nM_50nM_Runx2_10nM_Bound_filtered_normalized.txt",col.names = c("seq","runx1_10nM","runx1_50nM","runx2_10nM"))
data_fit = data_runx
lab.tf1_x = "runx1_10nM"
lab.tf1_y = "runx1_50nM"
lab.tf2_y = "runx2_10nM"
###finish hand typing input###

#The result of varRep function on two replicas of the same TF
res_varRep <- varRep(data_fit[,lab.tf1_x], data_fit[,lab.tf1_y], plot = T,
                     xlab = lab.tf1_x, ylab = lab.tf1_y)
#The result of classTF function on two different TFs
res_pref <- classTF(data_fit[,lab.tf1_x], data_fit[,lab.tf2_y], plot = T, 
                       xlab = lab.tf1_x, ylab = lab.tf2_y,
                       t = res_varRep$t, sigma = res_varRep$sigma)
#Write output file
data_output <- data_fit[,c(lab.tf1_x, lab.tf2_y)]
data_output$seq <- data_fit$seq
data_output$pref <- res_pref

#output the updated data, with the last column indicating preference
write.table(data_output, file = paste0("Preference_calls_",lab.tf1_x,"_vs_",lab.tf2_y,"_0.99CI.txt"), sep="\t",quote=F, row.names = F)

#Draw colored scatterplot
library(RColorBrewer)
library(classInt)
pref_tf1 = data_output[which(data_output$pref > 0), ] ## TF1 preferred sites
pref_tf2 = data_output[which(data_output$pref < 0), ] ## TF2 preferred sites
plotvar1 = pref_tf1$pref
plotclr1 = brewer.pal(9,"Reds")
class1 = classIntervals(plotvar1,9,syle="quantile")
colcode1=findColours(class1,plotclr1)
plotvar2 = pref_tf2$pref
plotclr2 = rev(brewer.pal(9,"Blues"))
class2 = classIntervals(plotvar2,9,syle="quantile")
colcode2 = findColours(class2,plotclr2)
plot(pref_tf1[,lab.tf1_x],pref_tf1[,lab.tf2_y],cex=.3,pch=16,col=colcode1,xlim=c(0,1),ylim=c(0,1))
points(data_output[which(data_output$pref==0),lab.tf1_x],data_output[which(data_output$pref==0),lab.tf2_y],cex=.3,pch=16,col="plum1")
points(pref_tf2[,lab.tf1_x],pref_tf2[,lab.tf2_y],cex=.3,pch=16,col=colcode2)

################################################################
# Part 2. Identify differentially preferred sites for genomic prediction matrix
################################################################
## Example using ETS
# Read in prediction files for all Elk1 ChIP-seq peaks using Elk1 (and Ets1) SVR model
Pred_elk1_allpeak_elk1 = read.table("GenomicPredictions/Elk1_k562_idr0.05_allpeak_201bp_elk1_100nM_bound_filtered_normalized_transformed_20bp_sepcores_allvalues.txt")
Pred_elk1_allpeak_ets1 = read.table("GenomicPredictions/Elk1_k562_idr0.05_allpeak_201bp_ets1_100nM_bound_filtered_normalized_transformed_20bp_sepcores_allvalues.txt")
Pred_matrix_tf1 = Pred_elk1_allpeak_elk1
Pred_matrix_tf2 = Pred_elk1_allpeak_ets1
## Estimate variance structure from replicate data, Fit regression line from TF pair PBM data, and calculate pref scores for each predicted values after model fitting
## Get preference matrix for genomic coordinates
data_fit = data_ets
lab.tf1_x <- "elk1_100nM"  ## Used in Replicate estimate (VarRep) and TF pair estimate (ClassTF)
lab.tf1_y <- "elk1_50nM"   ## Used in Replicate estimate
lab.tf2_y <- "ets1_100nM"  ## Used in TF pair estimate
res_varRep <- varRep(data_fit[,lab.tf1_x], data_fit[,lab.tf1_y], plot = F)
pred_matrix_pref = class_predTF(data_fit[,lab.tf1_x],data_fit[,lab.tf2_y],Pred_matrix_tf1,Pred_matrix_tf2,t=res_varRep$t,sigma=res_varRep$sigma,level=0.99,plot=F)

#### An additional function to check if preferred sites are also above NegCtrl range for the TF preferred
verify_neg = function(matrix_pref,matrix_occ1,matrix_occ2,neg1=0.2244,neg2=0.2128){ #neg1:Elk1 neg2:Ets1
  matrix_occ1 = matrix_occ1[,-(1:3)]
  matrix_occ2 = matrix_occ2[,-(1:3)]
  for(i in 1:length(matrix_pref[,1])){
    for(j in 1:length(matrix_pref[1,])){
      if(matrix_pref[i,j] == -1 && matrix_occ1[i,j] <= neg1){
        matrix_pref[i,j] = 0
      }else if(matrix_pref[i,j] == 1 && matrix_occ2[i,j] <= neg2){
        matrix_pref[i,j] = 0
      }
    }
  }
  return(matrix_pref)
}

Pred_elk1_allpeak_pref_negup = verify_neg(Pred_elk1_allpeak_pref,Pred_elk1_allpeak_elk1,Pred_elk1_allpeak_ets1) ## Test Elk1 ChIP-seq peak regions
Pred_ets1_allpeak_pref_negup = verify_neg(Pred_ets1_allpeak_pref,Pred_ets1_allpeak_elk1,Pred_ets1_allpeak_ets1) ## Test Ets1 ChIP-seq peak regions
