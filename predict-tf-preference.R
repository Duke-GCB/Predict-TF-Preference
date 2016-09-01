#!/usr/bin/env Rscript
# Load the WLS functions
source("Rcode/wls_functions.R")

PREDICTION_COLUMN = 4 # 1-indexed column in the bed file where the fractional (0.0 - 1.0) scores are

dan_predTF <- function(scores, labels, data_fit) {
  res_varRep <- varRep(data_fit[,labels.tf1_x], data_fit[,labels.tf1_y], plot = F)
  return(class_predTF(data_fit[,labels.tf1_x],data_fit[,labels.tf2_y],scores.x, scores.y,t=res_varRep$t,sigma=res_varRep$sigma,level=0.99,plot=T))
}

extract_scores <- function(bed_file) {
  # Read the bedfile as a table
  bed_data <- read.table(bed_file, header=FALSE, stringsAsFactors = FALSE)
  # Return only the scores from the prediction column
  return(bed_data[, PREDICTION_COLUMN])
}

load_fit_data <- function(family) {
  # hard-coding the column names for the PBM data here
  data = switch (family,
    ets = read.table("PBMdata/Combined_Ets1_100nM_Elk1_100nM_50nM_Gabpa_100nM_Bound_filtered_normalized.txt",col.names=c("seq","ets1_100nM","elk1_100nM","elk1_50nM","gabpa_100nM")),
    e2f = read.table("PBMdata/Combined_E2f1_250nM_200nM_E2f3_250nM_E2f4_500nM_800nM_Bound_filtered_normalized.txt",col.names = c("seq","e2f1_250nM","e2f1_200nM","e2f3_250nM","e2f4_500nM","e2f4_800nM")),
    bHLH = read.table("PBMdata/Combined_Max_Myc_Mad_MadL_Bound_filtered_normalized.txt",col.names = c("seq","max","myc","mad","madL")),
    runx = read.table("PBMdata/Combined_Runx1_10nM_50nM_Runx2_10nM_Bound_filtered_normalized.txt",col.names = c("seq","runx1_10nM","runx1_50nM","runx2_10nM"))
  )
  return(data)
}

args = commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop("Usage: predict-tf-preference.R tf1.bed tf2.bed family tf1_x tf1_y tf2_y")
}

# read the scores - arguments 1 and 2 have the bed files
tf_predictions_1_bed_file <- args[1]
tf_predictions_2_bed_file <- args[2]
# argument 3 has the family, which identifies the PBMData and its columns
family <- args[3]

# arguments 4, 5, 6 have the labels tf1_x, tf1_y, and tf2_y
labels.tf1_x <- args[4]  ## Used in Replicate estimate (VarRep) and TF pair estimate (ClassTF)
labels.tf1_y <- args[5]  ## Used in Replicate estimate
labels.tf2_y <- args[6]  ## Used in TF pair estimate

scores.x = extract_scores(tf_predictions_1_bed_file)
scores.y = extract_scores(tf_predictions_2_bed_file)
data_fit = load_fit_data(family)

result = dan_predTF(scores, labels, data_fit)
print(result)
