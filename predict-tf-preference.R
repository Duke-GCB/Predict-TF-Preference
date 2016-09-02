#!/usr/bin/env Rscript
# Load the WLS functions

# Hack to source wls_functions.R relatively
# This is necessary because cwl changes the Docker workdir
# via http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- paste(sep="/", script.basename, "Rcode/wls_functions.R")
source(other.name)

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

write_preferences <- function(prefs, template_bedfile, output_file) {
  # Writes out a BED file with preferences. Starts by reading the template_bedfile
  # and replaces the values in PREDICTION_COLUMN with prefs
  bed_data <- read.table(template_bedfile, header=FALSE, stringsAsFactors = FALSE)
  bed_data[, PREDICTION_COLUMN] = prefs
  write.table(bed_data, output_file, quote = F, col.names = F, row.names = F)
}

load_fit_data <- function(family) {
  # hard-coding the column names for the PBM data here
  data = switch (family,
    ets = read.table(paste(sep="/", script.basename,"PBMdata/Combined_Ets1_100nM_Elk1_100nM_50nM_Gabpa_100nM_Bound_filtered_normalized.txt"),col.names=c("seq","ets1_100nM","elk1_100nM","elk1_50nM","gabpa_100nM")),
    e2f = read.table(paste(sep="/", script.basename,"PBMdata/Combined_E2f1_250nM_200nM_E2f3_250nM_E2f4_500nM_800nM_Bound_filtered_normalized.txt"),col.names = c("seq","e2f1_250nM","e2f1_200nM","e2f3_250nM","e2f4_500nM","e2f4_800nM")),
    bHLH = read.table(paste(sep="/", script.basename,"PBMdata/Combined_Max_Myc_Mad_MadL_Bound_filtered_normalized.txt"),col.names = c("seq","max","myc","mad","madL")),
    runx = read.table(paste(sep="/", script.basename,"PBMdata/Combined_Runx1_10nM_50nM_Runx2_10nM_Bound_filtered_normalized.txt"),col.names = c("seq","runx1_10nM","runx1_50nM","runx2_10nM"))
  )
  return(data)
}

args = commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  stop("Usage: predict-tf-preference.R tf1.bed tf2.bed output.bed family tf1_x tf1_y tf2_y")
}

# read the scores - arguments 1 and 2 have the bed files
tf_predictions_1_bed_file <- args[1]
tf_predictions_2_bed_file <- args[2]
# Output file name in argument 3
output_preferences_bed_file <- args[3] # output bed file

# argument 4 has the family, which identifies the PBMData and its columns
family <- args[4]

# arguments for labels tf1_x, tf1_y, and tf2_y
labels.tf1_x <- args[5]  ## Used in Replicate estimate (VarRep) and TF pair estimate (ClassTF)
labels.tf1_y <- args[6]  ## Used in Replicate estimate
labels.tf2_y <- args[7]  ## Used in TF pair estimate

scores.x = extract_scores(tf_predictions_1_bed_file)
scores.y = extract_scores(tf_predictions_2_bed_file)
data_fit = load_fit_data(family)

prefs = dan_predTF(scores, labels, data_fit)
write_preferences(prefs, tf_predictions_1_bed_file, output_preferences_bed_file)
