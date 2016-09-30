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
    ETS = read.table(paste(sep="/", script.basename,"PBMdata/Combined_Ets1_100nM_Elk1_100nM_50nM_Gabpa_100nM_Bound_filtered_normalized.txt"),col.names=c("seq","ets1_100nM","elk1_100nM","elk1_50nM","gabpa_100nM")),
    E2F = read.table(paste(sep="/", script.basename,"PBMdata/Combined_E2f1_250nM_200nM_E2f3_250nM_E2f4_500nM_800nM_Bound_filtered_normalized.txt"),col.names = c("seq","e2f1_250nM","e2f1_200nM","e2f3_250nM","e2f4_500nM","e2f4_800nM")),
    bHLH = read.table(paste(sep="/", script.basename,"PBMdata/Combined_Max_Myc_Mad_MadL_Bound_filtered_normalized.txt"),col.names = c("seq","max","myc","mad","madL")),
    RUNX = read.table(paste(sep="/", script.basename,"PBMdata/Combined_Runx1_10nM_50nM_Runx2_10nM_50nM_Bound_filtered_normalized_V2.txt"),col.names = c("seq","runx1_10nM","runx1_50nM","runx2_10nM"))
  )
  return(data)
}

family_from_tfs <- function(tf1, tf2) {
  family = switch(paste(tf1, tf2),
                  "Elk1 Ets1" = "ETS",
                  "E2f1 E2f4" = "E2F",
                  "E2f1 E2f3" = "E2F",
                  "Mad1 c-Myc" = "bHLH",
                  "Runx1 Runx2" = "RUNX")
  return(family)
}

labels_from_tfs <- function(tf1, tf2) {
  labels = switch(paste(tf1,tf2),
                  "Elk1 Ets1" = list(tf1_x="elk1_100nM", tf1_y="elk1_50nM", tf2_y="ets1_100nM"),
                  "E2f1 E2f4" = list(tf1_x="e2f1_250nM", tf1_y="e2f1_200nM", tf2_y="e2f4_500nM"),
                  "E2f1 E2f3" = list(tf0_x="e2f4_500nM", tf0_y="e2f4_800nM", tf1_x="e2f1_250nM", tf2_y="e2f3_250nM"),
                  "Mad1 c-Myc" = list(tf1_x="mad", tf1_y="madL", tf2_y="myc"),
                  "Runx1 Runx2" = list(tf1_x="runx1_10nM", tf1_y="runx1_50nM", tf2_y="runx2_10nM")
                  )
  return(labels)
}

calc_varRep <- function(data_fit, labels) {
  if(is.character(labels$tf0_x)) {
    return(varRep(data_fit[,labels$tf0_x], data_fit[,labels$tf0_y], plot = F,
                         xlab = labels$tf0_x, ylab = labels$tf0_y))
  } else {
    return(varRep(data_fit[,labels$tf1_x], data_fit[,labels$tf1_y], plot = F,
                  xlab = labels$tf1_x, ylab = labels$tf1_y))
  }
}

calc_prefs <- function(data_fit, res_varRep, labels, scores) {
  return(class_predTF(data_fit[,labels$tf1_x], data_fit[,labels$tf2_y], scores.x, scores.y,
                      t=res_varRep$t, sigma=res_varRep$sigma, level=0.99, plot=F))
}

args = commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: predict-tf-preference.R tf1 tf2 tf1.bed tf2.bed output.bed")
}

tf1 <- args[1]
tf2 <- args[2]
tf_predictions_1_bed_file <- args[3]
tf_predictions_2_bed_file <- args[4]
output_preferences_bed_file <- args[5] # output bed file
labels <- labels_from_tfs(tf1, tf2)
if(is.null(labels)) {
  stop("TFs ", tf1, " and ", tf2, " are not a valid combination")
}

family = family_from_tfs(tf1, tf2)
scores.x = extract_scores(tf_predictions_1_bed_file)
scores.y = extract_scores(tf_predictions_2_bed_file)
data_fit = load_fit_data(family)
res_varRep = calc_varRep(data_fit, labels)
prefs = calc_prefs(data_fit, res_varRep, labels, scores)
write_preferences(prefs, tf_predictions_1_bed_file, output_preferences_bed_file)
