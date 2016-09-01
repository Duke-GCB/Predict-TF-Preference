#!/bin/bash

# This is just a test to make sure arguments are parsed correctly
# The named bed data actually was not generated with the provided model
./predict-tf-preference.R \
  ~/Data/tf-dna-predictions/hardac-results/hg19-E2F1-chr1-CCGC-predictions.bed \
  ~/Data/tf-dna-predictions/hardac-results/hg19-E2F4-chr1-CCGC-predictions.bed \
  e2f \
  e2f1_250nM \
  e2f1_200nM \
  e2f4_500nM

