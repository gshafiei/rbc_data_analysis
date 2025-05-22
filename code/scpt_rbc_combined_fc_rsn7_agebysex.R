
rm(list = ls())
# Load libraries and helper functions
library(stringr)
library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)
library(ggseg)
library(ggsegSchaefer)
library(paletteer)
library(pals)
library(ggseg3d)
library(ggplot2)
library(scales)

source("/Users/gshafiei/Desktop/rbc_data_analysis/code/func_GAM_rbc.R")

# Set paths
project_path <- '/Users/gshafiei/Desktop/RBC/'
data_path <- paste(project_path, 'data/dataR/', sep = "")
outpath <- paste(project_path, 'results/revision/agebysex/function/', sep = "")


# # file names
# combined_df_withinbetween_fcrsn7_noqc
# combined_df_withinbetween_fcrsn7_artifact
# combined_df_withinbetween_fcrsn7_artifact_harmonized

# Define modeling parameters
dataset <- 'combined' # combined data across studies
qc_version <- 'artifact' # can be 'artifact' or 'noqc'
gamtype <- 'agebysex'
harmonized <- TRUE # whether to use harmonized data

# Dynamically build filename for the dataset based on flags
# Selects appropriate TSV input depending on predictors and harmonization status
# Output: dtype (dataset label string)
dtype <- sprintf('%s_df_withinbetween_fcrsn7_%s', dataset,
                 qc_version)
if (harmonized == TRUE){
  dtype <- sprintf('%s_df_withinbetween_fcrsn7_%s_harmonized', dataset, qc_version)
}

# Load and prepare the dataset for GAMs
# Schaefer-400 parcellation -- functional connectivity data
netpair.schaefer400.all <- read.csv(paste(data_path,
                                          sprintf('%s.tsv', dtype),
                                          sep = ""),
                                    sep = '\t')

# will use dataframe as a covariate
netpair.schaefer400.all$dataset <- as.factor(netpair.schaefer400.all$study)
# will use sex as a covariate
netpair.schaefer400.all$sex <- as.factor(netpair.schaefer400.all$sex)

# Filter out rows where sex is not "female" or "male"
netpair.schaefer400.all <- netpair.schaefer400.all %>%
  filter(sex %in% c("Female", "Male"))

# ordered
netpair.schaefer400.all$oSex <- factor(netpair.schaefer400.all$sex,
                                       levels = c("Female", "Male"),
                                       ordered = TRUE)

smooth_var <- 'age'
covars <- 'oSex + medianFD'
model_var <- 'oSex'

##############################
# Section 1: GAM Model Fitting
# Fit GAMs for Network Pairs
##############################
# For each network pair (28 network pairs), fit GAMs using helper functions
# gam.fit.smooth() for age-related analysis: age as smooth term
# gam.fit.linear() for linear p-factor modeling
# Save p-values, partial R2, and onset/peak values if age
# Adjust p-values (FDR correction) and export result tables as CSV

if(dataset == 'combined'){
  #### Run gam.fit.smooth in all networks: with k=3
  # # list of regions to run gam.fit.smooth function on below
  # netpair_labels <- colnames(netpair.schaefer400.all)[0:28]
  netpair_labels <- names(netpair.schaefer400.all[0:28]) %>% 
    as.data.frame() %>% set_names("netpair")
  
  # network connectivity GAMs
  #empty matrix to save gam.fit output to
  gam.variable <- matrix(data=NA, nrow=28, ncol=3)
  #for each network pair
  for(row in c(1:nrow(netpair_labels))){
    netpair <- netpair_labels$netpair[row]
    #run the gam.fit.smooth function
    GAM.RESULTS <- gam.factorsmooth.interaction(measure = 'netpair', 
                                                atlas = 'schaefer400',
                                                dataset = 'all', 
                                                region = netpair,
                                                smooth_var = smooth_var, 
                                                int_var = model_var,
                                                covariates = covars,
                                                knots = 3, set_fx = FALSE)
    gam.variable[row,] <- GAM.RESULTS}
  
  gam.variable <- as.data.frame(gam.variable)
  colnames(gam.variable) <- c("netpair","GAM.variable.Fvalue","GAM.variable.pvalue")
  
  # pvalues
  pvalues = gam.variable$GAM.variable.pvalue
  pvaluesfdrs<-p.adjust(pvalues, method="BH")
  gam.variable$GAM.variable.pvaluefdr <- pvaluesfdrs
  
  # make sure values are numerical
  cols = c(2:4)
  gam.variable[,cols] = apply(gam.variable[,cols], 2,
                              function(x) as.numeric(as.character(x)))
  write.csv(gam.variable, paste(outpath,
                                sprintf('csvFiles/%s_%s_statistics.csv',
                                        dtype, gamtype),
                                sep = ""),
            row.names = F, quote = F)
  rm(gam.variable)
  gc()
}
