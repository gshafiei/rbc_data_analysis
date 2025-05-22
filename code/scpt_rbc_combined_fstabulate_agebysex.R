
rm(list = ls())

# Load required libraries and helper GAM functions
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
outpath <- paste(project_path, 'results/revision/agebysex/structure/', sep = "")

# # file names
# combined_df_ct_noqc
# combined_df_ct_artifact
# combined_df_ct_artifact_harmonized

# Define modeling parameters
dataset <- 'combined' # combined data across studies
qc_version <- 'artifact' # can be 'artifact' or 'noqc'
gamtype <- 'agebysex'
harmonized <- TRUE # whether to use harmonized data

# Metric of interest (e.g., 'ct', 'sa', 'gv', 'lgi') and associated map
metric <- 'sa'
corticalmap <- 'meanVal'

# Dynamically build filename for the dataset based on flags
# Selects appropriate TSV input depending on predictors and harmonization status
# Output: dtype (dataset label string)
if (harmonized == TRUE){corticalmap <- 'meanValHarmonized'}

dtype <- sprintf('%s_df_%s_%s', dataset, metric, qc_version)
if (harmonized == TRUE){
  dtype <- sprintf('%s_df_%s_%s_harmonized', dataset, metric, qc_version)
}

# Load and prepare the dataset for GAMs
# Schaefer-400 parcellation -- structural data
metric.schaefer400.all <- read.csv(paste(data_path, sprintf('%s.tsv', dtype),
                                         sep = ""),
                                   sep = '\t')
# will use dataframe as a covariate
metric.schaefer400.all$dataset <- as.factor(metric.schaefer400.all$study)
# will use sex as a covariate
metric.schaefer400.all$sex <- as.factor(metric.schaefer400.all$sex)
# add sex as ordered factor
# metric.schaefer400.all$oSex <- ordered(metric.schaefer400.all$sex, 
#                                        levels = c('Female','Male', 'Other'))

# Filter out rows where sex is not "female" or "male"
metric.schaefer400.all <- metric.schaefer400.all %>%
  filter(sex %in% c("Female", "Male"))

metric.schaefer400.all$oSex <- factor(metric.schaefer400.all$sex, 
                                     levels = c("Female", "Male"), 
                                     ordered = TRUE)

smooth_var <- 'age'
covars <- 'oSex + euler'
model_var <- 'oSex'

###############################
# Section 1: Overall Fit for Mean Value
# Fit a gam for combined mean values
###############################
# Fit a GAM using mean value (e.g., meanVal or meanValHarmonized)
# gam.fit.smooth() for age-related analysis: age as smooth term
# gam.fit.linear() for linear p-factor modeling
# Predict fitted values with confidence intervals using custom predict functions
# Plot curve with ggplot2 and save as SVG

#### PREDICT GAM SMOOTH FITTED VALUES ####

GAM.RESULTS <- gam.factorsmooth.interaction(measure = 'metric', 
                                            atlas = 'schaefer400',
                                            dataset = 'all', 
                                            region = corticalmap,
                                            smooth_var = smooth_var, 
                                            int_var = model_var,
                                            covariates = covars,
                                            knots = 3, set_fx = FALSE)
gam.results <- as.data.frame(GAM.RESULTS)

# predict smooth fitted values for 'Female'
ctx.predicted.metric.F <- gam.smooth.predict.covariateinteraction(measure = 'metric',
                                                                  atlas = 'schaefer400',
                                                                  dataset = 'all',
                                                                  region = corticalmap,
                                                                  smooth_var = smooth_var,
                                                                  int_var = model_var,
                                                                  int_var.predict = 'Female',
                                                                  covariates = covars,
                                                                  knots = 3, 
                                                                  set_fx = FALSE,
                                                                  increments = 200)
ctx.predicted.metric.F$sex <- 'Female'

# predict smooth fitted values for 'Male'
ctx.predicted.metric.M <- gam.smooth.predict.covariateinteraction(measure = 'metric',
                                                                  atlas = 'schaefer400',
                                                                  dataset = 'all',
                                                                  region = corticalmap,
                                                                  smooth_var = smooth_var,
                                                                  int_var = model_var,
                                                                  int_var.predict = 'Male',
                                                                  covariates = covars,
                                                                  knots = 3, 
                                                                  set_fx = FALSE,
                                                                  increments = 200)
ctx.predicted.metric.M$sex <- 'Male'

# combine
ctx.predicted.metric <- bind_rows(ctx.predicted.metric.F, ctx.predicted.metric.M)


# plot age
if(metric == 'ct'){
  lolim <- 2.0
  hilim <- 3.5
} # ct
if(metric == 'sa'){
  lolim <- 110
  hilim <- 700
} # sa
if(metric == 'gv'){
  lolim <- 350
  hilim <- 2200
} # gv
if(metric == 'lgi'){
  lolim <- 2.0
  hilim <- 4.0
} # lgi

ggplot(metric.schaefer400.all, aes(x = age,
                                   if (harmonized == TRUE){y = meanValHarmonized}
                                   else{y = meanVal})) +
  geom_point(aes(color = sex), alpha = 0.25, size = 1) +
    geom_ribbon(data = ctx.predicted.metric, 
              aes(x = age, y = .fitted, ymin = .lower_ci, ymax = .upper_ci, fill = sex), 
              alpha = 0.5, inherit.aes = FALSE) +
    geom_line(data = ctx.predicted.metric, 
            aes(x = age, y = .fitted, color = sex), 
            size = 1.2, inherit.aes = FALSE) +
  labs(x='\nage', y=sprintf('%s\n', corticalmap),
       title=sprintf('F-stat=%s\npval=%s\n',
                     gam.results$gam.int.F, gam.results$gam.int.pvalue),
       color = "Sex", fill = "Sex") +
  theme_classic() +
  theme(
    axis.text = element_text(size=12, family = "Arial", color = c("black")),
    axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
    axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
  theme(legend.position="right") +
  theme(aspect.ratio=1) +
  scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), limits = c(6,22), expand = c(0.05,.05)) +
  ylim(lolim, hilim)

ggsave(paste(outpath, sprintf('%s_%s.svg', dtype, gamtype), sep = ""),
       dpi = 300,
       plot = last_plot())

##############################
# Section 2: Parcel-wise GAM Analyses
##############################
