
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
  filter(sex %in% c("Male", "Female"))

# ordered
netpair.schaefer400.all$oSex <- factor(netpair.schaefer400.all$sex,
                                       levels = c("Male", "Female"),
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
  # also estimate a sign for F-stat by predicting outputs for F and M and finding the average difference
  gam.variable.sign <- matrix(data=NA, nrow=28, ncol=1)
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
    gam.variable[row,] <- GAM.RESULTS
    
    # predict smooth fitted values for 'Female'
    predicted.metric.F <- gam.smooth.predict.covariateinteraction(measure = 'netpair',
                                                                  atlas = 'schaefer400',
                                                                  dataset = 'all',
                                                                  region = netpair,
                                                                  smooth_var = smooth_var,
                                                                  int_var = model_var,
                                                                  int_var.predict = 'Female',
                                                                  covariates = covars,
                                                                  knots = 3, 
                                                                  set_fx = FALSE,
                                                                  increments = 200)
    predicted.metric.F$sex <- 'Female'
    
    # predict smooth fitted values for 'Male'
    predicted.metric.M <- gam.smooth.predict.covariateinteraction(measure = 'netpair',
                                                                  atlas = 'schaefer400',
                                                                  dataset = 'all',
                                                                  region = netpair,
                                                                  smooth_var = smooth_var,
                                                                  int_var = model_var,
                                                                  int_var.predict = 'Male',
                                                                  covariates = covars,
                                                                  knots = 3, 
                                                                  set_fx = FALSE,
                                                                  increments = 200)
    predicted.metric.M$sex <- 'Male'
    
    # differences
    diff <-predicted.metric.F$.fitted - predicted.metric.M$.fitted
    mean_diff <- mean(diff, na.rm = TRUE)
    
    # append sign to output df
    gam.variable.sign[row,] <- mean_diff
    
    # combine
    predicted.metric <- bind_rows(predicted.metric.F, predicted.metric.M)
    
    }
  
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
  
  # add sign
  gam.variable.sign <- as.data.frame(gam.variable.sign)
  colnames(gam.variable.sign) <- c("meandiff")
  
  signed_F <- gam.variable$GAM.variable.Fvalue * sign(gam.variable.sign$meandiff)
  gam.variable$GAM.variable.Fvalue.signed <- signed_F
  
  write.csv(gam.variable, paste(outpath,
                                sprintf('csvFiles/%s_%s_statistics_signed.csv',
                                        dtype, gamtype),
                                sep = ""),
            row.names = F, quote = F)
  rm(gam.variable)
  gc()
}

###############################
# Section 2: Single Example Fits + Plot
# Fit a gam for certain network pairs as examples
###############################
# Fit and visualize GAM curve for one network pair (e.g., Cont.Default)
# Generate ribbon plots and save to SVG using ggplot2
withinrsn <- FALSE
#### PREDICT GAM SMOOTH FITTED VALUES ####
if(withinrsn == TRUE){
  networkpair <- 'Cont.Cont'
  } else{networkpair <- 'Cont.Default'} # SalVentAttn.Vis or Default.Cont

# predict smooth fitted values for 'Female'
ctx.predicted.metric.F <- gam.smooth.predict.covariateinteraction(measure = 'netpair',
                                                                  atlas = 'schaefer400',
                                                                  dataset = 'all',
                                                                  region = networkpair,
                                                                  smooth_var = smooth_var,
                                                                  int_var = model_var,
                                                                  int_var.predict = 'Female',
                                                                  covariates = covars,
                                                                  knots = 3, 
                                                                  set_fx = FALSE,
                                                                  increments = 200)
ctx.predicted.metric.F$sex <- 'Female'

# predict smooth fitted values for 'Male'
ctx.predicted.metric.M <- gam.smooth.predict.covariateinteraction(measure = 'netpair',
                                                                  atlas = 'schaefer400',
                                                                  dataset = 'all',
                                                                  region = networkpair,
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


# plot
ggplot(data = netpair.schaefer400.all, aes(x = age,
                                           if (withinrsn == TRUE){y = Cont.Cont}
                                           else{y = .data[[networkpair]]})) +
  geom_point(aes(color = sex), alpha = 0.25, size = 1) +
  geom_ribbon(data = ctx.predicted.metric, 
              aes(x = age, y = .fitted, ymin = .lower_ci, ymax = .upper_ci, fill = sex), 
              alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = ctx.predicted.metric, 
            aes(x = age, y = .fitted, color = sex), 
            size = 1.2, inherit.aes = FALSE) +
  labs(x='\nage', y=sprintf('%s\n', networkpair)) +
  theme_classic() +
  theme(
    axis.text = element_text(size=12, family = "Arial", color = c("black")),
    axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
    axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
  theme(legend.position="none") +
  theme(aspect.ratio=1) +
  scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), limits = c(6,22), expand = c(0.05,.05)) # +
  # (if (withinrsn == TRUE){ylim(0.03, 0.8)}
  #  else{ylim(-0.40, 0.75)})

ggsave(paste(outpath, sprintf('%s_%s_%s_v2.svg', dtype, gamtype, networkpair), sep = ""),
       dpi = 300,
       plot = last_plot())
