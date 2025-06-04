
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
outpath <- paste(project_path, 'results/revision/sex/structure/control_for_meanORtiv/', sep = "")

# # file names
# combined_df_ct_noqc
# combined_df_ct_artifact
# combined_df_ct_artifact_harmonized

# Define modeling parameters
dataset <- 'combined' # combined data across studies
qc_version <- 'artifact' # can be 'artifact' or 'noqc'
gamtype <- 'sex'
harmonized <- TRUE # whether to use harmonized data
controlformean <- FALSE # whether to include mean map as covariate
controlforTIV <- TRUE # whether to include mean map as covariate

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

if (controlforTIV == TRUE){
  dtype <- sprintf('%s_eTIV', dtype)
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
  filter(sex %in% c("Male", "Female"))

metric.schaefer400.all$sex <- factor(metric.schaefer400.all$sex, 
                                     levels = c("Male", "Female"), 
                                     ordered = FALSE)

# # ordered
# metric.schaefer400.all$oSex <- factor(metric.schaefer400.all$sex,
#                                      levels = c("Male", "Female"),
#                                      ordered = TRUE)

# If enabled, adjust covariates to control for mean cortical value
if (controlformean == TRUE){
  if (harmonized == TRUE){
    meanmapcovar <- metric.schaefer400.all$meanValHarmonized
    covars <- 'euler + meanmapcovar'}
  else if (harmonized == FALSE){
    meanmapcovar <- metric.schaefer400.all$meanVal
    covars <- 'euler + meanmapcovar'}
} else if (controlforTIV == TRUE){
  tivmapcovar <- metric.schaefer400.all$eTIV
  covars <- 'euler + tivmapcovar'
} else if (controlformean == FALSE && controlforTIV == FALSE){
  covars <- 'euler'
}

smooth_var <- 'age'
# in this analysis, ordered sex or categorical sex (unordered) will give 
# identical results (because i only have 2 levels)
model_var <- 'sex'

##############################
# Section 1: Parcel-wise GAM Analyses for Sex Effects
##############################
# Loop through 400 Schaefer parcels and run region-specific GAMs
# Collect GAM stats (F, p, partial R2, etc.)
# Save output as .csv

# Merge with SA rank and adjust p-values (FDR)
# Plot histograms (partial R2 dist) and brain maps:
#   - partial R2 values
#   - FDR-significant p-values
#   - R2-based rank and rank significance maps

if(dataset == 'combined'){
  #### Region-wise GAM Statistics and Derivative-based Temporal Developmental Properties ####
  #list of regions to run gam.fit.smooth function on below
  schaefer.regions <- names(metric.schaefer400.all[0:400]) %>% as.data.frame() %>% set_names("region")
  
  # fit GAMs for sex
  gam.variable.schaefer <- matrix(data=NA, nrow=400, ncol=5)
  #for each schaefer region
  for(row in c(1:nrow(schaefer.regions))){
    region <- schaefer.regions$region[row]
    GAM.RESULTS <- gam.fit.linear(measure = "metric", atlas = "schaefer400",
                                  dataset = "all",
                                  region = region,
                                  smooth_var = smooth_var,
                                  linear_var = model_var,
                                  covariates = covars,
                                  knots = 3, set_fx = FALSE)
    #and append results to output df
    gam.variable.schaefer[row,] <- GAM.RESULTS}
  
  gam.variable.schaefer <- as.data.frame(gam.variable.schaefer)
  colnames(gam.variable.schaefer) <- c("region","GAM.variable.Fvalue","GAM.variable.pvalue",
                                       "GAM.variable.partialR2","Anova.variable.pvalue")
  cols = c(2:5)
  gam.variable.schaefer[,cols] = apply(gam.variable.schaefer[,cols], 2,
                                       function(x) as.numeric(as.character(x)))

  
  write.csv(gam.variable.schaefer, paste(outpath, sprintf('csvFiles/%s_%s_statistics.csv',
                                                          dtype, gamtype),
                                         sep=""),
            row.names = F, quote = F)
  rm(gam.variable.schaefer)
  gc()
  
  # re-read the results from above and compare with SA
  #GAM age smooth statistics, generated with fitGAMs_fluctuationamplitude_age.R
  gam.variable.schaefer <- read.csv(paste(outpath, sprintf('csvFiles/%s_%s_statistics.csv',
                                                           dtype, gamtype),
                                          sep=""))
  
  # SA axis
  sa.schaefer400 <- read.csv(paste(project_path,
                                   'data/SArank_schaefer400_7Networks.csv',
                                   sep = ""))
  
  # sa.schaefer400 <- sa.schaefer400 %>% select(-X)
  colnames(sa.schaefer400) <- c("SA.rank", "region")
  
  gam.variable.schaefer$region <- gsub("X", "", gam.variable.schaefer$region)
  
  gam.variable.schaefer <- merge(gam.variable.schaefer, sa.schaefer400, by = "region")
  
  csvR2 <- data.frame(gam.variable.schaefer$region)
  csvR2$partialR2 <- gam.variable.schaefer$GAM.variable.partialR2
  # pvlues
  # GAMs
  pvalues = gam.variable.schaefer$GAM.variable.pvalue
  GAMpvaluesfdrs<-p.adjust(pvalues, method="BH")
  
  # Anova
  pvalues = gam.variable.schaefer$Anova.variable.pvalue
  Anovapvaluesfdrs<-p.adjust(pvalues, method="BH")
  
  csvR2$anovaPvaluefdr <- Anovapvaluesfdrs
  csvR2$gamPvaluefdr <- GAMpvaluesfdrs
  outputPath <- paste(outpath, sprintf('csvFiles/%s_%s_r2.csv', dtype, gamtype),
                      sep="")
  write.csv(csvR2, outputPath, row.names=FALSE)
  
  # Effect size
  # histogram
  ggplot(gam.variable.schaefer, aes(x = GAM.variable.partialR2)) +
    geom_histogram(binwidth=.01, fill="darkcyan", color="#e9ecef", alpha=0.9) +
    theme_bw()
  
  ggsave(filename = paste(outpath, sprintf('%s_%s_histogram_partialR2.svg', dtype,
                                           gamtype),
                          sep = ""),
         dpi = 300, width = 3 , height = 2)
  
  # brain
  maxval <- max(abs(gam.variable.schaefer$GAM.variable.partialR2))
  
  ggseg(.data = gam.variable.schaefer, atlas = "schaefer7_400",
        mapping=aes(fill = GAM.variable.partialR2, colour=I("#e9ecef"),
                    size=I(.03)), position = c("stacked")) +
    theme_void() +
    paletteer::scale_fill_paletteer_c("pals::warmcool",
                                      na.value="transparent", direction = -1,
                                      limits = c(-maxval, maxval),
                                      oob = squish)
  
  ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_partialR2.svg', dtype,
                                           gamtype),
                          sep = ""),
         dpi = 300, width = 3 , height = 2)
  
  
  # Replace GAM.variable.partialR2 with ranks robustly
  # (1) Rank the partial R2 values
  gam.variable.schaefer$GAM.variable.rankR2 <- rank(
    gam.variable.schaefer$GAM.variable.partialR2, 
    ties.method = "average"
  )
  
  # (2) Find the nearest negative value, if it exists
  neg_vals <- gam.variable.schaefer$GAM.variable.partialR2[
    gam.variable.schaefer$GAM.variable.partialR2 < 0
  ]
  
  if (length(neg_vals) > 0) {
    # (3) Center ranks around the largest negative value (i.e., nearest to 0)
    nearestNeg <- max(neg_vals)
    nearestNegIdx <- which(gam.variable.schaefer$GAM.variable.partialR2 == nearestNeg)
    nearestNegRank <- gam.variable.schaefer$GAM.variable.rankR2[nearestNegIdx]
    
    gam.variable.schaefer$GAM.variable.rankR2 <- 
      gam.variable.schaefer$GAM.variable.rankR2 - (nearestNegRank + 1)
    
    centered <- TRUE
  } else {
    # No negatives — do not shift the ranks
    centered <- FALSE
  }
  
  # (4) Max absolute value (for plotting scale)
  maxval <- max(abs(gam.variable.schaefer$GAM.variable.rankR2), na.rm = TRUE)
  
  
  # brain ranks
  ggseg(.data = gam.variable.schaefer, atlas = "schaefer7_400",
        mapping=aes(fill = GAM.variable.rankR2, colour=I("#e9ecef"), size=I(.03)),
        position = c("stacked")) + theme_void() +
    paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent",
                                      direction = -1,
                                      limits = c(-maxval, maxval),
                                      oob = squish)
  
  
  ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_rank_partialR2.svg',
                                           dtype, gamtype),
                          sep = ""),
         dpi = 300, width = 3 , height = 2)
  
  # brain pvalues
  # Anova
  pvalues = gam.variable.schaefer$Anova.variable.pvalue
  pvaluesfdrs<-p.adjust(pvalues, method="BH")
  
  Anovasignumber = sum(pvaluesfdrs < 0.05, na.rm=TRUE)
  pvaluesfdrs[pvaluesfdrs >= 0.05] <- NA
  gam.variable.schaefer$Anova.variable.pvaluefdr <- pvaluesfdrs
  
  ggseg(.data = gam.variable.schaefer, atlas = "schaefer7_400",
        mapping=aes(fill = Anova.variable.pvaluefdr, colour=I("#e9ecef"), size=I(.03)),
        position = c("stacked")) + theme_void() + ggtitle(Anovasignumber) +
    paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent",
                                      direction = -1,
                                      limits = c(0, 0.05),
                                      # limits = c(min(metric.regional.statistics$GAM.smooth.pvalue),
                                      #            max(metric.regional.statistics$GAM.smooth.pvalue)),
                                      oob = squish)
  
  ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_pval_partialR2.svg',
                                           dtype, gamtype),
                          sep = ""),
         dpi = 300, width = 3 , height = 2)
  
  # significant ranks
  pvalues = gam.variable.schaefer$Anova.variable.pvalue
  pvaluesfdrs <- p.adjust(pvalues, method="BH")
  rankR2sig <- gam.variable.schaefer$GAM.variable.rankR2
  rankR2sig[(pvaluesfdrs >= 0.05)] <- NA
  gam.variable.schaefer$GAM.variable.rankR2sig <- rankR2sig
  
  maxval <- max(abs(gam.variable.schaefer$GAM.variable.rankR2sig), na.rm=T)
  
  # brain significant ranks
  ggseg(.data = gam.variable.schaefer, atlas = "schaefer7_400",
        mapping=aes(fill = GAM.variable.rankR2sig, colour=I("#e9ecef"), size=I(.03)),
        position = c("stacked")) + theme_void() + ggtitle(Anovasignumber) +
    paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent",
                                      direction = -1,
                                      limits = c(-maxval, maxval),
                                      oob = squish)
  
  ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_ranksig_partialR2.svg',
                                           dtype, gamtype),
                          sep = ""),
         dpi = 300, width = 3 , height = 2)
}
