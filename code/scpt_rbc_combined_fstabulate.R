
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
outpath <- paste(project_path, 'results/revision/structure/', sep = "")

# # for pfactor
# combined_df_ct_noqc_pfactor_filter
# combined_df_ct_artifact_pfactor_filter
# combined_df_ct_artifact_pfactor_filter_harmonized
#
# # for age
# combined_df_ct_noqc
# combined_df_ct_artifact
# combined_df_ct_artifact_harmonized

# Define modeling parameters
dataset <- 'combined' # combined data across studies
qc_version <- 'artifact' # can be 'artifact' or 'noqc'
gamtype <- 'age' # can be 'age' or 'pfactor'
harmonized <- TRUE # whether to use harmonized data
controlformean <- FALSE # whether to include mean map as covariate

# Metric of interest (e.g., 'ct', 'sa', 'gv', 'lgi') and associated map
metric <- 'ct'
corticalmap <- 'meanVal'

# Dynamically build filename for the dataset based on flags
# Selects appropriate TSV input depending on predictors and harmonization status
# Output: dtype (dataset label string)
if (harmonized == TRUE){corticalmap <- 'meanValHarmonized'}

if (gamtype == 'age'){
  dtype <- sprintf('%s_df_%s_%s', dataset, metric, qc_version)
  if (harmonized == TRUE){
    dtype <- sprintf('%s_df_%s_%s_harmonized', dataset, metric, qc_version)
  }
}
if (gamtype == 'pfactor'){
  dtype <- sprintf('%s_df_%s_%s_pfactor_filter', dataset, metric, qc_version)
  if (harmonized == TRUE){
    dtype <- sprintf('%s_df_%s_%s_pfactor_filter_harmonized', dataset, metric, qc_version)
  }
}

if (gamtype == 'age'){smooth_var <- 'age'}
if (gamtype == 'age'){covars <- 'sex + euler'}

if (gamtype == 'pfactor'){smooth_var <- 'age'}
if (gamtype == 'pfactor'){covars <- 'sex + euler'}
if (gamtype == 'pfactor'){linear_var <- 'p_factor_mcelroy_harmonized_all_samples'}

# Load and prepare the dataset for GAMs
# Schaefer-400 parcellation -- structural data
metric.schaefer400.all <- read.csv(paste(data_path, sprintf('%s.tsv', dtype),
                                         sep = ""),
                                   sep = '\t')
# will use dataframe as a covariate
metric.schaefer400.all$dataset <- as.factor(metric.schaefer400.all$study)
# will use sex as a covariate
metric.schaefer400.all$sex <- as.factor(metric.schaefer400.all$sex)

# If enabled, adjust covariates to control for mean cortical value
if (controlformean == TRUE){
  if (harmonized == TRUE){
    meanmapcovar <- metric.schaefer400.all$meanValHarmonized
    covars <- 'sex + euler + meanmapcovar'}
  else if (harmonized == FALSE){
    meanmapcovar <- metric.schaefer400.all$meanVal
    covars <- 'sex + euler + meanmapcovar'}
}

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
if(gamtype == 'age'){
  GAM.RESULTS <- gam.fit.smooth(measure = 'metric', atlas = 'schaefer400',
                                dataset = 'all', region = corticalmap,
                                smooth_var = smooth_var, covariates = covars,
                                knots = 3, set_fx = FALSE, stats_only = FALSE)
  gam.results <- as.data.frame(GAM.RESULTS)

  ctx.predicted.metric <- gam.smooth.predict(measure = 'metric',
                                             atlas = 'schaefer400',
                                             dataset = 'all',
                                             region = corticalmap,
                                             smooth_var = smooth_var,
                                             covariates = covars,
                                             knots = 3, set_fx = FALSE,
                                             increments = 200)

  # get predicted.smooth df from function output
  ctx.predicted.metric <- as.data.frame(ctx.predicted.metric[3])

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
  geom_point(aes(color = dataset), size = 2) + # color = "#115c25"
  geom_ribbon(data = ctx.predicted.metric, aes(x = age, y = .fitted,
                                               ymin = .lower_ci,
                                               ymax = .upper_ci),
              alpha = .7, linetype = 0,) +
  geom_line(data = ctx.predicted.metric, aes(x = age, y = .fitted)) +
  labs(x='\nage', y=sprintf('%s\n', corticalmap),
       title=sprintf('partialR2=%s\nanovaPval=%s\n',
                     gam.results$partialRsq, gam.results$anova.smooth.pvalue)) +
  theme_classic() +
  theme(
    axis.text = element_text(size=12, family = "Arial", color = c("black")),
    axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
    axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
  theme(legend.position="none") +
  theme(aspect.ratio=1) +
  scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), limits = c(6,22), expand = c(0.05,.05)) +
  # xlim(6, 22) +
  ylim(lolim, hilim)

ggsave(paste(outpath, sprintf('%s_%s.svg', dtype, gamtype), sep = ""),
       dpi = 300,
       plot = last_plot())
}

if(gamtype == 'pfactor'){
  GAM.RESULTS <- gam.fit.linear(measure = 'metric', atlas = 'schaefer400',
                                dataset = 'all', region = corticalmap,
                                smooth_var = smooth_var,
                                linear_var = linear_var,
                                covariates = covars,
                                knots = 3, set_fx = FALSE)
  gam.results <- as.data.frame(GAM.RESULTS)

  ctx.predicted.metric <- gam.linear.predict(measure = 'metric',
                                             atlas = 'schaefer400',
                                             dataset = 'all',
                                             region = corticalmap,
                                             smooth_var = smooth_var,
                                             linear_var = linear_var,
                                             covariates = covars,
                                             knots = 3, set_fx = FALSE,
                                             increments = 200)

  # get predicted.smooth df from function output
  ctx.predicted.metric <- as.data.frame(ctx.predicted.metric[3])

  # plot p-factor
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
ggplot(data = metric.schaefer400.all, aes(x = p_factor_mcelroy_harmonized_all_samples,
                                          if (harmonized == TRUE){y = meanValHarmonized}
                                          else{y = meanVal})) +
  geom_point(aes(color = dataset), size = 2) + # color = "#115c25"
  geom_ribbon(data = ctx.predicted.metric, aes(x = p_factor_mcelroy_harmonized_all_samples,
                                               y = .fitted,
                                               ymin = .lower_ci, ymax = .upper_ci),
              alpha = .7, linetype = 0,) +
  geom_line(data = ctx.predicted.metric, aes(x = p_factor_mcelroy_harmonized_all_samples,
                                             y = .fitted)) +
  labs(x="\npfactor", y=sprintf('%s\n', corticalmap),
       title=sprintf('partialR2=%s\nanovaPval=%s\n',
                     gam.results$partialRsq, gam.results$anova.linear.pvalue)) +
  theme_classic() +
  theme(
    axis.text = element_text(size=12, family = "Arial", color = c("black")),
    axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
    axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
  theme(legend.position="none") +
  theme(aspect.ratio=1) +
  scale_x_continuous(breaks=c(-2, -1, 0, 1, 2, 3), limits = c(-2, 3), expand = c(0.05,0.05)) +
  # xlim(-2, 3) +
  ylim(lolim, hilim)

ggsave(paste(outpath, sprintf('%s_%s.svg', dtype, gamtype), sep = ""),
       dpi = 300,
       plot = last_plot())
}

##############################
# Section 2: Parcel-wise GAM Analyses
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

if(gamtype=='age'){
  gam.variable.schaefer <- matrix(data=NA, nrow=400, ncol=10)
  #for each schaefer region
  for(row in c(1:nrow(schaefer.regions))){
    region <- schaefer.regions$region[row]
    GAM.RESULTS <- gam.fit.smooth(measure = "metric", atlas = "schaefer400",
                                  dataset = "all",
                                  region = region, smooth_var = smooth_var,
                                  covariates = covars,
                                  knots = 3, set_fx = FALSE, stats_only = FALSE)
    #and append results to output df
    gam.variable.schaefer[row,] <- GAM.RESULTS}

  gam.variable.schaefer <- as.data.frame(gam.variable.schaefer)
  colnames(gam.variable.schaefer) <- c("region","GAM.variable.Fvalue","GAM.variable.pvalue",
                                       "GAM.variable.partialR2","Anova.variable.pvalue",
                                       "age.onsetchange", "age.peakchange",
                                       "minage.decrease","maxage.increase","age.maturation")
  cols = c(2:10)
  gam.variable.schaefer[,cols] = apply(gam.variable.schaefer[,cols], 2,
                                       function(x) as.numeric(as.character(x)))
}

if(gamtype=='pfactor'){
  gam.variable.schaefer <- matrix(data=NA, nrow=400, ncol=5)
  #for each schaefer region
  for(row in c(1:nrow(schaefer.regions))){
    region <- schaefer.regions$region[row]
    GAM.RESULTS <- gam.fit.linear(measure = "metric", atlas = "schaefer400",
                                  dataset = "all",
                                  region = region,
                                  smooth_var = smooth_var,
                                  linear_var = linear_var,
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
}

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

# # ranks
# rankR2 = rank(gam.variable.schaefer$GAM.variable.partialR2, ties.method = c("average"))
# gam.variable.schaefer$GAM.variable.rankR2 <- rankR2
# 
# nearestNeg <- max(gam.variable.schaefer$GAM.variable.partialR2[gam.variable.schaefer$GAM.variable.partialR2<0])
# nearestNegIdx <- which(gam.variable.schaefer$GAM.variable.partialR2 == nearestNeg)
# 
# nearestNegRank <- gam.variable.schaefer$GAM.variable.rankR2[nearestNegIdx]
# 
# gam.variable.schaefer$GAM.variable.rankR2 <- gam.variable.schaefer$GAM.variable.rankR2-(nearestNegRank+1)
# 
# maxval <- max(abs(gam.variable.schaefer$GAM.variable.rankR2))

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

###############################
# Section 3: Study-Specific Analyses
# Fit a gam for each study
###############################
# Iterate over each study (bhrc, ccnp, hbn, nki, pnc)
# Load data and re-fit GAMs on the same cortical metric
# Visualize predicted smooths overlaid for each study
# Save visualizations as SVG plots

# study-specific fits
data_labels <- c('bhrc', 'ccnp', 'hbn', 'nki', 'pnc') %>% as.data.frame() %>% set_names("data")

b=ggplot()
for(row in c(1:nrow(data_labels))){
  dataset <- data_labels$data[row]
  if(dataset=='bhrc'){ribboncolor <- '#F8766D'}
  if(dataset=='ccnp'){ribboncolor <- '#A3A500'}
  if(dataset=='hbn'){ribboncolor <- '#00BF7D'}
  if(dataset=='nki'){ribboncolor <- '#00B0F6'}
  if(dataset=='pnc'){ribboncolor <- '#E76BF3'}


  if (gamtype == 'age'){
    dtype <- sprintf('%s_df_%s_%s', dataset, metric, qc_version)
    outlabel <- sprintf('df_%s_%s', metric, qc_version)
    if (harmonized == TRUE){
      dtype <- sprintf('%s_df_%s_%s_harmonized', dataset, metric, qc_version)
      outlabel <- sprintf('df_%s_%s_harmonized', metric, qc_version)
    }
  }
  if (gamtype == 'pfactor'){
    dtype <- sprintf('%s_df_%s_%s_pfactor_filter', dataset, metric, qc_version)
    outlabel <- sprintf('df_%s_%s_pfactor_filter', metric, qc_version)
    if (harmonized == TRUE){
      dtype <- sprintf('%s_df_%s_%s_pfactor_filter_harmonized', dataset, metric, qc_version)
      outlabel <- sprintf('df_%s_%s_pfactor_filter_harmonized', metric, qc_version)
    }
  }

  if (gamtype == 'age'){smooth_var <- 'age'}
  if (gamtype == 'age'){covars <- 'sex + euler'}

  if (gamtype == 'pfactor'){smooth_var <- 'age'}
  if (gamtype == 'pfactor'){covars <- 'sex + euler'}
  if (gamtype == 'pfactor'){linear_var <- 'p_factor_mcelroy_harmonized_all_samples'}

  # prepare data for GAMs
  metric.schaefer400.all <- read.csv(paste(data_path, sprintf('%s.tsv', dtype),
                                           sep = ""),
                                     sep = '\t')
  # will use dataframe as a covariate
  metric.schaefer400.all$dataset <- as.factor(metric.schaefer400.all$study)
  # will use sex as a covariate
  metric.schaefer400.all$sex <- as.factor(metric.schaefer400.all$sex)

  if (controlformean == TRUE){
    if (harmonized == TRUE){
      meanmapcovar <- metric.schaefer400.all$meanValHarmonized
      covars <- 'sex + euler + meanmapcovar'}
    else if (harmonized == FALSE){
      meanmapcovar <- metric.schaefer400.all$meanVal
      covars <- 'sex + euler + meanmapcovar'}
  }

  if (gamtype == 'pfactor'){
    # fit gam and get fitted lines
    ctx.predicted.metric <- gam.linear.predict(measure = 'metric',
                                               atlas = 'schaefer400',
                                               dataset = 'all',
                                               region = corticalmap,
                                               smooth_var = smooth_var,
                                               linear_var = linear_var,
                                               covariates = covars,
                                               knots = 3, set_fx = FALSE,
                                               increments = 200)
    # get predicted.smooth df from function output
    ctx.predicted.metric <- as.data.frame(ctx.predicted.metric[3])

    # for p-factor
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

    b <- b +
      geom_ribbon(data = ctx.predicted.metric, aes(x = p_factor_mcelroy_harmonized_all_samples,
                                                   y = .fitted, ymin = .lower_ci, ymax = .upper_ci),
                  alpha = .3, linetype = 0, fill = c(ribboncolor)) +
      geom_line(data = ctx.predicted.metric, aes(x = p_factor_mcelroy_harmonized_all_samples,
                                                 y = .fitted), color = c(ribboncolor)) +
      labs(x='\npfactor', y=sprintf('%s\n', corticalmap)) +
      theme_classic() +
      theme(
        axis.text = element_text(size=12, family = "Arial", color = c("black")),
        axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
        axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
      theme(legend.position="none") +
      theme(aspect.ratio=1) +
      scale_x_continuous(breaks=c(-2, -1, 0, 1, 2, 3), limits=c(-2, 3), expand = c(0.05,.05)) +
      ylim(lolim, hilim)
  }

  if (gamtype == 'age'){
    # fit gam and get fitted lines
    ctx.predicted.metric <- gam.smooth.predict(measure = 'metric',
                                               atlas = 'schaefer400',
                                               dataset = 'all',
                                               region = corticalmap,
                                               smooth_var = smooth_var,
                                               covariates = covars,
                                               knots = 3, set_fx = FALSE,
                                               increments = 200)
    # get predicted.smooth df from function output
    ctx.predicted.metric <- as.data.frame(ctx.predicted.metric[3])

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

    b <- b +
      geom_ribbon(data = ctx.predicted.metric, aes(x = age,
                                                   y = .fitted, ymin = .lower_ci, ymax = .upper_ci),
                  alpha = .3, linetype = 0, fill = c(ribboncolor)) +
      geom_line(data = ctx.predicted.metric, aes(x = age,
                                                 y = .fitted), color = c(ribboncolor)) +
      labs(x='\nage', y=sprintf('%s\n', corticalmap)) +
      theme_classic() +
      theme(
        axis.text = element_text(size=12, family = "Arial", color = c("black")),
        axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
        axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
      theme(legend.position="none") +
      theme(aspect.ratio=1) +
      scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), limits = c(6,22), expand = c(0.05,.05)) +
      ylim(lolim, hilim)
  }

}
print(b)

ggsave(paste(outpath, sprintf('studyfits_%s_%s.svg', outlabel, gamtype), sep = ""),
       dpi = 300,
       plot = last_plot())

###############################
# Section 4: Site-Specific Analyses
# Updated with unique study_site labels from metric.schaefer400.all
###############################
# Iterate over each site
# Load data and re-fit GAMs on the same cortical metric
# Visualize predicted smooths overlaid for each study
# Save visualizations as SVG plots

# Define study sites and corresponding ribbon colors
study_site <- c(
  'PNC1'        = '#D11141',  # bright red
  'HBNsiteSI'   = '#00B159',  # bright green
  'HBNsiteRU'   = '#00AEDB',  # cyan
  'HBNsiteCBIC' = '#F37735',  # orange
  'HBNsiteCUNY' = '#FFC425',  # yellow
  'BHRC-1'      = '#8C564B',  # brown
  'BHRC-2'      = '#6A1B9A',  # purple
  'NKI'         = '#1F77B4',  # blue
  'Colornest'   = '#E377C2'   # pink
)

# Load combined dataset
if (gamtype == 'age') {
  dtype <- sprintf('combined_df_%s_%s', metric, qc_version)
  outlabel <- sprintf('df_%s_%s', metric, qc_version)
  if (harmonized == TRUE) {
    dtype <- sprintf('combined_df_%s_%s_harmonized', metric, qc_version)
    outlabel <- sprintf('df_%s_%s_harmonized', metric, qc_version)
  }
} else if (gamtype == 'pfactor') {
  dtype <- sprintf('combined_df_%s_%s_pfactor_filter', metric, qc_version)
  outlabel <- sprintf('df_%s_%s_pfactor_filter', metric, qc_version)
  if (harmonized == TRUE) {
    dtype <- sprintf('combined_df_%s_%s_pfactor_filter_harmonized', metric, qc_version)
    outlabel <- sprintf('df_%s_%s_pfactor_filter_harmonized', metric, qc_version)
  }
}

metric.schaefer400.all <- read.csv(paste(data_path, sprintf('%s.tsv', dtype), sep = ""), sep = '\t')
metric.schaefer400.all$sex <- as.factor(metric.schaefer400.all$sex)
metric.schaefer400.all$study_site <- as.factor(metric.schaefer400.all$study_site)

# Loop through each study site
b <- ggplot()
for (dataset in names(study_site)) {
  ribboncolor <- study_site[[dataset]]
  sitemetric.schaefer400.all <- metric.schaefer400.all %>% dplyr::filter(study_site == dataset)
  
  if (gamtype == 'age') {
    smooth_var <- 'age'
    covars <- 'sex + euler'
  } else if (gamtype == 'pfactor') {
    smooth_var <- 'age'
    covars <- 'sex + euler'
    linear_var <- 'p_factor_mcelroy_harmonized_all_samples'
  }
  
  if (controlformean == TRUE) {
    if (harmonized == TRUE) {
      meanmapcovar <- sitemetric.schaefer400.all$meanValHarmonized
    } else {
      meanmapcovar <- sitemetric.schaefer400.all$meanVal
    }
    covars <- 'sex + euler + meanmapcovar'
  }
  
  if (gamtype == 'pfactor') {
    ctx.predicted.metric <- gam.linear.predict(measure = 'sitemetric', 
                                               atlas = 'schaefer400', 
                                               dataset = 'all', 
                                               region = corticalmap, 
                                               smooth_var = smooth_var, 
                                               linear_var = linear_var, 
                                               covariates = covars, 
                                               knots = 3, 
                                               set_fx = FALSE, 
                                               increments = 200)
    ctx.predicted.metric <- as.data.frame(ctx.predicted.metric[3])
    
    if(metric == 'ct'){lolim <- 2.0; hilim <- 3.5}
    if(metric == 'sa'){lolim <- 110; hilim <- 700}
    if(metric == 'gv'){lolim <- 350; hilim <- 2200}
    if(metric == 'lgi'){lolim <- 2.0; hilim <- 4.0}
    
    b <- b +
      geom_ribbon(data = ctx.predicted.metric, aes(x = p_factor_mcelroy_harmonized_all_samples, 
                                                   y = .fitted, ymin = .lower_ci, ymax = .upper_ci), 
                  alpha = .3, linetype = 0, fill = ribboncolor) +
      geom_line(data = ctx.predicted.metric, aes(x = p_factor_mcelroy_harmonized_all_samples, 
                                                 y = .fitted), color = ribboncolor) +
      labs(x='\npfactor', y=sprintf('%s\n', corticalmap)) +
      theme_classic() +
      theme(axis.text = element_text(size=12, family = "Arial", color = "black"),
            axis.title.x = element_text(size=12, family ="Arial", color = "black"),
            axis.title.y = element_text(size=12, family ="Arial", color = "black"),
            legend.position="none",
            aspect.ratio=1) +
      scale_x_continuous(breaks=c(-2, -1, 0, 1, 2, 3), limits=c(-2, 3), expand = c(0.05,.05)) +
      ylim(lolim, hilim)
  }
  
  if (gamtype == 'age') {
    ctx.predicted.metric <- gam.smooth.predict(measure = 'sitemetric',
                                               atlas = 'schaefer400',
                                               dataset = 'all',
                                               region = corticalmap,
                                               smooth_var = smooth_var,
                                               covariates = covars,
                                               knots = 3, set_fx = FALSE,
                                               increments = 200)
    ctx.predicted.metric <- as.data.frame(ctx.predicted.metric[3])
    
    if(metric == 'ct'){lolim <- 2.0; hilim <- 3.5}
    if(metric == 'sa'){lolim <- 110; hilim <- 700}
    if(metric == 'gv'){lolim <- 350; hilim <- 2200}
    if(metric == 'lgi'){lolim <- 2.0; hilim <- 4.0}
    
    b <- b +
      geom_ribbon(data = ctx.predicted.metric, aes(x = age, y = .fitted, ymin = .lower_ci, ymax = .upper_ci), 
                  alpha = .3, linetype = 0, fill = ribboncolor) +
      geom_line(data = ctx.predicted.metric, aes(x = age, y = .fitted), color = ribboncolor) +
      labs(x='\nage', y=sprintf('%s\n', corticalmap)) +
      theme_classic() +
      theme(axis.text = element_text(size=12, family = "Arial", color = "black"),
            axis.title.x = element_text(size=12, family ="Arial", color = "black"),
            axis.title.y = element_text(size=12, family ="Arial", color = "black"),
            legend.position="none",
            aspect.ratio=1) +
      scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), limits = c(6,22), expand = c(0.05,.05)) +
      ylim(lolim, hilim)
  }
}

print(b)
ggsave(paste(outpath, sprintf('sitefits_%s_%s.svg', outlabel, gamtype), sep = ""), dpi = 300, plot = last_plot())
