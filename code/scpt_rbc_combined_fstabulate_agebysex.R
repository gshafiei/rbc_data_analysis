
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
controlformean <- FALSE # whether to include mean map as covariate
controlforTIV <- FALSE # whether to include mean map as covariate

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

metric.schaefer400.all$oSex <- factor(metric.schaefer400.all$sex, 
                                     levels = c("Male", "Female"), 
                                     ordered = TRUE)

if (controlformean == TRUE){
  if (harmonized == TRUE){
    meanmapcovar <- metric.schaefer400.all$meanValHarmonized
    covars <- 'oSex + euler + meanmapcovar'}
  else if (harmonized == FALSE){
    meanmapcovar <- metric.schaefer400.all$meanVal
    covars <- 'oSex + euler + meanmapcovar'}
} else if (controlforTIV == TRUE){
  tivmapcovar <- metric.schaefer400.all$eTIV
  covars <- 'euler + tivmapcovar'
} else if (controlformean == FALSE && controlforTIV == FALSE){
  covars <- 'oSex + euler'
}

smooth_var <- 'age'
model_var <- 'oSex'

# ###############################
# # Section 1: Overall Fit for Mean Value
# # Fit a gam for combined mean values
# ###############################
# # Fit a GAM using mean value (e.g., meanVal or meanValHarmonized)
# # gam.fit.smooth() for age-related analysis: age as smooth term
# # gam.fit.linear() for linear p-factor modeling
# # Predict fitted values with confidence intervals using custom predict functions
# # Plot curve with ggplot2 and save as SVG
# 
# #### PREDICT GAM SMOOTH FITTED VALUES ####
# 
# GAM.RESULTS <- gam.factorsmooth.interaction(measure = 'metric',
#                                             atlas = 'schaefer400',
#                                             dataset = 'all',
#                                             region = corticalmap,
#                                             smooth_var = smooth_var,
#                                             int_var = model_var,
#                                             covariates = covars,
#                                             knots = 3, set_fx = FALSE)
# gam.results <- as.data.frame(GAM.RESULTS)
# 
# # predict smooth fitted values for 'Female'
# ctx.predicted.metric.F <- gam.smooth.predict.covariateinteraction(measure = 'metric',
#                                                                   atlas = 'schaefer400',
#                                                                   dataset = 'all',
#                                                                   region = corticalmap,
#                                                                   smooth_var = smooth_var,
#                                                                   int_var = model_var,
#                                                                   int_var.predict = 'Female',
#                                                                   covariates = covars,
#                                                                   knots = 3,
#                                                                   set_fx = FALSE,
#                                                                   increments = 200)
# ctx.predicted.metric.F$sex <- 'Female'
# 
# # predict smooth fitted values for 'Male'
# ctx.predicted.metric.M <- gam.smooth.predict.covariateinteraction(measure = 'metric',
#                                                                   atlas = 'schaefer400',
#                                                                   dataset = 'all',
#                                                                   region = corticalmap,
#                                                                   smooth_var = smooth_var,
#                                                                   int_var = model_var,
#                                                                   int_var.predict = 'Male',
#                                                                   covariates = covars,
#                                                                   knots = 3,
#                                                                   set_fx = FALSE,
#                                                                   increments = 200)
# ctx.predicted.metric.M$sex <- 'Male'
# 
# # combine
# ctx.predicted.metric <- bind_rows(ctx.predicted.metric.F, ctx.predicted.metric.M)
# 
# 
# # plot age
# if(metric == 'ct'){
#   lolim <- 2.0
#   hilim <- 3.5
# } # ct
# if(metric == 'sa'){
#   lolim <- 110
#   hilim <- 700
# } # sa
# if(metric == 'gv'){
#   lolim <- 350
#   hilim <- 2200
# } # gv
# if(metric == 'lgi'){
#   lolim <- 2.0
#   hilim <- 4.0
# } # lgi
# 
# ggplot(metric.schaefer400.all, aes(x = age,
#                                    if (harmonized == TRUE){y = meanValHarmonized}
#                                    else{y = meanVal})) +
#   geom_point(aes(color = sex), alpha = 0.25, size = 1) +
#     geom_ribbon(data = predicted.metric,
#               aes(x = age, y = .fitted, ymin = .lower_ci, ymax = .upper_ci, fill = sex),
#               alpha = 0.5, inherit.aes = FALSE) +
#     geom_line(data = predicted.metric,
#             aes(x = age, y = .fitted, color = sex),
#             size = 1.2, inherit.aes = FALSE) +
#   labs(x='\nage', y=sprintf('%s\n', corticalmap),
#        title=sprintf('F-stat=%s\npval=%s\n',
#                      gam.results$gam.int.F, gam.results$gam.int.pvalue),
#        color = "Sex", fill = "Sex") +
#   theme_classic() +
#   theme(
#     axis.text = element_text(size=12, family = "Arial", color = c("black")),
#     axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
#     axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
#   theme(legend.position="right") +
#   theme(aspect.ratio=1) +
#   scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), limits = c(6,22), expand = c(0.05,.05)) +
#   ylim(lolim, hilim)
# 
# ggsave(paste(outpath, sprintf('%s_%s.svg', dtype, gamtype), sep = ""),
#        dpi = 300,
#        plot = last_plot())

##############################
# Section 2: Parcel-wise GAM Analyses
##############################
# Loop through 400 Schaefer parcels and run region-specific GAMs
# Collect GAM stats (F, p)
# Save output as .csv

# Merge with SA rank and adjust p-values (FDR)
# Plot brain maps:
#   - F values
#   - FDR-significant p-values
#   - F-based rank and rank significance maps

if(dataset == 'combined'){
  #### Region-wise GAM Statistics and Derivative-based Temporal Developmental Properties ####
  #list of regions to run gam.fit.smooth function on below
  schaefer.regions <- names(metric.schaefer400.all[0:400]) %>% as.data.frame() %>% set_names("region")
  
  # fit GAMs for sex
  gam.variable.schaefer <- matrix(data=NA, nrow=400, ncol=3)
  # also estimate a sign for F-stat by predicting outputs for F and M and finding the average difference
  gam.variable.schaefer.sign <- matrix(data=NA, nrow=400, ncol=1)
  
  # Create a PDF to store all plots
  pdf(paste(outpath, sprintf('%s_%s_schaefer400_gam_plots.pdf', dtype, gamtype), sep = ""), 
      width = 8, height = 10)
  # par(mfrow = c(5, 4), mar = c(4, 4, 2, 1))  # 20 plots per page
  
  #for each schaefer region
  for(row in c(1:nrow(schaefer.regions))){
    region <- schaefer.regions$region[row]
    GAM.RESULTS <- gam.factorsmooth.interaction(measure = 'metric', 
                                                atlas = 'schaefer400',
                                                dataset = 'all', 
                                                region = region,
                                                smooth_var = smooth_var, 
                                                int_var = model_var,
                                                covariates = covars,
                                                knots = 3, set_fx = FALSE)
    #and append results to output df
    gam.variable.schaefer[row,] <- GAM.RESULTS
  
    # predict smooth fitted values for 'Female'
    predicted.metric.F <- gam.smooth.predict.covariateinteraction(measure = 'metric',
                                                                  atlas = 'schaefer400',
                                                                  dataset = 'all',
                                                                  region = region,
                                                                  smooth_var = smooth_var,
                                                                  int_var = model_var,
                                                                  int_var.predict = 'Female',
                                                                  covariates = covars,
                                                                  knots = 3, 
                                                                  set_fx = FALSE,
                                                                  increments = 200)
    predicted.metric.F$sex <- 'Female'
    
    # predict smooth fitted values for 'Male'
    predicted.metric.M <- gam.smooth.predict.covariateinteraction(measure = 'metric',
                                                                  atlas = 'schaefer400',
                                                                  dataset = 'all',
                                                                  region = region,
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
    gam.variable.schaefer.sign[row,] <- mean_diff
    
    # combine
    predicted.metric <- bind_rows(predicted.metric.F, predicted.metric.M)
    
    # plot age
    p <- ggplot(metric.schaefer400.all, aes(x = age,
                                       y = .data[[region]])) +
      geom_point(aes(color = sex), alpha = 0.25, size = 1) +
      geom_ribbon(data = predicted.metric, 
                  aes(x = age, y = .fitted, ymin = .lower_ci, ymax = .upper_ci, fill = sex), 
                  alpha = 0.5, inherit.aes = FALSE) +
      geom_line(data = predicted.metric, 
                aes(x = age, y = .fitted, color = sex), 
                size = 1.2, inherit.aes = FALSE) +
      labs(x='\nage', y=sprintf('%s\n', region),
           title=sprintf('F-stat=%s\n',
                         as.numeric(GAM.RESULTS[2]) * sign(mean_diff)),
           color = "Sex", fill = "Sex") +
      theme_classic() +
      theme(
        axis.text = element_text(size=12, color = c("black")), # family = "Arial",
        axis.title.x = element_text(size=12, color = c("black")), # family = "Arial",
        axis.title.y = element_text(size=12, color = c("black"))) + # family = "Arial",
      theme(legend.position="right") +
      theme(aspect.ratio=1) +
      scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), limits = c(6,22), expand = c(0.05,.05))
    print(p)  # Send to PDF
  }
  dev.off()  # Close the PDF
  
 
  # prepare output csv results
  gam.variable.schaefer <- as.data.frame(gam.variable.schaefer)
  colnames(gam.variable.schaefer) <- c("region","GAM.variable.Fvalue","GAM.variable.pvalue")
  cols = c(2:3)
  gam.variable.schaefer[,cols] = apply(gam.variable.schaefer[,cols], 2,
                                       function(x) as.numeric(as.character(x)))

  # add sign
  gam.variable.schaefer.sign <- as.data.frame(gam.variable.schaefer.sign)
  colnames(gam.variable.schaefer.sign) <- c("meandiff")
  
  signed_F <- gam.variable.schaefer$GAM.variable.Fvalue * sign(gam.variable.schaefer.sign$meandiff)
  gam.variable.schaefer$GAM.variable.Fvalue.signed <- signed_F
  
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
  
  csvF <- data.frame(gam.variable.schaefer$region)
  csvF$Fvalue <- gam.variable.schaefer$GAM.variable.Fvalue
  csvF$signedFvalue <- gam.variable.schaefer$GAM.variable.Fvalue.signed
  # pvlues
  # GAMs
  pvalues = gam.variable.schaefer$GAM.variable.pvalue
  GAMpvaluesfdrs<-p.adjust(pvalues, method="BH")
  
  csvF$gamPvaluefdr <- GAMpvaluesfdrs
  outputPath <- paste(outpath, sprintf('csvFiles/%s_%s_F.csv', dtype, gamtype),
                      sep="")
  write.csv(csvF, outputPath, row.names=FALSE)
  
  # Effect size
  # brain
  maxval <- max(abs(gam.variable.schaefer$GAM.variable.Fvalue.signed))
  
  ggseg(.data = gam.variable.schaefer, atlas = "schaefer7_400",
        mapping=aes(fill = GAM.variable.Fvalue.signed, colour=I("#e9ecef"),
                    size=I(.03)), position = c("stacked")) +
    theme_void() +
    paletteer::scale_fill_paletteer_c("pals::warmcool",
                                      na.value="transparent", direction = -1,
                                      limits = c(-maxval, maxval),
                                      oob = squish)
  
  ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_Fvalue_signed.svg', dtype,
                                           gamtype),
                          sep = ""),
         dpi = 300, width = 3 , height = 2)
  
  
  # Replace GAM.variable.Fvalue.signed with ranks robustly
  # (1) Rank F values
  gam.variable.schaefer$GAM.variable.rankF <- rank(
    gam.variable.schaefer$GAM.variable.Fvalue.signed, 
    ties.method = "average"
  )
  
  # (2) Find the nearest negative value, if it exists
  neg_vals <- gam.variable.schaefer$GAM.variable.Fvalue.signed[
    gam.variable.schaefer$GAM.variable.Fvalue.signed < 0
  ]
  
  if (length(neg_vals) > 0) {
    # (3) Center ranks around the largest negative value (i.e., nearest to 0)
    nearestNeg <- max(neg_vals)
    nearestNegIdx <- which(gam.variable.schaefer$GAM.variable.Fvalue.signed == nearestNeg)
    nearestNegRank <- gam.variable.schaefer$GAM.variable.rankF[nearestNegIdx]
    
    gam.variable.schaefer$GAM.variable.rankF <- 
      gam.variable.schaefer$GAM.variable.rankF - (nearestNegRank + 1)
    
    centered <- TRUE
  } else {
    # No negatives — do not shift the ranks
    centered <- FALSE
  }
  
  # (4) Max absolute value (for plotting scale)
  maxval <- max(abs(gam.variable.schaefer$GAM.variable.rankF), na.rm = TRUE)
  
  
  # brain ranks
  ggseg(.data = gam.variable.schaefer, atlas = "schaefer7_400",
        mapping=aes(fill = GAM.variable.rankF, colour=I("#e9ecef"), size=I(.03)),
        position = c("stacked")) + theme_void() +
    paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent",
                                      direction = -1,
                                      limits = c(-maxval, maxval),
                                      oob = squish)
  
  
  ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_rank_Fvalue_signed.svg',
                                           dtype, gamtype),
                          sep = ""),
         dpi = 300, width = 3 , height = 2)
  
  # brain pvalues
  # GAM
  pvalues = gam.variable.schaefer$GAM.variable.pvalue
  pvaluesfdrs<-p.adjust(pvalues, method="BH")
  
  GAMsignumber = sum(pvaluesfdrs < 0.05, na.rm=TRUE)
  pvaluesfdrs[pvaluesfdrs >= 0.05] <- NA
  gam.variable.schaefer$GAM.variable.pvaluefdr <- pvaluesfdrs
  
  ggseg(.data = gam.variable.schaefer, atlas = "schaefer7_400",
        mapping=aes(fill = GAM.variable.pvaluefdr, colour=I("#e9ecef"), size=I(.03)),
        position = c("stacked")) + theme_void() + ggtitle(GAMsignumber) +
    paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent",
                                      direction = -1,
                                      limits = c(0, 0.05),
                                      # limits = c(min(metric.regional.statistics$GAM.smooth.pvalue),
                                      #            max(metric.regional.statistics$GAM.smooth.pvalue)),
                                      oob = squish)
  
  ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_pval_Fvalue_signed.svg',
                                           dtype, gamtype),
                          sep = ""),
         dpi = 300, width = 3 , height = 2)
  
  # significant ranks
  pvalues = gam.variable.schaefer$GAM.variable.pvalue
  pvaluesfdrs <- p.adjust(pvalues, method="BH")
  rankFsig <- gam.variable.schaefer$GAM.variable.rankF
  rankFsig[(pvaluesfdrs >= 0.05)] <- NA
  gam.variable.schaefer$GAM.variable.rankFsig <- rankFsig
  
  maxval <- max(abs(gam.variable.schaefer$GAM.variable.rankFsig), na.rm=T)
  
  # brain significant ranks
  ggseg(.data = gam.variable.schaefer, atlas = "schaefer7_400",
        mapping=aes(fill = GAM.variable.rankFsig, colour=I("#e9ecef"), size=I(.03)),
        position = c("stacked")) + theme_void() + ggtitle(GAMsignumber) +
    paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent",
                                      direction = -1,
                                      limits = c(-maxval, maxval),
                                      oob = squish)
  
  ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_ranksig_Fvalue_signed.svg',
                                           dtype, gamtype),
                          sep = ""),
         dpi = 300, width = 3 , height = 2)
}
