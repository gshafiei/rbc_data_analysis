
rm(list = ls())

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

source("/Users/gshafiei/Desktop/RBC/code/func_GAM_rbc.R")

project_path <- '/Users/gshafiei/Desktop/RBC/'
data_path <- paste(project_path, 'data/dataR/', sep = "")
outpath <- paste(project_path, 'results/structure/', sep = "")

# # for pfactor
# combined_df_ct_noqc_pfactor_filter
# combined_df_ct_artifact_pfactor_filter
# combined_df_ct_artifact_pfactor_filter_harmonized
# 
# # for age
# combined_df_ct_noqc
# combined_df_ct_artifact
# combined_df_ct_artifact_harmonized

dataset <- 'combined'
qc_version <- 'artifact' # 'artifact' or 'noqc'
gamtype <- 'pfactor' # 'age' or 'pfactor'
harmonized <- TRUE
controlformean <- FALSE

metric <- 'lgi'
corticalmap <- 'meanVal'
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

###############################
# Fit a gam for combined mean values
###############################
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

# ggsave(paste(outpath, sprintf('%s_%s.png', dtype, gamtype), sep = ""),
#        dpi = 300,
#        plot = last_plot())
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

# ggsave(paste(outpath, sprintf('%s_%s.png', dtype, gamtype), sep = ""),
#        dpi = 300,
#        plot = last_plot())
ggsave(paste(outpath, sprintf('%s_%s.svg', dtype, gamtype), sep = ""),
       dpi = 300,
       plot = last_plot())
}




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

# ggsave(filename = paste(outpath, sprintf('%s_%s_histogram_partialR2.png', dtype, 
#                                          gamtype), 
#                         sep = ""), 
#        dpi = 300, width = 3 , height = 2)
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
# grDevices::RdPu
# grDevices::PinkYl
# ggthemes::Red-Blue Diverging
# grDevices::Spectral
# pals::ocean.matter
# pals::coolwarm
# ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_partialR2.png', dtype, 
#                                          gamtype), 
#                         sep = ""), 
#        dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_partialR2.svg', dtype, 
                                         gamtype), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)

# ranks
rankR2 = rank(gam.variable.schaefer$GAM.variable.partialR2, ties.method = c("average"))
gam.variable.schaefer$GAM.variable.rankR2 <- rankR2

nearestNeg <- max(gam.variable.schaefer$GAM.variable.partialR2[gam.variable.schaefer$GAM.variable.partialR2<0])
nearestNegIdx <- which(gam.variable.schaefer$GAM.variable.partialR2 == nearestNeg)

nearestNegRank <- gam.variable.schaefer$GAM.variable.rankR2[nearestNegIdx]

gam.variable.schaefer$GAM.variable.rankR2 <- gam.variable.schaefer$GAM.variable.rankR2-(nearestNegRank+1)

maxval <- max(abs(gam.variable.schaefer$GAM.variable.rankR2))

# brain ranks
ggseg(.data = gam.variable.schaefer, atlas = "schaefer7_400", 
      mapping=aes(fill = GAM.variable.rankR2, colour=I("#e9ecef"), size=I(.03)), 
      position = c("stacked")) + theme_void() + 
  paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent", 
                                    direction = -1, 
                                    limits = c(-maxval, maxval), 
                                    oob = squish) 

# ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_rank_partialR2.png', 
#                                          dtype, gamtype), 
#                         sep = ""), 
#        dpi = 300, width = 3 , height = 2)
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

# nannumber = sum(metric.regional.statistics$GAM.smooth.pvalue < 0.05)
# metric.regional.statistics$GAM.smooth.pvalue[metric.regional.statistics$GAM.smooth.pvalue >= 0.05] <- NA

ggseg(.data = gam.variable.schaefer, atlas = "schaefer7_400", 
      mapping=aes(fill = Anova.variable.pvaluefdr, colour=I("#e9ecef"), size=I(.03)), 
      position = c("stacked")) + theme_void() + ggtitle(Anovasignumber) +
  paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent", 
                                    direction = -1, 
                                    limits = c(0, 0.05), 
                                    # limits = c(min(metric.regional.statistics$GAM.smooth.pvalue), 
                                    #            max(metric.regional.statistics$GAM.smooth.pvalue)), 
                                    oob = squish) 

# ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_pval_partialR2.png', 
#                                          dtype, gamtype), 
#                         sep = ""), 
#        dpi = 300, width = 3 , height = 2)
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

# ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_ranksig_partialR2.png', 
#                                          dtype, gamtype), 
#                         sep = ""), 
#        dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_ranksig_partialR2.svg', 
                                         dtype, gamtype), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
}

###############################
# Fit a gam for each study
###############################
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

# ggsave(paste(outpath, sprintf('studyfits_%s_%s.png', outlabel, gamtype), sep = ""),
#        dpi = 300,
#        plot = last_plot())
ggsave(paste(outpath, sprintf('studyfits_%s_%s.svg', outlabel, gamtype), sep = ""),
       dpi = 300,
       plot = last_plot())

# # brain sig
# pvalues = gam.age.schaefer$GAM.age.pvalue
# pvaluesfdrs <- p.adjust(pvalues, method="BH")
# allR2sig <- gam.age.schaefer$GAM.age.partialR2
# allR2sig[(pvaluesfdrs >= 0.05)] <- NA
# gam.age.schaefer$GAM.age.partialR2sig <- allR2sig
# 
# ggseg(.data = gam.age.schaefer, atlas = "schaefer7_400", 
#       mapping=aes(fill = GAM.age.partialR2sig, colour=I("#e9ecef"), 
#                   size=I(.03)), position = c("stacked")) + 
#   theme_void() + 
#   paletteer::scale_fill_paletteer_c("pals::ocean.matter", 
#                                     na.value="transparent", direction = -1, 
#                                     limits = c(min(gam.age.schaefer$GAM.age.partialR2), 
#                                                max(gam.age.schaefer$GAM.age.partialR2)), 
#                                     oob = squish) 
# 
# ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_sig_partialR2_%s_%s.png', dtype, smoothterm), 
#                         sep = ""), 
#        dpi = 300, width = 3 , height = 2)
# ggsave(filename = paste(outpath, sprintf('%s_%s_brainmap_sig_partialR2_%s_%s.svg', dtype, smoothterm), 
#                         sep = ""), 
#        dpi = 300, width = 3 , height = 2)

# # SA rank comparison
# rho = cor.test(gam.age.schaefer$GAM.age.partialR2, gam.age.schaefer$SA.rank, method = c("spearman"))
# 
# ggplot(gam.age.schaefer, aes(x = SA.rank, y = GAM.age.partialR2, fill = SA.rank)) + 
#   geom_point(aes(color = SA.rank), shape = 21, size = 2) +
#   scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", 
#                        guide = "colourbar", aesthetics = "fill", name = NULL, 
#                        midpoint = 200) +
#   scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", 
#                        guide = "colourbar", aesthetics = "color", name = NULL, 
#                        midpoint = 200) +
#   labs(x="\nS-A rank", y="partial R2\n") +
#   ggtitle(rho$estimate) +
#   geom_smooth(method = 'lm', se = TRUE, fill = alpha(c("gray70"),.7), col = "black", linewidth = 1) +
#   theme_classic() + 
#   theme(legend.position = "none") +
#   theme(axis.text = element_text(size = 12, family = "Arial", color = c("black")), 
#         axis.title = element_text(size = 12, family = "Arial", color = c("black")))
# 
# ggsave(filename = paste(outpath, sprintf('%s_%s_SArank_partialR2.png', dtype, smoothterm), 
#                         sep = ""), 
#        plot = last_plot())
# ggsave(filename = paste(outpath, sprintf('%s_%s_SArank_partialR2.svg', dtype, smoothterm), 
#                         sep = ""), 
#        plot = last_plot())
