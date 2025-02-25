
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
outpath <- paste(project_path, 'results/function/', sep = "")

# # for pfactor
# combined_df_withinbetween_fcrsn7_noqc_pfactor_filter
# combined_df_withinbetween_fcrsn7_artifact_pfactor_filter
# combined_df_withinbetween_fcrsn7_artifact_pfactor_filter_harmonized
# 
# # for age
# combined_df_withinbetween_fcrsn7_noqc
# combined_df_withinbetween_fcrsn7_artifact
# combined_df_withinbetween_fcrsn7_artifact_harmonized

dataset <- 'combined'
qc_version <- 'artifact' # 'artifact' or 'noqc'
gamtype <- 'pfactor' # 'age' or 'pfactor'
harmonized <- TRUE

if (gamtype == 'age'){
  dtype <- sprintf('%s_df_withinbetween_fcrsn7_%s', dataset, 
                   qc_version)
  if (harmonized == TRUE){
    dtype <- sprintf('%s_df_withinbetween_fcrsn7_%s_harmonized', dataset, qc_version)
  }
}
if (gamtype == 'pfactor'){
  dtype <- sprintf('%s_df_withinbetween_fcrsn7_%s_pfactor_filter', dataset, 
                   qc_version)
  if (harmonized == TRUE){
    dtype <- sprintf('%s_df_withinbetween_fcrsn7_%s_pfactor_filter_harmonized', 
                     dataset, qc_version)
  }
}

if (gamtype == 'age'){smooth_var <- 'age'}
if (gamtype == 'age'){covars <- 'sex + medianFD'}

if (gamtype == 'pfactor'){smooth_var <- 'age'}
if (gamtype == 'pfactor'){covars <- 'sex + medianFD'}
if (gamtype == 'pfactor'){linear_var <- 'p_factor_mcelroy_harmonized_all_samples'}

# prepare data for GAMs
netpair.schaefer400.all <- read.csv(paste(data_path, 
                                          sprintf('%s.tsv', dtype), 
                                          sep = ""), 
                                    sep = '\t')

# will use dataframe as a covariate
netpair.schaefer400.all$dataset <- as.factor(netpair.schaefer400.all$study)
# will use sex as a covariate
netpair.schaefer400.all$sex <- as.factor(netpair.schaefer400.all$sex)

##############################
# Fit GAMs for Network Pairs #
##############################
if(dataset == 'combined'){
#### Run gam.fit.smooth in all networks: with k=3
# #list of regions to run gam.fit.smooth function on below
# netpair_labels <- colnames(netpair.schaefer400.all)[0:28]
netpair_labels <- names(netpair.schaefer400.all[0:28]) %>% as.data.frame() %>% set_names("netpair")

# network connectivity GAMs
if(gamtype == 'age'){
  #empty matrix to save gam.fit output to
  gam.age <- matrix(data=NA, nrow=28, ncol=10)
  #for each network pair
  for(row in c(1:nrow(netpair_labels))){
    netpair <- netpair_labels$netpair[row]
    #run the gam.fit.smooth function
    GAM.RESULTS <- gam.fit.smooth(measure = "netpair", atlas = "schaefer400", 
                                  dataset = "all", region = netpair, 
                                  smooth_var = smooth_var, covariates = covars,
                                  knots = 3, set_fx = FALSE, stats_only = FALSE)
    #and append results to output df
    gam.age[row,] <- GAM.RESULTS}
  
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("netpair","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.partialR2",
                         "Anova.age.pvalue","age.onsetchange","age.peakchange",
                         "minage.decrease","maxage.increase","age.maturation")
  # pvalues
  pvalues = gam.age$GAM.age.pvalue
  pvaluesfdrs<-p.adjust(pvalues, method="BH")
  gam.age$GAM.age.pvaluefdr <- pvaluesfdrs
  
  pvalues = gam.age$Anova.age.pvalue
  pvaluesfdrs<-p.adjust(pvalues, method="BH")
  
  # nannumber = sum(pvaluesfdrs < 0.05)
  # pvaluesfdrs[pvaluesfdrs >= 0.05] <- NA
  gam.age$Anova.age.pvaluefdr <- pvaluesfdrs
  
  # make sure values are numerical
  cols = c(2:11)
  gam.age[,cols] = apply(gam.age[,cols], 2, 
                         function(x) as.numeric(as.character(x)))
  write.csv(gam.age, paste(outpath,
                           sprintf('csvFiles/%s_%s_statistics.csv', dtype, 
                                   gamtype),
                           sep = ""),
            row.names = F, quote = F)
  
  rm(gam.age)
  gc()
  }

if(gamtype == 'pfactor'){
  #empty matrix to save gam.fit output to
  gam.pfactor <- matrix(data=NA, nrow=28, ncol=5)
  #for each network pair
  for(row in c(1:nrow(netpair_labels))){
    netpair <- netpair_labels$netpair[row]
    #run the gam.fit.smooth function
    GAM.RESULTS <- gam.fit.linear(measure = "netpair", atlas = "schaefer400", 
                                  dataset = "all", region = netpair, 
                                  linear_var = linear_var, 
                                  smooth_var = smooth_var, 
                                  covariates = covars,
                                  knots = 3, set_fx = FALSE)
    #and append results to output df
    gam.pfactor[row,] <- GAM.RESULTS}
  
  gam.pfactor <- as.data.frame(gam.pfactor)
  colnames(gam.pfactor) <- c("netpair","GAM.pfactor.Fvalue","GAM.pfactor.pvalue",
                             "GAM.pfactor.partialR2",
                             "Anova.pfactor.pvalue")
  
  # pvalues
  pvalues = gam.pfactor$GAM.pfactor.pvalue
  pvaluesfdrs<-p.adjust(pvalues, method="BH")
  gam.pfactor$GAM.pfactor.pvaluefdr <- pvaluesfdrs
  
  pvalues = gam.pfactor$Anova.pfactor.pvalue
  pvaluesfdrs<-p.adjust(pvalues, method="BH")
  
  # nannumber = sum(pvaluesfdrs < 0.05)
  # pvaluesfdrs[pvaluesfdrs >= 0.05] <- NA
  gam.pfactor$Anova.pfactor.pvaluefdr <- pvaluesfdrs
  
  # make sure values are numerical
  cols = c(2:6)
  gam.pfactor[,cols] = apply(gam.pfactor[,cols], 2, 
                             function(x) as.numeric(as.character(x)))
  write.csv(gam.pfactor, paste(outpath,
                               sprintf('csvFiles/%s_%s_statistics.csv', 
                                       dtype, gamtype),
                               sep = ""),
            row.names = F, quote = F)
  rm(gam.pfactor)
  gc()
  }
}
###############################
# Fit a gam for certain network pairs as examples
###############################
#### PREDICT GAM SMOOTH FITTED VALUES ####
withinrsn <- FALSE
if(gamtype == 'pfactor'){
  networkpair <- 'Cont.Default'
  ctx.predicted.metric <- gam.linear.predict(measure = 'netpair', 
                                             atlas = 'schaefer400',
                                             dataset = 'all', 
                                             region = networkpair,
                                             smooth_var = smooth_var, 
                                             linear_var = linear_var,
                                             covariates = covars,
                                             knots = 3, set_fx = FALSE, 
                                             increments = 200)
  # get predicted.smooth df from function output
  ctx.predicted.metric <- as.data.frame(ctx.predicted.metric[3])
  
  # plot p-factor
  ggplot(data = netpair.schaefer400.all, aes(x = p_factor_mcelroy_harmonized_all_samples,
                                             y = Cont.Default)) +
    geom_point(aes(color = dataset), size = 2) + # color = "#115c25"
    geom_ribbon(data = ctx.predicted.metric, aes(x = p_factor_mcelroy_harmonized_all_samples,
                                                 y = .fitted,
                                                 ymin = .lower_ci, ymax = .upper_ci),
                alpha = .7, linetype = 0,) +
    geom_line(data = ctx.predicted.metric, aes(x = p_factor_mcelroy_harmonized_all_samples,
                                               y = .fitted)) +
    labs(x='\npfactor', y=sprintf('%s\n', networkpair)) +
    theme_classic() +
    theme(
      axis.text = element_text(size=12, family = "Arial", color = c("black")),
      axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
      axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
    theme(legend.position="none") +
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=c(-2, -1, 0, 1, 2, 3), limits=c(-2, 3), expand = c(0.05,.05)) +
    ylim(-0.25, 0.75) # between
  
  # ggsave(paste(outpath, sprintf('%s_%s_%s.png', dtype, gamtype, networkpair), sep = ""),
  #        dpi = 300,
  #        plot = last_plot())
  ggsave(paste(outpath, sprintf('%s_%s_%s.svg', dtype, gamtype, networkpair), sep = ""),
         dpi = 300,
         plot = last_plot())
}

if(gamtype == 'age'){
  if(withinrsn == TRUE){networkpair <- 'SalVentAttn.SalVentAttn'}
  else{networkpair <- 'Default.SalVentAttn'}
  ctx.predicted.metric <- gam.smooth.predict(measure = 'netpair', 
                                             atlas = 'schaefer400',
                                             dataset = 'all', 
                                             region = networkpair,
                                             smooth_var = smooth_var, 
                                             covariates = covars,
                                             knots = 3, set_fx = FALSE, 
                                             increments = 200)
  # get predicted.smooth df from function output
  ctx.predicted.metric <- as.data.frame(ctx.predicted.metric[3])

  # plot
  ggplot(data = netpair.schaefer400.all, aes(x = age, 
                                             if (withinrsn == TRUE){y = SalVentAttn.SalVentAttn}
                                             else{y = Default.SalVentAttn})) +
    geom_point(aes(color = dataset), size = 2) + # color = "#115c25"
    geom_ribbon(data = ctx.predicted.metric, aes(x = age, y = .fitted,
                                                 ymin = .lower_ci, ymax = .upper_ci),
                alpha = .7, linetype = 0,) +
    geom_line(data = ctx.predicted.metric, aes(x = age, y = .fitted)) +
    labs(x='\nage', y=sprintf('%s\n', networkpair)) +
    theme_classic() +
    theme(
      axis.text = element_text(size=12, family = "Arial", color = c("black")),
      axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
      axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
    theme(legend.position="none") +
    theme(aspect.ratio=1) +
    scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), limits = c(6,22), expand = c(0.05,.05)) +
    (if (withinrsn == TRUE){ylim(0.03, 0.8)}
     else{ylim(-0.40, 0.75)})
  
  # ggsave(paste(outpath, sprintf('%s_%s_%s.png', dtype, gamtype, networkpair), sep = ""),
  #        dpi = 300,
  #        plot = last_plot())
  ggsave(paste(outpath, sprintf('%s_%s_%s.svg', dtype, gamtype, networkpair), sep = ""),
         dpi = 300,
         plot = last_plot())
}

###############################
# Fit a gam for certain network pairs as examples
###############################
# study-specific fits
data_labels <- c('bhrc', 'ccnp', 'hbn', 'nki', 'pnc') %>% as.data.frame() %>% set_names("data")
withinrsn <- FALSE

b=ggplot()
for(row in c(1:nrow(data_labels))){ 
  dataset <- data_labels$data[row]
  if(dataset=='bhrc'){ribboncolor <- '#F8766D'}
  if(dataset=='ccnp'){ribboncolor <- '#A3A500'}
  if(dataset=='hbn'){ribboncolor <- '#00BF7D'}
  if(dataset=='nki'){ribboncolor <- '#00B0F6'}
  if(dataset=='pnc'){ribboncolor <- '#E76BF3'}
  
  if (gamtype == 'age'){
    dtype <- sprintf('%s_df_withinbetween_fcrsn7_%s', dataset, 
                     qc_version)
    outlabel <- sprintf('df_withinbetween_fcrsn7_%s', qc_version)
    if (harmonized == TRUE){
      dtype <- sprintf('%s_df_withinbetween_fcrsn7_%s_harmonized', 
                       dataset, qc_version)
      outlabel <- sprintf('df_withinbetween_fcrsn7_%s_harmonized', 
                          qc_version)
    }
  }
  if (gamtype == 'pfactor'){
    dtype <- sprintf('%s_df_withinbetween_fcrsn7_%s_pfactor_filter', dataset, 
                     qc_version)
    outlabel <- sprintf('df_withinbetween_fcrsn7_%s_pfactor_filter', qc_version)
    if (harmonized == TRUE){
      dtype <- sprintf('%s_df_withinbetween_fcrsn7_%s_pfactor_filter_harmonized', 
                       dataset, qc_version)
      outlabel <- sprintf('df_withinbetween_fcrsn7_%s_pfactor_filter_harmonized', 
                          qc_version)
    }
  }
  
  if (gamtype == 'age'){smooth_var <- 'age'}
  if (gamtype == 'age'){covars <- 'sex + medianFD'}
  
  if (gamtype == 'pfactor'){smooth_var <- 'age'}
  if (gamtype == 'pfactor'){covars <- 'sex + medianFD'}
  if (gamtype == 'pfactor'){linear_var <- 'p_factor_mcelroy_harmonized_all_samples'}

  # prepare data for GAMs
  netpair.schaefer400.all <- read.csv(paste(data_path, 
                                            sprintf('%s.tsv', dtype), 
                                            sep = ""), 
                                      sep = '\t')
  
  # will use dataframe as a covariate
  netpair.schaefer400.all$dataset <- as.factor(netpair.schaefer400.all$study)
  # will use sex as a covariate
  netpair.schaefer400.all$sex <- as.factor(netpair.schaefer400.all$sex)

  if (gamtype == 'pfactor'){
    # fit gam and get fitted lines
    networkpair <- 'Cont.Default'
    ctx.predicted.metric <- gam.linear.predict(measure = 'netpair', 
                                               atlas = 'schaefer400',
                                               dataset = 'all', 
                                               region = networkpair,
                                               smooth_var = smooth_var, 
                                               linear_var = linear_var,
                                               covariates = covars,
                                               knots = 3, set_fx = FALSE, 
                                               increments = 200)
    # get predicted.smooth df from function output
    ctx.predicted.metric <- as.data.frame(ctx.predicted.metric[3])
    
    b <- b +
      geom_ribbon(data = ctx.predicted.metric, aes(x = p_factor_mcelroy_harmonized_all_samples,
                                                   y = .fitted, ymin = .lower_ci, ymax = .upper_ci),
                  alpha = .3, linetype = 0, fill = c(ribboncolor)) +
      geom_line(data = ctx.predicted.metric, aes(x = p_factor_mcelroy_harmonized_all_samples,
                                                 y = .fitted), color = c(ribboncolor)) +
      labs(x='\npfactor', y=sprintf('%s\n', networkpair)) +
      theme_classic() +
      theme(
        axis.text = element_text(size=12, family = "Arial", color = c("black")),
        axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
        axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
      theme(legend.position="none") +
      theme(aspect.ratio=1) +
      scale_x_continuous(breaks=c(-2, -1, 0, 1, 2, 3), limits=c(-2, 3), expand = c(0.05,.05)) +
      ylim(-0.25, 0.75) # between
  }
  
  if (gamtype == 'age'){
    # fit gam and get fitted lines
    if(withinrsn == TRUE){networkpair <- 'SalVentAttn.SalVentAttn'}
    else{networkpair <- 'Default.SalVentAttn'}
    ctx.predicted.metric <- gam.smooth.predict(measure = 'netpair', 
                                               atlas = 'schaefer400',
                                               dataset = 'all', 
                                               region = networkpair,
                                               smooth_var = smooth_var, 
                                               covariates = covars,
                                               knots = 3, set_fx = FALSE, 
                                               increments = 200)
    # get predicted.smooth df from function output
    ctx.predicted.metric <- as.data.frame(ctx.predicted.metric[3])
    
    b <- b +
      geom_ribbon(data = ctx.predicted.metric, aes(x = age,
                                                   y = .fitted, ymin = .lower_ci, ymax = .upper_ci),
                  alpha = .3, linetype = 0, fill = c(ribboncolor)) +
      geom_line(data = ctx.predicted.metric, aes(x = age,
                                                 y = .fitted), color = c(ribboncolor)) +
      labs(x='\nage', y=sprintf('%s\n', networkpair)) +
      theme_classic() +
      theme(
        axis.text = element_text(size=12, family = "Arial", color = c("black")),
        axis.title.x = element_text(size=12, family ="Arial", color = c("black")),
        axis.title.y = element_text(size=12, family ="Arial", color = c("black"))) +
      theme(legend.position="none") +
      theme(aspect.ratio=1) +
      scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22), limits = c(6,22), expand = c(0.05,.05)) +
      (if(withinrsn == TRUE){ylim(0.03, 0.8)}
      else{ylim(-0.40, 0.75)})
  }
  
}
print(b)

# ggsave(paste(outpath, sprintf('studyfits_%s_%s_%s.png', outlabel, gamtype, networkpair), sep = ""),
#        dpi = 300,
#        plot = last_plot())
ggsave(paste(outpath, sprintf('studyfits_%s_%s_%s.svg', outlabel, gamtype, networkpair), sep = ""),
       dpi = 300,
       plot = last_plot())


# # Estimate GAM smooths based on model-predicted data and save out predicted y data**   
# gam.smooths <- matrix(data=NA, ncol=7) #empty matrix to save gam.predsmooth fits to
# colnames(gam.smooths) <- c("age","fit","se.fit","selo","sehi","index","label")
# 
# for(row in c(1:length(netpair_labels))){ #for each region
#   netpair <- netpair_labels[row] 
#   GAM.SMOOTH <- gam.predsmooth(measure = "netpair", atlas = "schaefer200", conn_type = "orig", region = netpair, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.predsmooth function
#   
#   preddata <- as.data.frame(GAM.SMOOTH[3]) #get predicted.smooth df from function output
#   preddata$index <- rep(x=row, 1000) #region index
#   preddata$label <- rep(x=GAM.SMOOTH[1], 1000) #label
#   gam.smooths <- rbind(gam.smooths, preddata)
#   
#   
# }
# gam.smooths <- gam.smooths[-1,] #remove empty initialization row
# gam.smooths$label <- as.character(gam.smooths$label)
# 
# 
# write.csv(gam.smooths, sprintf("%1$sGAMsmoothfits.networkpair.csv", gam_dir), row.names = F, quote = F)
