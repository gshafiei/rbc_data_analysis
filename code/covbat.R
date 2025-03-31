# Load necessary libraries
library(ComBatFamily)
library(dplyr)
library(mgcv)

#################
# structural data
#################
# # for pfactor
# combined_df_ct_artifact_pfactor_filter

# # for age
# combined_df_ct_artifact

# Define paths and parameters for structural data harmonization
project_path <- '/Users/gshafiei/Desktop/RBC/'
data_path <- paste(project_path, 'data/dataR/', sep = "")

metric <- 'lgi'  # can also be 'gv', 'ct', or 'sa'
# we only want to harmonize 'artifact' version
qc_version <- 'artifact'  # can be 'noqc', 'artifact'
pfactor_include <- 'yes' # whether to include p-factor in covariates

# Construct the input filename and read the data
if (pfactor_include == 'yes'){
  dtype <- sprintf('%s_%s_pfactor_filter', metric, qc_version)
}
if (pfactor_include == 'no'){
  dtype <- sprintf('%s_%s', metric, qc_version)
}

combined_data <- read.csv(paste(data_path, sprintf('combined_df_%s.tsv', dtype),
                                sep = ""),
                          sep = '\t')

# Extract covariates and features
age_vec <- combined_data$age
sex_vec <- as.factor(combined_data$sex)
euler_vec <- combined_data$euler
if (pfactor_include == 'yes'){
pfactor_vec <- combined_data$p_factor_mcelroy_harmonized_all_samples
}

structural_data = combined_data[0:400] # assumes features are in first 400 columns

# Prepare covariate dataframe
# Harmonize using covfam with GAM (including nonlinear age effects)
if (pfactor_include == 'yes'){
  covar_df <- bind_cols(combined_data$participant_id,
                        as.numeric(age_vec),
                        as.factor(sex_vec),
                        as.numeric(euler_vec),
                        as.numeric(pfactor_vec))
  covar_df <- dplyr::rename(covar_df,
                            participant_id=...1,
                            age = ...2,
                            sex = ...3,
                            euler = ...4,
                            pfactor = ...5)
  batch <- combined_data$study_site
  data.harmonized <- covfam(data=structural_data,
                            bat = as.factor(batch),
                            covar = covar_df,
                            gam,
                            y ~ s(age, k=3, fx=F) +
                              as.factor(sex) + euler + pfactor)
}

if (pfactor_include == 'no'){
  covar_df <- bind_cols(combined_data$participant_id,
                        as.numeric(age_vec),
                        as.factor(sex_vec),
                        as.numeric(euler_vec))
  covar_df <- dplyr::rename(covar_df,
                            participant_id=...1,
                            age = ...2,
                            sex = ...3,
                            euler = ...4)
  batch <- combined_data$study_site
  data.harmonized <- covfam(data=structural_data,
                            bat = as.factor(batch),
                            covar = covar_df,
                            gam,
                            y ~ s(age, k=3, fx=F) +
                              as.factor(sex) + euler)
}

# Save the harmonized structural data
data.harmonized_covbat <- data.frame(data.harmonized$dat.covbat)
outputPath <- paste(data_path, sprintf('combined_df_%s_harmonized.tsv', dtype),
                                sep = "")
write.csv(data.harmonized_covbat, outputPath, row.names=FALSE)

#################
# functional data
#################
# # for pfactor
# combined_df_withinbetween_fcrsn7_artifact_pfactor_filter

# # for age
# combined_df_withinbetween_fcrsn7_artifact

# Set parameters for functional harmonization
metric <- 'withinbetween_fcrsn7'
qc_version <- 'artifact'
pfactor_include <- 'no'

# Repeat same workflow for functional data
if (pfactor_include == 'yes'){
  dtype <- sprintf('%s_%s_pfactor_filter', metric, qc_version)
}
if (pfactor_include == 'no'){
  dtype <- sprintf('%s_%s', metric, qc_version)
}

combined_data <- read.csv(paste(data_path, sprintf('combined_df_%s.tsv', dtype),
                                sep = ""),
                          sep = '\t')

# Extract covariates
age_vec <- combined_data$age
sex_vec <- as.factor(combined_data$sex)
fd_vec <- combined_data$medianFD

if (pfactor_include == 'yes'){
  pfactor_vec <- combined_data$p_factor_mcelroy_harmonized_all_samples
}

functional_data = combined_data[0:28] # assumes features are in first 28 columns (within- between- rsn fc)

# Prepare covariate dataframe and harmonize
if (pfactor_include == 'yes'){
  covar_df <- bind_cols(combined_data$participant_id,
                        as.numeric(age_vec),
                        as.factor(sex_vec),
                        as.numeric(fd_vec),
                        as.numeric(pfactor_vec))
  covar_df <- dplyr::rename(covar_df,
                            participant_id=...1,
                            age = ...2,
                            sex = ...3,
                            fd = ...4,
                            pfactor = ...5)
  batch <- combined_data$study_site
  data.harmonized <- covfam(data=functional_data,
                            bat = as.factor(batch),
                            covar = covar_df,
                            gam,
                            y ~ s(age, k=3, fx=F) +
                              as.factor(sex) + fd + pfactor)
}

if (pfactor_include == 'no'){
  covar_df <- bind_cols(combined_data$participant_id,
                        as.numeric(age_vec),
                        as.factor(sex_vec),
                        as.numeric(fd_vec))
  covar_df <- dplyr::rename(covar_df,
                            participant_id=...1,
                            age = ...2,
                            sex = ...3,
                            fd = ...4)
  batch <- combined_data$study_site
  data.harmonized <- covfam(data=functional_data,
                            bat = as.factor(batch),
                            covar = covar_df,
                            gam,
                            y ~ s(age, k=3, fx=F) +
                              as.factor(sex) + fd)
}

# Save the harmonized functional data
data.harmonized_covbat <- data.frame(data.harmonized$dat.covbat)
outputPath <- paste(data_path, sprintf('combined_df_%s_harmonized.tsv', dtype),
                    sep = "")
write.csv(data.harmonized_covbat, outputPath, row.names=FALSE)
