
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

project_path <- '/Users/gshafiei/Desktop/RBC/'
data_path <- paste(project_path, 'data/dataR/', sep = "")

metric <- 'lgi'  # 'gv', 'ct', 'sa', 'lgi'
# we only want to harmonize 'artifact' version
qc_version <- 'artifact'  # 'noqc', 'artifact'
pfactor_include <- 'yes'

if (pfactor_include == 'yes'){
  dtype <- sprintf('%s_%s_pfactor_filter', metric, qc_version)
}
if (pfactor_include == 'no'){
  dtype <- sprintf('%s_%s', metric, qc_version)
}

combined_data <- read.csv(paste(data_path, sprintf('combined_df_%s.tsv', dtype),
                                sep = ""),
                          sep = '\t')

age_vec <- combined_data$age
sex_vec <- as.factor(combined_data$sex)
euler_vec <- combined_data$euler
if (pfactor_include == 'yes'){
pfactor_vec <- combined_data$p_factor_mcelroy_harmonized_all_samples
}

structural_data = combined_data[0:400]

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

# save output
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

project_path <- '/Users/gshafiei/Desktop/RBC/'
data_path <- paste(project_path, 'data/dataR/', sep = "")

metric <- 'withinbetween_fcrsn7'
# we only want to harmonize 'artifact' version
qc_version <- 'artifact'  # 'noqc', 'artifact'
pfactor_include <- 'no'

if (pfactor_include == 'yes'){
  dtype <- sprintf('%s_%s_pfactor_filter', metric, qc_version)
}
if (pfactor_include == 'no'){
  dtype <- sprintf('%s_%s', metric, qc_version)
}

combined_data <- read.csv(paste(data_path, sprintf('combined_df_%s.tsv', dtype),
                                sep = ""),
                          sep = '\t')

age_vec <- combined_data$age
sex_vec <- as.factor(combined_data$sex)
fd_vec <- combined_data$medianFD

if (pfactor_include == 'yes'){
  pfactor_vec <- combined_data$p_factor_mcelroy_harmonized_all_samples
}

functional_data = combined_data[0:28]

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

# save output
data.harmonized_covbat <- data.frame(data.harmonized$dat.covbat)
outputPath <- paste(data_path, sprintf('combined_df_%s_harmonized.tsv', dtype),
                    sep = "")
write.csv(data.harmonized_covbat, outputPath, row.names=FALSE)


###### not sure what this is -- check ######
# # combined_data <- read.csv('/Users/gshafiei/Desktop/TED_CBF/data/desc-standardForR_cbf_nki_age_2070.tsv', sep = '\t')
# 
# age_vec <- combined_data$age
# sex_vec <- as.factor(combined_data$sex)
# # euler_vec <- combined_data$euler
# euler_vec <- combined_data$meanFD
# pfactor_vec <- combined_data$p_factor_mcelroy_harmonized_all_samples
# 
# # ct_data = combined_data[4:5]
# # ct_data = combined_data[0:400]
# ct_data = combined_data[0:28]
# # ct_data = combined_data[0:2]

# covar_df <- bind_cols(as.numeric(age_vec), as.factor(sex_vec))
# covar_df <- dplyr::rename(covar_df,
#                           age = ...1,
#                           sex = ...2)
# batch <- combined_data$dataset_2
# 
# data.harmonized <- covfam(data=ct_data, bat = as.factor(batch), covar = covar_df, gam, y ~ s(age, k=3, fx=F) + as.factor(sex))
# data.harmonized_covbat <- data.frame(data.harmonized$dat.covbat)
# outputPath <- '/Users/gshafiei/Desktop/TED_CBF/data/desc-standardForR_cbf_nki_age_2070_harmonized.csv'
# write.csv(data.harmonized_covbat, outputPath, row.names=FALSE)