################
# prepare data for plotting and GAM in R
################
import numpy as np
import pandas as pd

datapath = '/Users/gshafiei/Desktop/RBC/'

qc_version = 'noqc'  # 'artifact' or 'noqc'
################
# Generate tsv files for R for RBC functional data data
################
dataset = ['bhrc', 'ccnp', 'hbn', 'nki', 'pnc']

for iDataset in dataset:
    if qc_version == 'noqc':
        fcdata = pd.read_csv((datapath + 'data/rsn_rbc/' +
                              'study-%s_task-rest_atlas-schaefer400n17' +
                              '_desc-fcwithinbetween-rsn7.csv')
                             % iDataset, low_memory=False)
    elif qc_version == 'artifact':
        fcdata = pd.read_csv((datapath + 'data/rsn_rbc_artifact/' +
                              'study-%s_task-rest_atlas-schaefer400n17' +
                              '_desc-fcwithinbetween-rsn7.csv')
                             % iDataset, low_memory=False)

    subj_list = [iSubj.split('-')[1]
                 for iSubj in list(fcdata['participant_id'])]
    fileNames = list(fcdata['full_fileName'])

    # use the longest file name to create an empty dataframe with all columns
    check_filename = [len(fileNames[iSubj].split('/')[-1].split('_'))
                      for iSubj in range(len(subj_list))]
    unique_namelength = np.unique(np.array(check_filename))
    maxidx = np.where(np.array(check_filename)
                      == unique_namelength.max())[0][0]

    # generate main df
    split_name = fileNames[maxidx].split('/')[-1].split('_')
    col_names_max = [split_title.split('-')[0]
                     for split_title in split_name[:-3]]
    df_main = pd.DataFrame(columns=col_names_max)

    # separate the between within connectivity from the rest
    fc_df = fcdata.iloc[:, :28]

    for iSubj in range(len(subj_list)):
        # use filenames to get the column info
        split_name = fileNames[iSubj].split('/')[-1].split('_')
        col_names = [split_title.split('-')[0]
                     for split_title in split_name[:-3]]
        df_temp = pd.DataFrame(columns=col_names)
        col_vals = [split_title.split('-')[1]
                    for split_title in split_name[:-3]]
        df_temp.loc[0] = col_vals

        df_main = pd.concat([df_main, df_temp], ignore_index=True)

    # combine network connectivity data with other info
    fc_df = pd.concat([fc_df, df_main], axis=1)

    # load demographocs
    demogs = pd.read_csv(datapath + 'data/demogs/'
                         'study-%s_desc-participants.tsv' % iDataset,
                         delimiter='\t')
    demogs['participant_id'] = demogs['participant_id'].astype(str)

    # load fc qc to get medianFD
    qc_data_fc = pd.read_csv(datapath + 'data/QC_data/' +
                             'study-%s_desc-functional_qc.tsv'
                             % iDataset.upper(),
                             delimiter='\t')
    qc_data_fc['participant_id'] = qc_data_fc['participant_id'].astype(str)

    if qc_version == 'artifact':
        # only keep individuals who've passed FC QC
        qc_data_fc = qc_data_fc.query('fmriExclude == 0')
        qc_data_fc.reset_index(inplace=True, drop=True)

    # get specific task data (here: task is "rest")
    qc_data_task = qc_data_fc[qc_data_fc['task'].str.contains('rest')]
    qc_data_task.reset_index(inplace=True, drop=True)

    # keep first session data if data is BHRC or NKI
    if iDataset == 'bhrc':
        # first session data: demogs
        demogs_ses1 = demogs.query('session_id == 1')
        demogs_ses1.reset_index(inplace=True, drop=True)
        del demogs
        demogs = demogs_ses1.copy()
        del demogs_ses1

        # first session data: fc qc
        qc_data_fc_ses1 = qc_data_task.query('session_id == 1')
        qc_data_fc_ses1.reset_index(inplace=True, drop=True)
        del qc_data_task
        qc_data_task = qc_data_fc_ses1.copy()
        del qc_data_fc_ses1

        # first session data: mean fc rsn
        temp_df_ses1 = fc_df[fc_df['ses'].str.contains('1')]
        temp_df_ses1.reset_index(inplace=True, drop=True)
        df_metric = temp_df_ses1.copy()
        del temp_df_ses1
    elif iDataset == 'nki':
        # first session data: demogs
        demogs_ses1 = demogs[demogs['session_id'].str.contains('BAS1')]
        demogs_ses1.reset_index(inplace=True, drop=True)
        del demogs
        demogs = demogs_ses1.copy()
        del demogs_ses1

        # first session data: fc qc
        qc_data_fc_ses1 = qc_data_task[
            qc_data_task['session_id'].str.contains('BAS1')]
        qc_data_fc_ses1.reset_index(inplace=True, drop=True)
        del qc_data_task
        qc_data_task = qc_data_fc_ses1.copy()
        del qc_data_fc_ses1

        # first session data: mean fc rsn
        temp_df_ses1 = fc_df[fc_df['ses'].str.contains('BAS1')]
        temp_df_ses1.reset_index(inplace=True, drop=True)
        df_metric = temp_df_ses1.copy()
        del temp_df_ses1
    else:
        df_metric = fc_df.copy()

    # group by particpant for mean FD and df metric (within above DFs)
    if iDataset != 'hbn':
        assert df_metric['ses'].unique().shape[0] == 1, 'multiple ses labels'
    assert df_metric['task'].unique().shape[0] == 1, 'multiple task labels'
    assert df_metric['atlas'].unique().shape[0] == 1, 'multiple atlas labels'
    assert df_metric['space'].unique().shape[0] == 1, 'multiple space labels'

    # now that we are sure 'ses', 'task', 'atlas', 'space' are the same across
    # the df_metric dataframe, we drop them.
    # also drop 'acq', we'll average across runs
    df_metric.drop(['ses', 'task', 'atlas', 'space', 'acq'],
                   axis=1, inplace=True)
    if iDataset in ['bhrc', 'hbn', 'ccnp']:
        df_metric.drop(['run'], axis=1, inplace=True)
    df_metric = df_metric.rename(columns={'sub': 'participant_id'})
    df_metric_mean = df_metric.groupby(['participant_id']).mean()

    qc_data_fc_subset = qc_data_task[['participant_id',
                                      'medianFD', 'meanFD']].copy()
    qc_data_mean = qc_data_fc_subset.groupby(['participant_id']).mean()

    df_metric_mergeQC = pd.merge(df_metric_mean, qc_data_mean,
                                 on=['participant_id'])
    df_metric_QCdemog = pd.merge(df_metric_mergeQC, demogs,
                                 on='participant_id')

    column_to_move = df_metric_QCdemog.pop('participant_id')
    df_metric_QCdemog.insert(loc=28, column='participant_id',
                             value=column_to_move)

    # check whether everyone has FD and age values
    check_nan_medianFD = df_metric_QCdemog['medianFD'].isnull().values.any()
    check_nan_meanFD = df_metric_QCdemog['meanFD'].isnull().values.any()
    check_nan_age = df_metric_QCdemog['age'].isnull().values.any()
    assert bool(check_nan_medianFD) == 0, 'median FD has NaNs!'
    assert bool(check_nan_meanFD) == 0, 'mean FD has NaNs!'
    assert bool(check_nan_age) == 0, 'age has NaNs!'

    df_final = df_metric_QCdemog.copy()
    df_final = df_final.query('6 <= age <= 22')
    df_final.reset_index(inplace=True, drop=True)

    if qc_version == 'noqc':
        output_filename = (datapath + 'data/dataR/' +
                           '%s_df_withinbetween_fcrsn7_noqc.tsv'
                           % iDataset)
    elif qc_version == 'artifact':
        output_filename = (datapath + 'data/dataR/' +
                           '%s_df_withinbetween_fcrsn7_artifact.tsv'
                           % iDataset)

    df_final.to_csv(output_filename, index=False, sep='\t')

    # filter p-factor and save
    p_nanidx = np.where(np.isnan(df_final['p_factor_mcelroy' +
                                          '_harmonized' +
                                          '_all_samples'].values))
    df_final.drop(np.array(p_nanidx[0]), inplace=True)
    df_final.reset_index(inplace=True, drop=True)
    output_filename = (output_filename.split('.')[0] +
                       '_pfactor_filter.tsv')

    df_final.to_csv(output_filename, index=False, sep='\t')

################
# combine all 5 datasets
################
dtype = ['withinbetween_fcrsn7_noqc',
         'withinbetween_fcrsn7_artifact',
         'withinbetween_fcrsn7_noqc_pfactor_filter',
         'withinbetween_fcrsn7_artifact_pfactor_filter']

for iType in dtype:
    pnc = pd.read_csv(datapath + 'data/dataR/pnc_df_%s.tsv' % iType,
                      delimiter='\t')
    hbn = pd.read_csv(datapath + 'data/dataR/hbn_df_%s.tsv' % iType,
                      delimiter='\t')
    bhrc = pd.read_csv(datapath + 'data/dataR/bhrc_df_%s.tsv' % iType,
                       delimiter='\t')
    nki = pd.read_csv(datapath + 'data/dataR/nki_df_%s.tsv' % iType,
                      delimiter='\t')
    ccnp = pd.read_csv(datapath + 'data/dataR/ccnp_df_%s.tsv' % iType,
                       delimiter='\t')

    combined_df = pd.concat([pnc, hbn, bhrc, nki, ccnp], axis=0)

    combined_df.reset_index(inplace=True, drop=True)

    combined_df.to_csv(datapath + 'data/dataR/combined_df_%s.tsv' % iType,
                       index=False, sep='\t')

################
# regenerate tsv files for harmonized data
################
# load harmonized data and regenerate tsv
dtype = ['withinbetween_fcrsn7_artifact',
         'withinbetween_fcrsn7_artifact_pfactor_filter']
dataset = ['bhrc', 'ccnp', 'hbn', 'nki', 'pnc']

for iType in dtype:
    combined_df = pd.read_csv(datapath + 'data/dataR/combined_df_%s.tsv'
                              % iType,
                              delimiter='\t',
                              low_memory=False)
    combined_metric_harmonized = pd.read_csv(datapath + 'data/dataR/' +
                                             'combined_df_%s_harmonized.tsv'
                                             % iType,
                                             low_memory=False)

    labels = list(combined_df.columns[0:28])
    combined_df_harmonized = pd.DataFrame(np.array(combined_metric_harmonized),
                                          columns=labels)
    final_df = pd.concat([combined_df_harmonized,
                          combined_df.iloc[:, 28:]], axis=1)
    final_df.to_csv(datapath + 'data/dataR/combined_df_%s_harmonized.tsv'
                    % iType,
                    index=False, sep='\t')

    for iDataset in dataset:
        tmp_name = iDataset.upper()
        if iDataset == 'ccnp':
            tmp_name = 'Colornest'
        df_dataset = final_df[final_df['study'].str.contains(tmp_name)]
        df_dataset.reset_index(inplace=True, drop=True)
        df_dataset.to_csv(datapath + 'data/dataR/%s_df_%s_harmonized.tsv'
                          % (iDataset, iType),
                          index=False, sep='\t')
