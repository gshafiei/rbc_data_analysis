################
# prepare data for plotting and GAM in R
################
import numpy as np
import glob
import pandas as pd

datapath = '/Users/gshafiei/Desktop/RBC/'

qc_version = 'artifact'  # 'artifact' or 'noqc'
################
# Generate tsv files for RBC structural data
################
dtype = ['gv', 'ct', 'sa', 'lgi']
dataset = ['bhrc', 'ccnp', 'hbn', 'nki', 'pnc']

for iDataset in dataset:
    if qc_version == 'noqc':
        fileNames = sorted(glob.glob(datapath + 'data/%s_FreeSurfer/'
                                     '*_regionsurfacestats.tsv'
                                     % iDataset.upper()))
    elif qc_version == 'artifact':
        fileNames = sorted(glob.glob(datapath + 'data/%s_FreeSurfer_artifact/'
                                     '*_regionsurfacestats.tsv'
                                     % iDataset.upper()))

    subj_list = [iFile.split('/')[-1].split('_')[0].split('-')[1]
                 for iFile in fileNames]

    gv = []
    ct = []
    sa = []
    lgi = []

    for i, iFile in enumerate(fileNames):
        temp_file = iFile
        brain_data = pd.read_csv(iFile, delimiter='\t')

        schaefer400 = brain_data.loc[brain_data['atlas'] ==
                                     'Schaefer2018_400Parcels_7Networks_order',
                                     ('StructName', 'GrayVol', 'ThickAvg',
                                      'SurfArea', 'Mean_piallgi')]
        schaefer400.reset_index(inplace=True, drop=True)
        medial_wall_idx = np.where(np.array(schaefer400['StructName']) ==
                                   'Background+FreeSurfer_Defined_Medial_Wall')
        medial_wall_idx = medial_wall_idx[0]
        schaefer400.drop(medial_wall_idx, inplace=True)
        schaefer400.reset_index(inplace=True, drop=True)

        region_names = list(schaefer400['StructName'])
        gv.append(np.array(schaefer400['GrayVol']))
        ct.append(np.array(schaefer400['ThickAvg']))
        sa.append(np.array(schaefer400['SurfArea']))
        lgi.append(np.array(schaefer400['Mean_piallgi']))
        print('\nFile %i/%i done!' % (i, len(fileNames)))

    for iType in dtype:
        demogs = pd.read_csv(datapath + 'data/demogs/' +
                             'study-%s_desc-participants.tsv'
                             % iDataset.upper(),
                             delimiter='\t')
        demogs['participant_id'] = np.array(demogs['participant_id']).astype(
            'str')

        qc_data = pd.read_csv(datapath + 'data/QC_data/' +
                              'study-%s_desc-T1_qc.tsv'
                              % iDataset.upper(),
                              delimiter='\t')
        qc_data['participant_id'] = np.array(qc_data['participant_id']).astype(
            'str')

        # generate metric dataframe
        if iType == 'ct':
            df_metric = pd.DataFrame(np.array(ct), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(ct), axis=1)
        elif iType == 'gv':
            df_metric = pd.DataFrame(np.array(gv), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(gv), axis=1)
        elif iType == 'sa':
            df_metric = pd.DataFrame(np.array(sa), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(sa), axis=1)
        elif iType == 'lgi':
            df_metric = pd.DataFrame(np.array(lgi), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(lgi), axis=1)

        if iDataset == 'bhrc':
            demogs_ses1 = demogs.query('session_id == 1')
            demogs_ses1.reset_index(inplace=True, drop=True)
            del demogs
            demogs = demogs_ses1.copy()
            del demogs_ses1

            qc_data_ses1 = qc_data.query('session_id == 1')
            qc_data_ses1.reset_index(inplace=True, drop=True)
            del qc_data
            qc_data = qc_data_ses1.copy()
            del qc_data_ses1

        elif iDataset == 'nki':
            demogs_ses1 = demogs[demogs['session_id'].str.contains('BAS1')]
            demogs_ses1.reset_index(inplace=True, drop=True)
            del demogs
            demogs = demogs_ses1.copy()
            del demogs_ses1

            qc_data_ses1 = qc_data[qc_data['session_id'].str.contains('BAS1')]
            qc_data_ses1.reset_index(inplace=True, drop=True)
            del qc_data
            qc_data = qc_data_ses1.copy()
            del qc_data_ses1

        # add demographics
        x = np.array(df_metric['participant_id'])
        y = np.array(demogs['participant_id'])
        xy, x_ind, y_ind = np.intersect1d(x, y, return_indices=True)

        df_demogs = df_metric.iloc[x_ind, :]
        df_demogs.reset_index(inplace=True, drop=True)

        demogs_shared = demogs.iloc[y_ind, :]
        demogs_shared.reset_index(inplace=True, drop=True)

        df_demogs_shared = pd.merge(df_demogs, demogs_shared,
                                    on='participant_id', how='left')

        # add T1 qc info
        x = np.array(df_demogs_shared['participant_id'])
        y = np.array(qc_data['participant_id'])
        xy, x_ind, y_ind = np.intersect1d(x, y, return_indices=True)

        df_qc = df_demogs_shared.iloc[x_ind, :]
        df_qc.reset_index(inplace=True, drop=True)

        qc_data_shared = qc_data.iloc[y_ind, :]
        qc_data_shared.reset_index(inplace=True, drop=True)

        df_qc['euler'] = qc_data_shared['euler'].values
        df_qc['qc_determination'] = qc_data_shared['qc_determination'].values

        # # check whether everyone has euler and age values
        # check_nan_euler = df_qc['euler'].isnull().values.any()
        # check_nan_age = df_qc['age'].isnull().values.any()
        # if check_nan_euler or check_nan_age:
        #     euler_nanidx = np.where(np.isnan(df_qc['euler'].values))
        #     age_nanidx = np.where(np.isnan(df_qc['age'].values))
        #     all_nanidx = list(age_nanidx[0]) + list(euler_nanidx[0])
        #     df_qc.drop(np.array(all_nanidx), inplace=True)
        #     df_qc.reset_index(inplace=True, drop=True)

        df_final = df_qc.copy()
        df_final = df_final.query('6 <= age <= 22')
        df_final.reset_index(inplace=True, drop=True)

        if qc_version == 'noqc':
            output_filename = (datapath + 'data/dataR/%s_df_%s_noqc.tsv'
                               % (iDataset, iType))
        elif qc_version == 'artifact':
            output_filename = (datapath + 'data/dataR/%s_df_%s_artifact.tsv'
                               % (iDataset, iType))

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
dtype = ['ct_noqc', 'gv_noqc', 'sa_noqc', 'lgi_noqc',
         'ct_artifact', 'gv_artifact', 'sa_artifact', 'lgi_artifact',
         'ct_noqc_pfactor_filter', 'gv_noqc_pfactor_filter',
         'sa_noqc_pfactor_filter', 'lgi_noqc_pfactor_filter',
         'ct_artifact_pfactor_filter', 'gv_artifact_pfactor_filter',
         'sa_artifact_pfactor_filter', 'lgi_artifact_pfactor_filter']

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

    if 'lgi' in iType:
        iType_data = np.array(combined_df.iloc[:, 0:400])
        subj_nanidx = np.where(np.isnan(iType_data).all(axis=1))[0]
        # subj_nanidx = np.unique(np.where(np.isnan(iType_data))[0])
        combined_df.drop(subj_nanidx, inplace=True)
        combined_df.reset_index(inplace=True, drop=True)

    combined_df.to_csv(datapath + 'data/dataR/combined_df_%s.tsv' % iType,
                       index=False, sep='\t')

################
# regenerate tsv files for harmonized data
################
# load harmonized data and regenerate tsv
dtype = ['ct_artifact', 'gv_artifact', 'sa_artifact', 'lgi_artifact',
         'ct_artifact_pfactor_filter', 'gv_artifact_pfactor_filter',
         'sa_artifact_pfactor_filter', 'lgi_artifact_pfactor_filter']
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

    region_names = list(combined_df.columns[0:400])
    combined_df_harmonized = pd.DataFrame(np.array(combined_metric_harmonized),
                                          columns=region_names)
    harmonized_meanVal = np.mean(np.array(combined_metric_harmonized), axis=1)
    final_df = pd.concat([combined_df_harmonized,
                          combined_df.iloc[:, 400:]], axis=1)
    final_df['meanValHarmonized'] = harmonized_meanVal
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
