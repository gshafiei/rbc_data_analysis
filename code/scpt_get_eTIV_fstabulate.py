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
dataset = ['bhrc', 'ccnp', 'hbn', 'nki', 'pnc']

for iDataset in dataset:
    if qc_version == 'noqc':
        fileNames = sorted(glob.glob(datapath + 'data/%s_FreeSurfer/'
                                     '*_brainmeasures.tsv'
                                     % iDataset.upper()))
    elif qc_version == 'artifact':
        fileNames = sorted(glob.glob(datapath + 'data/%s_FreeSurfer_artifact/'
                                     '*_brainmeasures.tsv'
                                     % iDataset.upper()))

    subj_list = [iFile.split('/')[-1].split('_')[0].split('-')[1]
                 for iFile in fileNames]

    # extract EstimatedTotalIntraCranialVol_eTIV from freesurfer data
    # gray matter volume
    eTIV = []

    for i, iFile in enumerate(fileNames):
        temp_file = iFile
        brain_data = pd.read_csv(iFile, delimiter='\t')

        eTIV.append(brain_data['EstimatedTotalIntraCranialVol_eTIV'].values)
        print('\nFile %i/%i done!' % (i, len(fileNames)))

    df_eTIV = pd.DataFrame(np.array(eTIV), columns=['eTIV'])
    df_eTIV['participant_id'] = subj_list

    df_eTIV.to_csv(datapath + 'data/eTIV/%s_df_eTIV_artifact.tsv' % iDataset,
                   index=False, sep='\t')

################
# combine all 5 datasets
################
dtype = ['eTIV_artifact']

for iType in dtype:
    pnc = pd.read_csv(datapath + 'data/eTIV/pnc_df_%s.tsv' % iType,
                      delimiter='\t')
    hbn = pd.read_csv(datapath + 'data/eTIV/hbn_df_%s.tsv' % iType,
                      delimiter='\t')
    bhrc = pd.read_csv(datapath + 'data/eTIV/bhrc_df_%s.tsv' % iType,
                       delimiter='\t')
    nki = pd.read_csv(datapath + 'data/eTIV/nki_df_%s.tsv' % iType,
                      delimiter='\t')
    ccnp = pd.read_csv(datapath + 'data/eTIV/ccnp_df_%s.tsv' % iType,
                       delimiter='\t')

    combined_df = pd.concat([pnc, hbn, bhrc, nki, ccnp], axis=0)

    combined_df.reset_index(inplace=True, drop=True)

    combined_df.to_csv(datapath + 'data/eTIV/combined_df_%s.tsv' % iType,
                       index=False, sep='\t')

################
# load combined gv and sa R data and eTIV data to them
################
dtype = ['artifact_harmonized']

for iType in dtype:
    sa_data = pd.read_csv(datapath + 'data/dataR/combined_df_sa_%s.tsv'
                          % iType,
                          delimiter='\t',
                          low_memory=False)
    gv_data = pd.read_csv(datapath + 'data/dataR/combined_df_gv_%s.tsv'
                          % iType,
                          delimiter='\t',
                          low_memory=False)
eTIV_data = pd.read_csv(datapath + 'data/eTIV/combined_df_eTIV_artifact.tsv',
                        delimiter='\t',
                        low_memory=False)

# add eTIV
x = np.array(sa_data['participant_id'])
y = np.array(eTIV_data['participant_id'])
xy, x_ind, y_ind = np.intersect1d(x, y, return_indices=True)

df_sa_shared = sa_data.iloc[x_ind, :]
df_sa_shared.reset_index(inplace=True, drop=True)

eTIV_shared = eTIV_data.iloc[y_ind, :]
eTIV_shared.reset_index(inplace=True, drop=True)

df_sa_eTIV = pd.merge(df_sa_shared, eTIV_shared,
                      on='participant_id', how='left')

df_sa_eTIV.to_csv(datapath + 'data/dataR/combined_df_sa_%s_eTIV.tsv'
                  % iType, index=False, sep='\t')

# add eTIV
x = np.array(gv_data['participant_id'])
y = np.array(eTIV_data['participant_id'])
xy, x_ind, y_ind = np.intersect1d(x, y, return_indices=True)

df_gv_shared = gv_data.iloc[x_ind, :]
df_gv_shared.reset_index(inplace=True, drop=True)

eTIV_shared = eTIV_data.iloc[y_ind, :]
eTIV_shared.reset_index(inplace=True, drop=True)

df_gv_eTIV = pd.merge(df_gv_shared, eTIV_shared,
                      on='participant_id', how='left')

df_gv_eTIV.to_csv(datapath + 'data/dataR/combined_df_gv_%s_eTIV.tsv'
                  % iType, index=False, sep='\t')
