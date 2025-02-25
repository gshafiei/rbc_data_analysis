
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# load data
datapath = '/Users/gshafiei/Desktop/RBC/'

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Myriad Pro']
plt.rcParams['font.size'] = 18.0


def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "***"
    elif pvalue >= 0.0001 and pvalue <= 0.001:
        return "**"
    elif pvalue >= 0.001 and pvalue <= 0.05:
        return "*"
    return ""


####################################
# age effects in intrinsic networks
####################################
metric = 'pfactor'

if metric == 'age':
    dtype = ['noqc',
             'artifact',
             'artifact_harmonized']
    effect_label = 'GAM.age.partialR2'
    pval_label = 'Anova.age.pvaluefdr'
elif metric == 'pfactor':
    dtype = ['noqc_pfactor_filter',
             'artifact_pfactor_filter',
             'artifact_pfactor_filter_harmonized']
    effect_label = 'GAM.pfactor.partialR2'
    pval_label = 'Anova.pfactor.pvaluefdr'

for iType in dtype:
    age_eff = pd.read_csv(datapath + 'results/function/csvFiles/' +
                          'combined_df_withinbetween_fcrsn7' +
                          '_%s_%s_statistics.csv'
                          % (iType, metric))

    # mask the upper triangle
    mask = np.mask_indices(7, np.triu, 0)
    temp_rsn = np.zeros((7, 7))
    temp_rsn[mask] = age_eff[effect_label].values
    temp_rsn = temp_rsn + temp_rsn.T
    # np.fill_diagonal(temp_rsn, 1)

    temp_pvalfdr = np.zeros((7, 7))
    temp_pvalfdr[mask] = age_eff[pval_label].values
    temp_pvalfdr = temp_pvalfdr + temp_pvalfdr.T

    temp_pvalfdr_astr = [convert_pvalue_to_asterisks(pval)
                         for pval in age_eff[pval_label].values]
    temp_pvalfdr_astr = np.array(temp_pvalfdr_astr)
    temp_pvalfdr_astr_mat = np.empty((7, 7), dtype='<U3')
    temp_pvalfdr_astr_mat[mask] = temp_pvalfdr_astr
    i_lower = np.tril_indices(7, -1)
    temp_pvalfdr_astr_mat[i_lower] = temp_pvalfdr_astr_mat.T[i_lower]

    rsn_labels = [netpair.split('.')[1] for netpair
                  in list(age_eff['netpair'])]
    rsn_labels = rsn_labels[:7]
    new_order = ['Default', 'Cont', 'Limbic', 'SalVentAttn', 'DorsAttn',
                 'SomMot', 'Vis']
    reorder_idx = [rsn_labels.index(rsnname) for rsnname in new_order]
    reordered_rsn_ageeff = temp_rsn[np.ix_(reorder_idx, reorder_idx)]
    reordered_rsn_pvalfdr = temp_pvalfdr[np.ix_(reorder_idx, reorder_idx)]
    reordered_rsn_pvalfdr[reordered_rsn_pvalfdr > 0.05] = np.nan
    reordered_rsn_pvalfdr_astr = temp_pvalfdr_astr_mat[np.ix_(reorder_idx,
                                                              reorder_idx)]

    # maxvalue = np.nanmax(np.abs(temp_rsn))
    # print('\nmax value: %s' % str(maxvalue))
    if metric == 'age':
        maxvalue = 0.25
    elif metric == 'pfactor':
        maxvalue = 0.02

    plt.ion()
    ax = sns.heatmap(reordered_rsn_ageeff, cmap='coolwarm',
                     yticklabels=new_order, xticklabels=new_order,
                     annot=reordered_rsn_pvalfdr_astr, fmt="",
                     linewidth=0.5,
                     vmin=-maxvalue, vmax=maxvalue)
    # ax.axes.set_title('Partial R2 - %s' % iType)
    ax.figure.set_figwidth(5)  # 4.5 7
    ax.figure.set_figheight(4)  # 4 6
    plt.tight_layout()
    plt.show()
    plt.savefig(datapath + 'results/function/' +
                'combined_df_withinbetween_fcrsn7_' +
                '%s_%s_statistics_asterisk.svg'
                % (iType, metric),
                bbox_inches='tight', dpi=300,
                transparent=True)
    plt.close()
