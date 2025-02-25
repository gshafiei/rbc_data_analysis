################
# calculate average between and within network connectivity for RBC data
################
import numpy as np
import glob
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components

proj_path = '/cbica/projects/developmental_gradients/'
outpath = proj_path + 'data_postproc/rsn_rbc_artifact/'

labelinfo7 = np.loadtxt(proj_path + 'schaefer/' +
                        'Schaefer2018_400Parcels_7Networks_order_info.txt',
                        dtype='str', delimiter='\t')


# copied from Gael's sklearn.manifold.spectral_embedding_
# https://ibex.readthedocs.io/en/latest/_modules/sklearn/manifold/spectral_embedding_.html
def _graph_connected_component(graph, node_id):
    """Find the largest graph connected components that contains one
    given node

    Parameters
    ----------
    graph : array-like, shape: (n_samples, n_samples)
        adjacency matrix of the graph, non-zero weight means an edge
        between the nodes

    node_id : int
        The index of the query node of the graph

    Returns
    -------
    connected_components_matrix : array-like, shape: (n_samples,)
        An array of bool value indicating the indexes of the nodes
        belonging to the largest connected components of the given query
        node
    """
    n_node = graph.shape[0]
    if sparse.issparse(graph):
        # speed up row-wise access to boolean connection mask
        graph = graph.tocsr()
    connected_nodes = np.zeros(n_node, dtype=bool)
    nodes_to_explore = np.zeros(n_node, dtype=bool)
    nodes_to_explore[node_id] = True
    for _ in range(n_node):
        last_num_component = connected_nodes.sum()
        np.logical_or(connected_nodes, nodes_to_explore, out=connected_nodes)
        if last_num_component >= connected_nodes.sum():
            break
        indices = np.where(nodes_to_explore)[0]
        nodes_to_explore.fill(False)
        for i in indices:
            if sparse.issparse(graph):
                neighbors = graph[i].toarray().ravel()
            else:
                neighbors = graph[i]
            np.logical_or(nodes_to_explore, neighbors, out=nodes_to_explore)
    return connected_nodes


def _graph_is_connected(graph):
    """ Return whether the graph is connected (True) or Not (False)

    Parameters
    ----------
    graph : array-like or sparse matrix, shape: (n_samples, n_samples)
        adjacency matrix of the graph, non-zero weight means an edge
        between the nodes

    Returns
    -------
    is_connected : bool
        True means the graph is fully connected and False means not
    """
    if sparse.isspmatrix(graph):
        # sparse graph, find all the connected components
        n_connected_components, _ = connected_components(graph)
        return n_connected_components == 1
    else:
        # dense graph, find all connected components start from node 0
        return _graph_connected_component(graph, 0).sum() == graph.shape[0]


# ###############
# FC within between rsn: 7 x 7 (rsn x rsn)
# ###############
# use schafer atlas info to assign each network edge a "connected" (1) vs
# "not connected" (0) status
rsnlabels = []
for row in range(0, len(labelinfo7), 2):
    rsnlabels.append(labelinfo7[row].split('_')[2])

edge_status = np.zeros((len(rsnlabels), len(rsnlabels)))
for i in range(len(rsnlabels)):
    for j in range(len(rsnlabels)):
        if rsnlabels[i] == rsnlabels[j]:
            edge_status[i, j] = 1
        else:
            edge_status[i, j] = 0

# network reordering
# CPAC FC data is in Schaefer400-17NetworkOrdering. I'll be reordering
# the nodes to get the networks in Schaefer400-7NetworkOrdering using
# the ordering index below
reorder_dir = proj_path + '/schaefer/schaefer_ordering_mapper/'
mapped = pd.read_csv(reorder_dir + 'Schaefer_400-17_mappedto_400-7.csv')

reorder_idx = mapped['mapped_indices'].values

# loop through datasets and estimate 7x7 mean rsn connectivity
dataset = ['pnc', 'hbn', 'bhrc', 'ccnp', 'nki']
for iDataset in dataset:
    datapath = proj_path + 'data_pmacs/%s_CPAC_artifact/' % iDataset.upper()
    fileNames = sorted(glob.glob(datapath + 'cpac_RBCv0/sub-*/ses-*/func/' +
                                 '*_task-rest*_atlas-Schaefer2018p400n17' +
                                 '_space-MNI152NLin6ASym' +
                                 '_reg-36Parameter' +
                                 '_desc-PearsonNilearn_correlations.tsv'))

    subj_list = [iFile.split('/')[-1].split('_')[0] for iFile in fileNames]
    masked_rsnmean_all = []
    rsnconn_labels_all = []

    print('\n---------------------- ')
    print('\nDataset %s starting:' % iDataset)
    print('\n---------------------- ')

    for i, iFile in enumerate(fileNames):
        # load each data file
        avgFCmri = pd.read_csv(iFile, delimiter='\t', header=None)
        avgFCmri = np.array(avgFCmri)

        # matrix permuted to 7-network ordering
        avgFCmri = avgFCmri[np.ix_(reorder_idx, reorder_idx)]

        if _graph_is_connected(avgFCmri):
            avgFC_copy = avgFCmri.copy()
            np.fill_diagonal(avgFC_copy, np.nan)

            tempdf1 = pd.DataFrame(avgFC_copy)
            tempdf1['rsn'] = rsnlabels

            tempdf1_mean = tempdf1.groupby(['rsn']).mean()

            tempdf2 = tempdf1_mean.T
            tempdf2['rsn'] = rsnlabels

            rsn_means = tempdf2.groupby(['rsn']).mean()
        else:
            print('\nSubj %s FC not fully connected!' % (subj_list[i]))

            # if graph isn't fully connected, replace zeros with NaNs before
            # finding the mean
            zeroidx = np.where(np.sum(avgFCmri, axis=1) == 1)[0]
            avgFCmri[zeroidx, :] = np.nan
            avgFCmri[:, zeroidx] = np.nan

            avgFC_copy = avgFCmri.copy()
            np.fill_diagonal(avgFC_copy, np.nan)

            tempdf1 = pd.DataFrame(avgFC_copy)
            tempdf1['rsn'] = rsnlabels

            tempdf1_mean = tempdf1.groupby(['rsn']).mean()

            tempdf2 = tempdf1_mean.T
            tempdf2['rsn'] = rsnlabels

            rsn_means = tempdf2.groupby(['rsn']).mean()

        # mask the upper triangle of 7x7 rsn means (including the diagonal)
        mask = np.mask_indices(7, np.triu, 0)
        masked_rsnmean = rsn_means.values[mask]

        # get the labels
        rsnmean_labels = list(rsn_means.columns)

        check_nan_mean = rsn_means.isnull().values.any()
        if not check_nan_mean:
            rsnconn_labels = []
            for r, rsnmean in enumerate(masked_rsnmean):
                labelidx = np.where(rsn_means.values == rsnmean)
                rsnconn_labels.append(rsnmean_labels[labelidx[0][0]] + '-' +
                                      rsnmean_labels[labelidx[1][0]])

        masked_rsnmean_all.append(masked_rsnmean)
        rsnconn_labels_all.append(rsnconn_labels)

        print('\nFile %i/%i done!' % (i, len(fileNames)))

    final_labels = list(np.unique(np.array(rsnconn_labels_all)))
    final_df = pd.DataFrame(np.array(masked_rsnmean_all),
                            columns=final_labels)
    final_df['participant_id'] = subj_list
    final_df['full_fileName'] = fileNames

    final_df.to_csv((outpath +
                     'study-%s_task-rest_atlas-schaefer400n17' +
                     '_desc-fcwithinbetween-rsn7.csv') % iDataset,
                    index=False)
