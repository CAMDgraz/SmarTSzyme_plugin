"""
@authors: Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sanchez-Murcia [pedro.murcia@medunigraz.at]
"""

import glob
import sys
import numpy as np
import pandas as pd
import heapq
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pickle


def find_matrices(path: str) -> np.ndarray:
    """
    Return a list with the matrices id in the path

    Parameters
    ----------
        path: str
            Path to the interaction matrices
    
    Return
    ------
        matrices_id: np .ndarray
            Numpy array with the matrices id in path
    """
    es_matrices = glob.glob(f'{path}/es_*.pickle')
    pts_matrices = glob.glob(f'{path}/pts_*.pickle')

    if len(es_matrices) != len(pts_matrices):
        print('\nError!!! Number of matrices for the es different than for the pts')
        sys.exit()
    
    matrices_id = np.unique([int(matrix.split('.')[-2]) for matrix in
                             es_matrices])
    return matrices_id

def load_pickle(file):
    with open(file, "rb") as f:
        data = pickle.load(f)
    return data

def score(coupling_path, matrices_path, nres, batch):

    # Check interaction matrices and coupling files ================================
    # Get number of couplings
    coupling_stab = glob.glob(f'{coupling_path}/stab_flux.*.csv')
    coupling_dest = glob.glob(f'{coupling_path}/dest_flux.*.csv')

    # Get interaction matrices 
    matrices_id = find_matrices(matrices_path)

    if len(coupling_dest) != len(coupling_stab):
        print('\nError!!! The number of stabilizing and destabilizing files differ')
        return
    elif len(coupling_stab) != len(matrices_id):
        print('\nError!!! The number of interaction matrices and couplings differ')
        return

    batches = np.arange(batch, len(coupling_stab) + batch, batch)
    print(f'No. of couplings file:          {len(coupling_stab)}')
    print(f'No. Batches:                    {len(batches)}')

    # ==========================================================================

    # Main analysis ============================================================
    score_batches = np.zeros((len(batches), nres))

    print(f'\n*** Analysing results ***')
    for batch_id, batch in enumerate(batches):
        # Variables for stab and destab contributions
        stab_dest = np.zeros((batch, nres))
        flux_stab_batch = np.zeros((batch, nres))
        flux_dest_batch = np.zeros((batch, nres))

        for job_id, job in enumerate(matrices_id[:batch]):
            results_stab = pd.read_csv(f'{coupling_path}/stab_flux.{job}.csv',
                                    sep=',')
            results_dest = pd.read_csv(f'{coupling_path}/dest_flux.{job}.csv',
                                    sep=',')
            flux_stab_batch[job] = np.asarray(results_stab.iloc[:, 1])
            flux_dest_batch[job] = np.asarray(results_dest.iloc[:, 1])

            # get es and pts interactions
            es_int = np.zeros((nres, nres))
            pts_int = np.zeros((nres, nres))

            es_files = glob.glob(f'{matrices_path}/es_*.{job}.pickle')
            pts_files = glob.glob(f'{matrices_path}/pts_*.{job}.pickle')

            for es_file in es_files:
                es_int += load_pickle(es_file)
            for pts_file in pts_files:
                pts_int += load_pickle(pts_file)
            

            diff_matrix = pts_int - es_int
            stab_dest[job_id] = diff_matrix.sum(axis=0)/abs(diff_matrix.sum(axis=0))

        dscore = (stab_dest > 0).sum(axis=0)/batch
        sscore = (stab_dest < 0).sum(axis=0)/batch
        dest_score = (dscore*flux_dest_batch).sum(axis=0) 
        stab_score = (sscore*flux_stab_batch).sum(axis=0)
        score_batches[batch_id] = dest_score - stab_score
        print(f'Batch {batch_id + 1}: Done')
    
    return score_batches, batches

def get_top(score_batches, topn, nres):
        # Get topn
    print(f'\n*** Getting top {topn} residues ***')
    top_stab_residues = np.zeros(score_batches.shape)
    top_dest_residues = np.zeros(score_batches.shape)

    for batch_id, score in enumerate(score_batches):
        ordered = []
        heapq.heapify(ordered)
        for res in range(nres):
            heapq.heappush(ordered, (score[res], res))
        
        nlargest = [i[1] for i in heapq.nlargest(topn, ordered)]
        nsmallest = [i[1] for i in heapq.nsmallest(topn, ordered)]

        top_dest_residues[batch_id, nlargest] = 1
        top_stab_residues[batch_id, nsmallest] = 1
    return top_stab_residues, top_dest_residues

def plot_score(score_batches, batches, topn, nres, ax, plot_type):
    top_stab_residues, top_dest_residues = get_top(score_batches, topn, nres)

    print(f'Maximum mutatibility score:')
    for batch_id, res_batch in enumerate(top_stab_residues):
        print(f'Batch {batch_id}: {np.where(res_batch == 1)[0] + 1}')

    print(f'Minimum mutatibility score:')
    for batch_id, res_batch in enumerate(top_dest_residues):
        print(f'Batch {batch_id}: {np.where(res_batch == 1)[0] + 1}')
    
    top_stab = np.unique(np.where(top_stab_residues == 1)[1])
    top_dest = np.unique(np.where(top_dest_residues == 1)[1])

    # Maximum mutability
    if plot_type == "Top Max":
        sns.heatmap(top_dest_residues[:, top_dest], cmap='cubehelix_r', vmin=0, vmax=1,
                    center=0.5, cbar=False, linewidths=0.5, linecolor='grey', ax=ax,
                    square=True)
        ax.set_yticks(np.arange(0.5, len(batches) + 0.5, 1))
        ax.set_yticklabels(batches, rotation=0, size=10)
        ax.set_xticks(np.arange(0.5, len(top_dest) + 0.5, 1))
        ax.set_xticklabels(top_dest + 1, rotation=45, size=10)

        ax.set_ylabel(r'Number of Trajectories')
        ax.set_xlabel(r'Residue id')

    # Minimum mutability
    if plot_type == "Top Min":
        sns.heatmap(top_stab_residues[:, top_stab], cmap='cubehelix_r', vmin=0, vmax=1,
                    center=0.5, cbar=False, linewidths=0.5, linecolor='grey', ax=ax,
                    square=True)
        ax.set_yticks(np.arange(0.5, len(batches) + 0.5, 1))
        ax.set_yticklabels(batches, rotation=0, size=10)
        ax.set_xticks(np.arange(0.5, len(top_stab) + 0.5, 1))
        ax.set_xticklabels(top_stab + 1, rotation=45, size=10)

        ax.set_ylabel(r'Number of Trajectories')
        ax.set_xlabel(r'Residue id')

def plot_experimental(expres, score_batches, batches, ax):
      
    exp_res = np.asarray(expres, dtype=int) - 1

    max_val = score_batches[:, exp_res].max()
    min_val = score_batches[:, exp_res].min()

    vmax = max_val if max_val > abs(min_val) else abs(min_val)

    sns.heatmap(score_batches[:, exp_res], cmap='vlag', center=0,
                ax=ax, vmin=-vmax, vmax=vmax, square=True)
    ax.set_yticks(np.arange(0.5, len(batches) + 0.5, 1))
    ax.set_yticklabels(batches, rotation=0, size=10)
    ax.set_xticks(np.arange(0.5, len(exp_res) + 0.5, 1))
    ax.set_xticklabels(exp_res + 1, rotation=45, size=10)

    ax.set_ylabel(r'Number of Trajectories')
    ax.set_xlabel(r'Residue id')
    return