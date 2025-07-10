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


def calculate_coupling(path_matrices, nres, interactions, catalytic_residues,
                       coupling_out):
    print('\n*** Loading interaction matrices and calculating coupling ***')
    matrices_id = find_matrices(path_matrices)

    # Calculate coupling for each reaction mechanism
    for job in matrices_id:
        # es matrix
        matrix_es = np.zeros((nres, nres))
        int_es = []
        for interaction in interactions:
            int_es.extend(glob.glob(f'{path_matrices}/es_{interaction}.{job}.pickle'))
        for matrix in int_es:
            matrix_es += load_pickle(matrix)

        # pts matrix
        matrix_pts = np.zeros((nres, nres))
        int_pts = []
        for interaction in interactions:
            int_pts.extend(glob.glob(f'{path_matrices}/pts_{interaction}.{job}.pickle'))
        for matrix in int_pts:
            matrix_pts += load_pickle(matrix)
        
        matrix_diff = matrix_pts - matrix_es

        del matrix_pts
        del matrix_es

        # calculate edge-to-edge
        edge_to_edge, edges = edge_transfer_matrix(matrix_diff)
        dweights = get_weights(matrix_diff, edges)

        flux_stab = np.zeros(nres)
        flux_dest = np.zeros(nres)

        # get edges involved in the catalytic residues
        edges_catal = []
        for cr in catalytic_residues:
            for edge_id, vertices in enumerate(edges):
                if cr == vertices[0] or cr == vertices[1]:
                    edges_catal.append(edge_id)

        # do not consider the diagonal elements (or embedness)
        for diag in range(edge_to_edge.shape[0]):
            edge_to_edge[diag][diag] = 0
        
        for residue in range(nres):
            if residue in catalytic_residues:
                continue
            edges_residue_stab = []
            edges_residue_dest = []
            for edge_id, vertices in enumerate(edges):
                if residue == vertices[0] or residue == vertices[1]:
                    if dweights[edge_id][edge_id] < 0:
                        edges_residue_stab.append(edge_id)
                    elif dweights[edge_id][edge_id] > 0:
                        edges_residue_dest.append(edge_id)
        
            # stab fluxes
            edges_catal_tile = np.tile(edges_catal, len(edges_residue_stab))
            edges_catal_tile = edges_catal_tile.astype(np.int16)
            edges_residue_repeat = np.repeat(edges_residue_stab, len(edges_catal))
            edges_residue_repeat = edges_residue_repeat.astype(np.int16)
            flux_stab[residue] = np.sum(abs(edge_to_edge[edges_residue_repeat,
                                                        edges_catal_tile]))

            # dest fluxes
            edges_catal_tile = np.tile(edges_catal, len(edges_residue_dest))
            edges_catal_tile = edges_catal_tile.astype(np.int16)
            edges_residue_repeat = np.repeat(edges_residue_dest, len(edges_catal))
            edges_residue_repeat = edges_residue_repeat.astype(np.int16)
            flux_dest[residue] = np.sum(abs(edge_to_edge[edges_residue_repeat,
                                                        edges_catal_tile]))
        # normalize fluxes
        total_flux = abs(edge_to_edge[:, edges_catal]).sum()

        flux_stab /= total_flux
        flux_dest /= total_flux

        # Print stab and destab ====================================================
        with open(f'{coupling_out}/coupling/stab_flux.{job}.csv', 'w') as f:
            f.write('residue,flux\n')
            for res, flux in enumerate(flux_stab):
                f.write(f'{res},{flux}\n')

        with open(f'{coupling_out}/coupling/dest_flux.{job}.csv', 'w') as f:
            f.write('residue,flux\n')
            for res, flux in enumerate(flux_dest):
                f.write(f'{res},{flux}\n')
    print('Done')
    return

def edge_transfer_matrix(matrix: np.ndarray) -> tuple[np.ndarray, list]:
    """
    Calculate the edge to edge transfer matrix used to select the most important
    residues based of the flux.

    Parameters
    ----------
        matrix: numpy.array
            2D array with the interactions of the system.
    
    Return
    ------
        edge_to_edge: numpy.array
            2D array with the edge_to_edge propensity values.
        esges : list
            List with the tuples containing the nodes in each edge.
    """
    edges, incidence = get_incidence_matrix(matrix)
    dweights = get_weights(matrix, edges)

    # Calculate the laplacian matrix
    laplacian_matrix = incidence @ dweights
    laplacian_matrix = laplacian_matrix @ incidence.T

    # Monroe-Penrose pseudoinverse
    pseudoinverse = np.linalg.pinv(laplacian_matrix)

    # Edge-To-Edge
    edge_to_edge = dweights @ incidence.T
    edge_to_edge = edge_to_edge @ pseudoinverse
    edge_to_edge = edge_to_edge @ incidence

    return edge_to_edge, edges

def get_incidence_matrix(matrix: np.ndarray) -> tuple[list, np.ndarray]:
    """
    Construct the incidence matrix of the interaction graph

    Parameters
    ----------
        matrix : numpy.array
            2D numpy array matrix of the graph.
    
    Return
    ------
        edges : list
            List of tuples (node1,node2) containing the nodes involved in the
            each edge.
        incidence_matrix : np.ndarray
            2D Numpy array with the incidence matrix. The direction is assign 
            from the node with the lower index to the highest.
    """
    # Get edges from the matrix
    edges = []
    connections = matrix.nonzero()
    for node1, node2 in zip(connections[0], connections[1]):
        if (node1,node2) not in edges and (node2,node1) not in edges:
            edges.append((node1,node2))
    
    incidence_matrix = np.zeros((len(matrix), len(edges)))
    for edge, vertices in enumerate(edges):
        incidence_matrix[vertices[0]][edge] = 1
        incidence_matrix[vertices[1]][edge] = -1
    return edges, incidence_matrix

def get_weights(matrix: np.ndarray, edges: list) -> np.ndarray:
    """
    Get diagonal matrix of weights

    Parameters
    ----------
        matrix : numpy.array
            2D interaction matrix
        edges : list
            List with the tuples containing the nodes in each edge.

    Return
    ------
        dweights : numpy.array
            Diagonal matrix with the weights of each edge.
    """
    weights = np.zeros(len(edges))
    for edge, vertices in enumerate(edges):
        weights[edge] = matrix[vertices[0]][vertices[1]]
    dweights = np.diag(weights)
    return dweights

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

def score(coupling_path, matrices_path, nres, batch, interactions):

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
    print(f'Coupling path:                  {coupling_path}')
    print(f'No. of couplings file:          {len(coupling_stab)}')
    print(f'No. Batches:                    {len(batches)}')
    print(f'Interactions:                   {interactions}')
    

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

            es_files = []
            pts_files = []
            for interaction in interactions:
                es_files.extend(glob.glob(f'{matrices_path}/es_{interaction}.{job}.pickle'))
                pts_files.extend(glob.glob(f'{matrices_path}/pts_{interaction}.{job}.pickle'))

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