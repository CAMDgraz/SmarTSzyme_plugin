"""
Created on Mon Nov 18 09:39:43 2024

@author: Aliaa []
         Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
         Pedro A. Sanchez-Murcia [pedro.murcia@medunigraz.at]
"""

import mdtraj as md
import os
import pandas as pd
import numpy as np
import seaborn as sns
import math

from PyQt5.QtWidgets import (QMessageBox, QFileDialog)


# Load systems =================================================================
valid_tops = set(['pdb', 'pdb.gz', 'h5', 'lh5', 'prmtop', 'parm7', 'prm7',
                  'psf', 'mol2', 'hoomdxml', 'gro', 'arc', 'hdf5', 'gsd'])
valid_trajs = set(['arc', 'dcd', 'binpos', 'xtc', 'trr', 'hdf5', 'h5', 'ncdf',
                   'netcdf', 'nc', 'pdb.gz', 'pdb', 'lh5', 'crd', 'mdcrd',
                   'inpcrd', 'restrt', 'rst7', 'ncrst', 'lammpstrj', 'dtr',
                   'stk', 'gro', 'xyz.gz', 'xyz', 'tng', 'xml', 'mol2',
                   'hoomdxml', 'gsd'])

def get_traj_top(type, edit_line):
        initial_dir = "./"
        ext_filter = "All files (*)"

        
        # Get path or open file dialog
        if edit_line.text() == "":
            file, _ = QFileDialog.getOpenFileName(None,
                                                  "Load Trajectory File",
                                                  initial_dir,
                                                  ext_filter)
            edit_line.setText(os.path.basename(file))
        else:
            file = edit_line.text()
        
        if file != "":
            if type == 'traj':
                # Check valid trajectry extension
                if is_valid_traj(file):
                    filename = file
                else:
                    edit_line.setText("")
                    return
            elif type == 'top':
                # Check valid trajectry extension
                if is_valid_top(file):
                    filename = file
                else:
                    edit_line.setText("")
                    return
        return filename

def is_valid_traj(traj, valid_trajs=valid_trajs):
    traj_ext = os.path.basename(traj).split('.')[-1] 
    if traj_ext not in valid_trajs:
        pop_error('Trajectory file',
                  f'Extension provided "{traj_ext}" is not valid')
        return False
    return True

def is_valid_top(top_file, valid_tops=valid_tops):
    top_ext = os.path.basename(top_file).split('.')[-1]
    if top_ext not in valid_tops:
        pop_error('Topology file',
                  f'Extension provided "{top_ext}" is not valid')
        return False
    return True

def need_for_top(traj_file):
    traj_ext = os.path.basename(traj_file).split('.')[-1]
    if traj_ext in ['h5', 'lh5', 'pdb']:
        return False
    return True

def load_traj(traj_file, top_file):
    if traj_file == '':
        pop_error('Trajectory file',
                  f'No Trajectory file supplied')
        return None
    
    if need_for_top(traj_file) and not top_file:
        traj_ext = os.path.basename(traj_file).split('.')[-1]
        pop_error('Topology file',
            f'Trajectory file with extension {traj_ext} needs a topology file')
        return None

    elif need_for_top(traj_file) and top_file:
        traj = md.load(traj_file, top=top_file)
        traj.center_coordinates()
        return traj

    elif not need_for_top(traj_file):
        traj = md.load(traj_file)
        traj.center_coordinates()
        return traj

# Atom and frame selection =====================================================
def atoms_sel(traj, selection):
    try:
        sel_idx = traj.topology.select(selection)
    except Exception:
        pop_error('Selection Error',
                  f'The selection "{selection}" is not valid in MDTraj')
        return np.asarray([])

    if sel_idx.size == 0:
        pop_error('Selection Error',
                  f'The provided selection correspond to no atoms')
        return False

    return sel_idx

# Plots ========================================================================
def validate_sel(traj, edit_line, natoms=None):
    selection = f'index {edit_line.text()}'
    sel_idx = atoms_sel(traj, selection)
    if sel_idx.any():
        if natoms:
            if len(sel_idx) != natoms:
                pop_error("Selection Error",
                        f"The number of atoms is incorrect")
    return sel_idx

def measure(type, traj, atoms):
    if type == 'distance':
        atoms = np.reshape(atoms, (-1, 2))
        measure = md.compute_distances(traj, atoms)
        measure *= 10
    elif type == 'angle':
        atoms = np.reshape(atoms, (-1, 3))
        measure = md.compute_angles(traj, atoms)
    elif type == 'dihedral':
        atoms = np.reshape(atoms, (-1, 4))
        measure = md.compute_dihedrals(traj, atoms)
    elif type == 'rmsd':
        measure = md.rmsd(traj, traj, frame=0, atom_indices=atoms,
                            precentered=True)
        measure *= 10

    return measure

def range_traj(traj, first, last, stride):
    n_frames = traj.n_frames
    first_range = range(0, n_frames)
    last_range = range(first + 1, n_frames+1)

    if last is not None:
        stride_range = range(1, last - first)
    else:
        stride_range = range(1, n_frames - first)

    if first not in first_range:
        raise ValueError('\n\n>>> First frame must be in the interval'
                         ' [{}, {}]'.format(first_range.start + 1,
                                            first_range.stop))
    if last and (last not in last_range):
        raise ValueError('\n\n>>> Last frame must be in the interval'
                         ' [{}, {}]'.format(last_range.start + 1,
                                            last_range.stop))
    if stride not in stride_range:
        raise ValueError('\n\n>>> Stride must be in the interval'
                         ' [{}, {}]'.format(stride_range.start,
                                            stride_range.stop))
    sliced = slice(first, last, stride)
    if sliced not in [slice(0, n_frames, 1), slice(0, None, 1)]:
        return traj.slice(sliced)

    return traj

# QM/MM ========================================================================

# Output =======================================================================
def write_cv(traj, frame_id, cv_dict, outdir):
    with open(f'{outdir}/cv.in', 'w') as f:
        for cv_id, cv in enumerate(cv_dict['type']):
            f.write('&colvar\n')
            f.write(f" cv_type='{cv_dict['type'][cv_id]}',\n")
            f.write(f' cv_ni={len(cv_dict['atoms'][cv_id])},\n')
            f.write(f' cv_i=')
            for atom in cv_dict['atoms'][cv_id]:
                f.write(f'{atom},')
            f.write(f'\n')
            f.write(f' npath=2,\n')
            f.write(f" path_mode='LINES'\n")
            if cv == 'DISTANCE':
                ipath = measure('distance', traj[frame_id], cv_dict['atoms'][cv_id])
            elif cv == 'ANGLE':
                ipath = measure('angle', traj[frame_id], cv_dict['atoms'][cv_id])
            elif cv == 'LCOD':
                ipath = 0
                f.write(f' cv_nr={len(cv_dict['coeff'])}\n')
                for r_id, r in enumerate(cv_dict['coeff'][cv_id]):
                    measure_ = measure('distance', traj[frame_id],
                                       cv_dict['atoms'][cv_id][r_id*2:r_id*2+2])
                    ipath += r*measure_
                    f.write(f' cv_r={r}\n') 

            f.write(f' path={ipath[0][0]},{cv_dict['fpath'][cv_id]},\n')
            f.write(f' NHARM=1,\n')
            f.write(f' HARM={cv_dict['harm'][cv_id]},\n/\n')
                
def write_qmmm(qmmm_dict, outdir):
    with open(f'{outdir}/qmmm.in', 'w') as f:
        f.write('&cntrl\n')
        f.write(f' ntx = 1,\n')                                       
        f.write(f' irest = 0,\n')                                   
        f.write(f' ntxo = 1,\n')                                      
        f.write(f' ntpr = 100,\n')                                    
        f.write(f' ntwx = 100,\n')                                     
        f.write(f' ntwv =-1,\n')                                     
        f.write(f' ntf = 1,\n')                                     
        f.write(f' ntb = 2,\n')                                    
        f.write(f' dielc = 1.0,\n')                                   
        f.write(f' cut = 10.,\n')                                     
        f.write(f' nsnb = 10,\n')                                     
        f.write(f' imin = 0,\n')                                      
        f.write(f' ibelly = 0,\n')                                    
        f.write(f' iwrap = 1,\n')                                     
        f.write(f' nstlim = {qmmm_dict['steps']},\n')                                
        f.write(f' dt = {qmmm_dict['timestep']},\n')                                   
        f.write(f' temp0 = 300.0,\n')                                 
        f.write(f' tempi = 300.0,\n')                                 
        f.write(f' ntt = 3,\n')                                       
        f.write(f' gamma_ln=1.0,\n')                                  
        f.write(f' vlimit = 20.0,\n')                                 
        f.write(f' ntp = 1,\n')                                 
        f.write(f' ntc = 1,\n')                                       
        f.write(f' tol = 0.00001,\n')                            
        f.write(f' pres0=1,\n')                                       
        f.write(f' comp=44.6,\n')                                     
        f.write(f' jfastw=0,\n')                                      
        f.write(f' nscm=1000,\n')                              
        f.write(f' ifqnt=1,\n')                                       
        f.write(f' infe=1,\n')                                        
        f.write(f' /\n')
        f.write(f' &qmmm\n')                   
        f.write(f" qmmask = '{qmmm_dict['mask']}',\n")
        f.write(f' qmcharge={qmmm_dict['charge']},\n')                                    
        f.write(f" qm_theory='{qmmm_dict['theory']}',\n")                             
        f.write(f' qmshake=0,\n')                                    
        f.write(f' writepdb=1,\n')                                     
        f.write(f' verbosity=0,\n')                                    
        f.write(f' qmcut=10.,\n')                                      
        f.write(f' dftb_telec=100,\n')                                 
        f.write(f' printcharges = 1,\n')                               
        f.write(f' printdipole = 1,\n')                                
        f.write(f' peptide_corr = 0,\n')                               
        f.write(f" dftb_slko_path='/usr/local/amber20/dat/slko/3ob-3-1',\n") 
        f.write(f' /\n')         
        f.write(f' &smd\n')                                                 
        f.write(f" output_file = 'smd_out.txt'\n")                      
        f.write(f' output_freq = 50\n')                                  
        f.write(f" cv_file= 'cv.in'\n")                                 
        f.write(f' /\n')

    


# Vizualisation ================================================================
def normalize_flux(csv_file):
    
    """
    Normalizes residue dependent flux values from csv file
    
    Parameters
    ----------
    csv file: str
        Path to csv file 
        contains residues (int), fluxes (float)
                    
    Return
    ------
    df['flux_norm']: DataFrame
        flux normalized with min-max normalization
        rounded to 4 decimal places
    """

    range_min = 0
    range_max = 1

    df = pd.read_csv(csv_file)
    df.loc[df['flux']!= 0, 'flux_norm'] = ((range_min) + ((df['flux'] - df['flux'].min()) * (range_max - (range_min))/ (df['flux'].max() - df['flux'].min()))).round(4)
    df.loc[df['flux_norm'].isna(), 'flux_norm'] = (range_max - (- range_min)) / 2 #middle value

    
    
    return df['flux_norm']
        
def read_pdb(pdb_file):
    
    """
    Reads pdb file as a nested list
    
    Parameters
    ----------
    pdb file: str
        Path to pdb file 
        
    Return
    ------
    pdb_list: nested list with specified coloumns
    """ 
    
    pdb_list = []
    with open(pdb_file) as pdbfile:
        for line in pdbfile:
            if line[:4] == 'ATOM' or line[:6] == "HETATM":
            
                # Split the line
                split_line = [line[:6],             # ATOM / HETATM (str)
                              int(line[6:11]),      # Atom serial number (int)
                              line[12:16],          # Atom name (str)
                              line[17:20],          # Residue name (str)
                              int(line[22:26]),     # Residue sequence number (int)
                              float(line[30:38]),   # X orthogonal Å coordinate (float)
                              float(line[38:46]),   # Y orthogonal Å coordinate (float)
                              float(line[46:54]),   # Z orthogonal Å coordinate (float)
                              float(line[54:60]),   # Occupancy (float)
                              float(line[60:66]),   # Temperature factor (float)
                              line[76:78]]          # Element symbol (str)                          
                
                pdb_list.append(split_line)
                
    return pdb_list             
    

def results_pdb(pdb_list,csv_flux):

    """
    Replaces beta-factors from original pdb 
    with previously created pdb nested list and calculated normalized flux values
    
    Parameters
    ----------
    pdb_list: nested list
        
    flux_norm: list of str
                
    Return
    ------
    mod_pdb: str 
        generates pdb file with modified beta-factors 
        depending on residue sequence number
    """
     
    for i in range(len(pdb_list)): 
        for j in range(len(csv_flux)):
            if j+1 == pdb_list[i][4]:
                pdb_list[i][9] = csv_flux[j]
    
    with open('results.pdb', 'w') as mod_pdb: 
        for i in range(len(pdb_list)):
            pdb_list[i][0] = pdb_list[i][0].ljust(6)
            pdb_list[i][1] = str(pdb_list[i][1]).rjust(5)
            pdb_list[i][2] = pdb_list[i][2].ljust(4)
            pdb_list[i][3] = pdb_list[i][3].rjust(3)
            pdb_list[i][4] = str(pdb_list[i][4]).rjust(6)
            pdb_list[i][5] = f"{pdb_list[i][5]:8.3f}".rjust(8)
            pdb_list[i][6] = f"{pdb_list[i][6]:8.3f}".rjust(8)
            pdb_list[i][7] = f"{pdb_list[i][7]:8.3f}".rjust(8)
            pdb_list[i][8] = f"{pdb_list[i][8]:6.2f}".rjust(6)
            pdb_list[i][9] = str(pdb_list[i][9]).rjust(6)
            pdb_list[i][10] = pdb_list[i][10].ljust(2)
            mod_pdb.write("%s%s %s %s%s    %s%s%s%s%s          %s\n"
                          % (pdb_list[i][0],pdb_list[i][1],pdb_list[i][2],
                             pdb_list[i][3],pdb_list[i][4],pdb_list[i][5],
                             pdb_list[i][6],pdb_list[i][7],pdb_list[i][8],
                             pdb_list[i][9],pdb_list[i][10]))
    
        
    return mod_pdb

# Errors =======================================================================
def pop_error(title, message):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setWindowTitle(title)
        msg.setText(message)
        msg.exec_()