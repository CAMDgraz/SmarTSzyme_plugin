#!/bin/python

"""
@author: Aliaa []
         Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
         Pedro A. Sanchez-Murcia [pedro.murcia@medunigraz.at]
"""

# Imports ======================================================================
# PyQt5
from matplotlib.figure import Figure
from PyQt5.QtWidgets import (QMessageBox, QFileDialog)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

# General
import os
import heapq
import numpy as np
import mdtraj as md
import pandas as pd
from pymol import cmd

# MPL canvas ===================================================================
class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, dpi=100): 
        self.fig = Figure(dpi=100)
        super().__init__(self.fig)
        

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
        return None, None
    
    if need_for_top(traj_file) and not top_file:
        traj_ext = os.path.basename(traj_file).split('.')[-1]
        pop_error('Topology file',
            f'Trajectory file with extension {traj_ext} needs a topology file')
        return None, None

    elif need_for_top(traj_file) and top_file:
        traj = md.load(traj_file, top=top_file)
        traj.center_coordinates()
        return traj, True

    elif not need_for_top(traj_file):
        traj = md.load(traj_file)
        traj.center_coordinates()
        return traj, False

# Selection ====================================================================
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

def get_mask_sel(mask_str, selection, idx):
    qmmm_mask = []
    mask_selection = 'index '
    residues = set([atom.resi for atom in cmd.get_model(f'{mask_str}').atom])

    if selection == 'all':
        mask_selection = 'resid '
        for resid in residues:
            qmmm_mask_ = f':{str(resid)}'
            qmmm_mask.append(qmmm_mask_)
            mask_selection += f'{resid}+'

    elif selection == 'backbone':
        for resid in residues:
            qmmm_mask_ = f':{str(resid)}@'
            for atom in cmd.get_model(f'resid {resid}').atom:
                if atom.name in ['N', 'H', 'CA', 'HA', 'C', 'O']:
                    mask_selection += f'{str(atom.id)}+'
                    qmmm_mask_ += f'{atom.name},'
            qmmm_mask.append(qmmm_mask_)
        
    elif selection == 'sidechain':
        for resid in residues:
            qmmm_mask_ = f':{str(resid)}@'
            for atom in cmd.get_model(f'resid {resid}').atom:
                if atom.name not in ['N', 'H', 'CA', 'HA', 'C', 'O']:
                    mask_selection += f'{str(atom.id)}+'
                    qmmm_mask_ += f'{atom.name},'
        qmmm_mask.append(qmmm_mask_)
    
    cmd.select(f'mask_{idx}', f'{mask_selection}')
    cmd.show('sticks', f'mask_{idx}')

    return qmmm_mask

def validate_sel(traj, selection, natoms=None):
    sel_idx = atoms_sel(traj, selection)
    if sel_idx.any():
        if natoms:
            if len(sel_idx) != natoms:
                pop_error("Selection Error",
                        f"The number of atoms is incorrect")
                sel_idx = []
    return sel_idx

# Plots ========================================================================
def measure(type, traj, atoms):
    if type == 'distance':
        atoms = np.reshape(atoms, (-1, 2))
        measure = md.compute_distances(traj, atoms)
        measure *= 10
    elif type == 'angle':
        atoms = np.reshape(atoms, (-1, 3))
        measure = md.compute_angles(traj, atoms)
        measure = np.rad2deg(measure)
    elif type == 'dihedral':
        atoms = np.reshape(atoms, (-1, 4))
        measure = md.compute_dihedrals(traj, atoms)
        measure = np.rad2deg(measure)
    elif type == 'rmsd':
        measure = md.rmsd(traj, traj, frame=0, atom_indices=atoms)
        measure *= 10

    return measure

# Output =======================================================================
def write_cv(traj, frame, cv_dict, outdir):
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
                atoms = cv_dict['atoms'][cv_id] - 1
                ipath = measure('distance', traj[frame - 1],
                                atoms)
            elif cv == 'ANGLE':
                atoms = cv_dict['atoms'][cv_id] - 1
                ipath = measure('angle', traj[frame - 1],
                                atoms)
            elif cv == 'LCOD':
                ipath = 0
                f.write(f' cv_nr={len(cv_dict['coeff'][cv_id])}\n')
                f.write(f' cv_r=')
                for r_id, r in enumerate(cv_dict['coeff'][cv_id]):
                    atoms = cv_dict['atoms'][cv_id][r_id*2:r_id*2+2] - 1
                    measure_ = measure('distance', traj[frame - 1],
                                       atoms)
                    ipath += r*measure_
                    f.write(f'{r},')
                f.write(f'\n')
            if not cv_dict['fpath'][cv_id]:
                f.write(f' path={ipath[0][0]:.4f},{-ipath[0][0]:.4f},\n')
            else:
                f.write(f' path={ipath[0][0]:.4f},{cv_dict['fpath'][cv_id]:.4f},\n')
            f.write(f' NHARM=1,\n')
            f.write(f' HARM={cv_dict['harm'][cv_id]},\n/\n')
                
def write_qmmm(qmmm_dict, outdir):
    with open(f'{outdir}/qmmm.in', 'w') as f:
        f.write('QMMM sMD\n')
        f.write('&cntrl\n')
        f.write(f' ntx=1,\n')                                       
        f.write(f' irest=0,\n')                                   
        f.write(f' ntxo=1,\n')                                      
        f.write(f' ntpr=100,\n')                                    
        f.write(f' ntwx=100,\n')                                     
        f.write(f' ntwv=-1,\n')                                     
        f.write(f' ntf=1,\n')                                     
        f.write(f' ntb=2,\n')                                    
        f.write(f' dielc=1.0,\n')                                   
        f.write(f' cut=10.,\n')                                     
        f.write(f' nsnb=10,\n')                                     
        f.write(f' imin=0,\n')                                      
        f.write(f' ibelly=0,\n')                                    
        f.write(f' iwrap=1,\n')                                     
        f.write(f' nstlim={qmmm_dict['steps']},\n')                                
        f.write(f' dt={qmmm_dict['timestep']},\n')                                   
        f.write(f' temp0=300.0,\n')                                 
        f.write(f' tempi=300.0,\n')                                 
        f.write(f' ntt=3,\n')                                       
        f.write(f' gamma_ln=1.0,\n')                                  
        f.write(f' vlimit=20.0,\n')                                 
        f.write(f' ntp=1,\n')                                 
        f.write(f' ntc=1,\n')                                       
        f.write(f' tol=0.00001,\n')                            
        f.write(f' pres0=1,\n')                                       
        f.write(f' comp=44.6,\n')                                     
        f.write(f' jfastw=0,\n')                                      
        f.write(f' nscm=1000,\n')                              
        f.write(f' ifqnt=1,\n')                                       
        f.write(f' infe=1,\n')                                        
        f.write(f'/\n')
        f.write(f'&qmmm\n')
        f.write(f" qmmask='")
        for mask in qmmm_dict['mask']:
            if mask == qmmm_dict['mask'][-1]:
                f.write(f"{mask}")
            else:
                f.write(f"{mask}|")
        f.write(f"',\n")
        f.write(f' qmcharge={qmmm_dict['charge']},\n')                                    
        f.write(f" qm_theory='{qmmm_dict['theory']}',\n")                             
        f.write(f' qmshake=0,\n')                                    
        f.write(f' writepdb=1,\n')                                     
        f.write(f' verbosity=0,\n')                                    
        f.write(f' qmcut=10.,\n')                                         
        f.write(f' printcharges=1,\n')                               
        f.write(f' printdipole=1,\n')                                
        f.write(f' peptide_corr=0,\n')                               
        if qmmm_dict['theory'] == 'DFTB3':
            f.write(f' dftb_telec=100,\n')
            f.write(f" dftb_slko_path='$AMBERHOME/dat/slko/3ob-3-1',\n")
        elif qmmm_dict['theory'] == 'EXTERN':
            f.write(f' qm_ewald=0,\n') 
        f.write(f'/\n')         
        f.write(f'&smd\n')                                                 
        f.write(f" output_file='smd_qmmm.txt'\n")                      
        f.write(f' output_freq=50\n')                                  
        f.write(f" cv_file='cv.in'\n")                                 
        f.write(f'/\n')
        if qmmm_dict['theory'] == 'EXTERN':
            f.write(f'&EXTERNTHEORY\n')
            f.write(f' Please write the corresponding parameters here\n/')

def write_sh(frames, outdir):
    with open(f'{outdir}/smd_jobs.txt', 'w') as f:
        for frame in frames:
            if frame != frames[-1]:
                f.write(f'{outdir}/qmmm_{frame}\n')
            elif frame == frames[-1]:
                f.write(f'{outdir}/qmmm_{frame}')

    with open(f'{outdir}/prepare.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('while IFS= read -r line; do\n')
        f.write('cp qmmm.in ${line}/\n')
        f.write('cp run_job.sh ${line}/\n')
        f.write('done < smd_jobs.txt\n')

    with open(f'{outdir}/run.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('while IFS= read -r line; do\n')
        f.write('cd ${line}/\n')
        f.write('frame=$(echo "${line}" | grep -o "[^_]*$")\n')
        f.write('sbatch run_job.sh ${frame}\n')
        f.write('cd ..\ndone < smd_jobs.txt\n')

def write_run(qmmm_dict, outdir):
    with open(f'{outdir}/run_job.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write(f'export amberpath={qmmm_dict['amberpath']}\n')
        f.write(f'export openmpi={qmmm_dict['mpipath']}\n')
        f.write('source ${amberpath}/amber.sh\n')
        f.write('export SANDER=${amberpath}/bin/sander.MPI\n')
        f.write('export PATH={openmpi}/bin${PATH:+:${PATH}}\n')
        f.write('export LD_LIBRARY_PATH=${openmpi}/')
        f.write('lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}\n\n')
        f.write('# Run sMD\n')
        f.write('mpirun -np 2 ${SANDER} -O -i qmmm.in -o frame_$1.out')
        f.write(' -p top_qmmm.top -c frame_$1.rst')
        f.write(' -r frame_$1.qmmm.rst -x traj_qmmm.nc')
        f.write(' -ref frame_$1.rst')

def save_rst(traj, frame, outdir):
    traj[frame - 1].save_amberrst7(f'{outdir}/frame_{frame}.rst')

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
    
    try:
        os.remove('results.pdb')
    except OSError:
        pass

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
    
        
    return 'results.pdb'

# sMD results ==================================================================
def exponential_average(works, time_vect, njobs, T, kb):
    """
    Perform exponential average

    Parameters:
        works : np.array
            Array containing the pulling works of the trajectories
        time_vect : np.array
            Array containing the time values (0:time:time_step)
        njobs : int
            Number of sMD trajectories
        T : float
            Temperature of the simulations
        kb : float
            Boltzmann constant
    Return:
        energy : np.array
            Array containing the energies resulting from JE
        avg_work : np.array
            Array containing the average work for each time step

    """
    avg_work = np.zeros(len(time_vect))
    exp_work = np.zeros(len(time_vect))
    for work_vect in works:
        exp = np.exp(-work_vect/(kb*T))
        exp_work += exp
        avg_work += work_vect

    avg_work = avg_work/njobs
    avg_exp_work = exp_work/njobs
    energy = -kb*T*np.log(avg_exp_work)

    return energy, avg_work


def cumulant_expansion(works, time_vect, njobs, T, kb):
    """
    Unbiased 2nd order cumulant expansion

    Parameters:
        works : np.array
            Array containing the pulling works of the trajectories
        njobs : int
            Number of sMD trajectories
        T : float
            Temperature of the simulations
        kb : float
            Boltzmann constant

    Return:
        energy : np.array
            Array containing the energies resulting from JE
        avg_work : np.array
            Array containing the average work for each time step

    """
    M = np.shape(works)[0]
    mean = np.mean(works, axis=0)                                   # <W>
    beta = 1/(kb*T)                                                 # ß
    mean_cuad = np.mean(works**2, axis=0)                           # <W^2>
    energy = mean - (beta/2)*(M/(M-1))*(mean_cuad - mean**2)        # F

    avg_work = np.zeros(len(time_vect))
    for work_vect in works:
        avg_work += work_vect

    avg_work = avg_work/njobs

    return energy, avg_work


def find_maximum(energy, time_vect):
    """
    Find energy maximum value using a gradient approach

    Parameters:
        energy : np.array
            Array containing the enrgies resulting from JE
        time_vect : np.array
            Array containing the time values (0:time:time_step)

    Return:
        m_energy : float
            Maximum energy value
        m_step : float
            Time step where the m_energy is found

    """
    maximum = []
    heapq.heapify(maximum)

    gradients = np.gradient(energy, time_vect)
    for grad in range(0, len(gradients) - 1):
        if (gradients[grad] > 0) and (gradients[grad+1] < 0):
            heapq.heappush(maximum, (energy[grad], grad))

    try:
        m_energy, m_step = heapq.nlargest(1, maximum)[0]
    except IndexError:
        m_energy, m_step = None, None
    return m_energy, m_step


def check_jobs(jobs, time, time_step):
    """
    Function for check incomplete or missing jobs

    Parameters:
        jobs : list
            List containing the path to the output files of the sMD
        time : float
            Simulation time in ps
        time_step : float
            Time step used in the simulations

    Return:
        incomplete_jobs : list
            List containing the path to the incomplete or missing output files
            of the sMD

    """
    length_QMMM = len(np.arange(0, time, time_step))

    incomplete_jobs = []
    for job in jobs:
        try:
            job_ = pd.read_csv(job, delim_whitespace=True, skiprows=3,
                               header=None, nrows=length_QMMM)
        except pd.errors.EmptyDataError:
            incomplete_jobs.append(job)
        except FileNotFoundError:
            incomplete_jobs.append(job)
        if len(job_) != length_QMMM:
            incomplete_jobs.append(job)

    return incomplete_jobs

# Errors and messages ==========================================================
def pop_error(title, message):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setWindowTitle(title)
        msg.setText(message)
        msg.exec_()

def pop_message(title, message):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setWindowTitle(title)
        msg.setText(message)
        msg.exec_()