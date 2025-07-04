#!/bin/python

"""
@authors: Aliaa Abd Elhalim [aliaa.abdelhalim@edu.fh-joanneum.at]
          Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sanchez-Murcia [pedro.murcia@medunigraz.at]
"""
import numpy as np
import pandas as pd
import heapq
import seaborn as sns

def write_cv(path, cv_dict, measure_dict, measures_pd, frame):
    cv_keys = cv_dict.keys()
    with open(path, "w") as f:
        pass

    for cv_key in cv_keys:
        cv_info = cv_dict[cv_key]
        with open(path, "a") as f:
            f.write(f"# CV: {cv_key}\n")
            f.write(f"&colvar\n")
            f.write(f"cv_type='{cv_info[0]}',\n")
            if cv_info[0] == 'LCOD':
                atoms = []
                for measure in cv_info[1]:
                    atoms_ = measure_dict[measure][2]
                    atoms.extend(atoms_)
            else:
                atoms = measure_dict[cv_info[1][0]][2]
            f.write(f"cv_ni={len(atoms)},\n")
            f.write(f"cv_i=")
            for atom in atoms:
                f.write(f"{atom},")
            if cv_info[0] == 'LCOD':
                f.write(f"\ncv_nr={len(cv_info[2])},\n")
                f.write(f"cv_r=")
                for coeff in cv_info[2]:
                    f.write(f"{coeff},")
            f.write(f"\nnpath=2,\n")
            f.write(f"path_mode='LINES',\n")
            if cv_info[0] == 'LCOD':
                coefficients = [int(coeff) for coeff in cv_info[2]]
                path_init = (measures_pd[cv_info[1]]*coefficients).sum(axis=1)
                path_init = np.asarray(path_init)[frame]
            else:
                path_init = np.asarray(measures_pd[cv_info[1][0]])[frame]
            f.write(f"path={path_init:.2f},{cv_info[3]},\n")
            f.write(f"NHARM=1,\n")
            f.write(f"HARM={cv_info[4]},\n/\n")
    return

def write_qmmm(path, qmmask, theory, steps, time_step, charge,
               amber_path, temp=300):
    with open(path, "w") as f:
        f.write(f"QMMM sMD\n")
        f.write(f"&cntrl\n")
        f.write(f" ntx=1,\n")
        f.write(f" irest=0,\n")
        f.write(f" ntxo=1,\n")
        f.write(f" ntpr=100,\n")
        f.write(f" ntwx=100,\n")
        f.write(f" ntwv=-1,\n")
        f.write(f" ntf=1,\n")
        f.write(f" ntb=2,\n")
        f.write(f" dielc=1.0,\n")
        f.write(f" cut=10.,\n")
        f.write(f" nsnb=10,\n")
        f.write(f" imin=0,\n")
        f.write(f" ibelly=0,\n")
        f.write(f" iwrap=1,\n")
        f.write(f" nstlim={steps},\n")
        f.write(f" dt={time_step},\n")
        f.write(f" temp0={temp},\n")
        f.write(f" tempi={temp},\n")
        f.write(f" ntt=3,\n")
        f.write(f" gamma_ln=1.0,\n")
        f.write(f" vlimit=20.0,\n")
        f.write(f" ntp=1,\n")
        f.write(f" ntc=1,\n")
        f.write(f" tol=0.00001,\n")
        f.write(f" pres0=1,\n")
        f.write(f" comp=44.6,\n")
        f.write(f" jfastw=0,\n")
        f.write(f" nscm=1000,\n")
        f.write(f" ifqnt=1,\n")
        f.write(f" infe=1,\n")
        f.write(f"/\n")
        f.write(f"&qmmm\n")
        f.write(f" qmmask='")
        for mask in qmmask:
            if mask != qmmask[-1]:
                f.write(f"{mask} |")
            else:
                f.write(f"{mask}")

        f.write(f"'\n qmcharge={charge},\n")
        f.write(f" qm_theory='{theory}',\n")
        f.write(f" qmshake=0,\n")
        f.write(f" writepdb=1,\n")
        f.write(f" verbosity=0,\n")
        f.write(f" qmcut=10.,\n")
        f.write(f" printcharges=1,\n")
        f.write(f" printdipole=1,\n")
        f.write(f" peptide_corr=0,\n")
        if theory == "DFTB3":
            f.write(f" dftb_telec=100,\n")
            f.write(f" dftb_slko_path='{amber_path}/dat/slko/3ob-3-1',\n")
        elif theory == "EXTERN":
            f.write(f" qm_ewald=0,\n")
        f.write(f"/\n")
        f.write(f"&smd\n")
        f.write(f" output_file='smd_qmmm.txt',\n")
        f.write(f" output_freq=50,\n")
        f.write(f" cv_file='cv.in',\n")
        f.write(f"/\n")
        if theory == "EXTERN":
            f.write(f"&EXTERNTHEORY\n")
            f.write(f" Write corresponding parameters here\n")
    return

def write_jobrun(path, amber_path, openmpi_path):
    with open(path, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write(f"export amberpath={amber_path}\n")
        f.write(f"export openmpi={openmpi_path}\n")
        f.write("source ${amber_path}/amber.sh\n\n")
        f.write("export PATH=${openmpi}/bin/${PATH:+:${PATH}}\n")
        f.write("export SANDER=${amberpath}/bin/sander.MPI\n")
        f.write("export MPIRUN=${openmpi}/bin/mpirun -np 4\n\n")
        f.write("$mpirun $SANDER -O -i qmmm.in -o qmmm.out -p top.top -c frame.rst -r qmmm.rst -x qmmm.nc -ref frame.rst\n")
    return

def smd_results(smd_list, time, time_step, type, temp, ax):
    print("Reading smd files ...")
    # Setup axes
    ax.set(ylabel=r'PMF $(kcal/mol)$', xlabel=r'Time $(ps)$')

    kb = 0.001982923700 # boltzman constant
    time_vect = np.arange(0, time, time_step)
    len_smd = len(time_vect)

    with open(smd_list, 'r') as f:
        jobs = f.read().splitlines()

    pulling_works = []
    for job in jobs:
        try:
            job_df = pd.read_csv(job, comment='#', nrows=len_smd, engine='python', 
                                header=None, sep=r'\s+')
            sns.lineplot(x=time_vect, y=job_df.iloc[:, -1], ax=ax, alpha=0.8)
        except:
            print(f'Skipping {job}')
            continue
        pulling_works.append(np.asarray(job_df.iloc[:, -1]))

    pulling_works = np.asarray(pulling_works)
    if type == "Cumulative":
        energy, avg_work = cumulant_expansion(pulling_works, time_vect,
                                              pulling_works.shape[0], temp, kb)
    elif type == "Exponential":
        energy, avg_work = exponential_average(pulling_works, time_vect,
                                              pulling_works.shape[0], temp, kb)
    
    m_energy, m_step = find_maximum(energy, time_vect)
    sns.lineplot(x=time_vect, y=energy, ax=ax, color='black')
    return m_energy, time_vect[m_step]

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
