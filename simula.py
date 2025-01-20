#!/bin/python

"""
@author: Aliaa []
         Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
         Pedro A. Sanchez-Murcia [pedro.murcia@medunigraz.at]
"""
# Imports
from PyQt5 import uic
from PyQt5.QtWidgets import (QWidget, QButtonGroup,QFileDialog, QLineEdit, QShortcut)
from PyQt5.QtGui import QKeySequence

from . import functions as fc
import numpy as np
import os
from pymol import cmd
import pandas as pd
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT


class SimulaWindow(QWidget):
    def __init__(self, traj):
        super().__init__()
        uifile = os.path.join(os.path.dirname(__file__),
                              'ui_files/simula_input.ui')
        uic.loadUi(uifile, self)   

        # Init variables =======================================================
        self.conditions_dict = {'type':[], 'atoms':[], 'symbol': [], 'value':[]}
        self.cv_dict = {'type':[], 'atoms':[], 'coeff':[], 'fpath':[],
                        'harm':[]}
        self.traj = traj
        self.traj_file = ''
        self.top_file = ''
        self.frames_sel = None
        self.mask_list = []
        self.qmmm_list = ''

        # Init setup ===========================================================
        self.checkBox_all.setChecked(True)
        self.combo_cvtype.addItems(['DISTANCE', 'ANGLE', 'LCOD'])
        self.combo_theory.addItems(['DFTB3', 'EXTERN'])
        self.combo_filter.addItems(['distance', 'angle', 'dihedral'])
        self.combo_cond.addItems(['>', '<', '>=', '<='])
        self.combo_mask.addItems(cmd.get_names('selections', 0))
        self.combo_atoms1.addItems(cmd.get_names('selections', 0))
        self.combo_atoms2.addItems(cmd.get_names('selections', 0))

        # Refresh shortcut
        self.refresh = QShortcut(QKeySequence('F5'), self)
        self.refresh.activated.connect(self.update_pymol_names)

        # Group check_boxes
        self.group_checked = QButtonGroup()
        self.group_checked.addButton(self.checkBox_automatic)
        self.group_checked.addButton(self.checkBox_manual)
        self.group_checked.setExclusive(True)

        self.group_select = QButtonGroup()
        self.group_select.addButton(self.checkBox_all)
        self.group_select.addButton(self.checkBox_backbone)
        self.group_select.addButton(self.checkBox_sidechain)
        self.group_select.setExclusive(True)  

        self.group_jarz = QButtonGroup()
        self.group_jarz.addButton(self.check_cumulant)
        self.group_jarz.addButton(self.check_exponential) 
        self.group_jarz.setExclusive(True)

        # Canvas
        self.canvas1 = fc.MplCanvas(self, width=5, height=4, dpi=100)
        toolbar1 = NavigationToolbar2QT(self.canvas1, self)
        self.verticalLayout_4.addWidget(toolbar1)
        self.verticalLayout_4.addWidget(self.canvas1)
        self.ax1 = self.canvas1.fig.add_subplot(111)  
                
        # Connections ==========================================================
        # Tab sMD
        self.button_add1.clicked.connect(self.add_cv)
        self.button_clear1.clicked.connect(self.clear_cv)
        self.button_addmask.clicked.connect(self.add_mask)
        self.button_clearmask.clicked.connect(self.clear_mask)

        # Tab Frames
        self.button_generate.clicked.connect(self.generate)
        self.checkBox_automatic.stateChanged.connect(self.onStateChanged)
        self.checkBox_manual.stateChanged.connect(self.onStateChanged)
        self.button_add2.clicked.connect(self.add_filter)
        self.button_clear2.clicked.connect(self.clear_filters)
        self.button_frames.clicked.connect(self.select_frames)
        self.button_clear_frames.clicked.connect(self.clear_frames)
        self.button_originalTop.clicked.connect(self.load_file)
        self.button_originalTraj.clicked.connect(self.load_file)
        self.button_output.clicked.connect(self.sel_output)

        # Tab sMD results
        self.button_qmmmlist.clicked.connect(self.get_qmmmlist)
        self.button_calculate.clicked.connect(self.calc_jarz)

    # Functions  ===============================================================
    # General
    def update_pymol_names(self):
        self.combo_mask.clear()
        self.combo_mask.addItems(cmd.get_names('selections', 0))
        self.combo_atoms1.clear()
        self.combo_atoms1.addItems(cmd.get_names('selections', 0))
        self.combo_atoms2.clear()
        self.combo_atoms2.addItems(cmd.get_names('selections', 0))
    
    # Tab sMD
    def add_cv(self):
        cv_type = self.combo_cvtype.currentText()
        coeff_flag = False
        fpath_flag = True

        selection_raw = self.combo_atoms1.currentText()
        try:
            selection_str = 'index '
            if selection_raw in cmd.get_names('selections', 0):
                atoms = np.asarray([atom.id for atom in
                                    cmd.get_model(selection_raw).atom],
                                    dtype=int)
                for atom in atoms:
                    selection_str += f'{str(int(atom) - 1)} '
            else:
                for atom in selection_raw.strip().split():
                    selection_str += f'{str(int(atom) - 1)} '
            sel_idx = fc.validate_sel(self.traj,
                                      f'{selection_str}', None)

            if cv_type == 'DISTANCE' and len(sel_idx) != 2:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            elif cv_type == 'ANGLE' and len(sel_idx) != 3:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            elif cv_type == 'DIHEDRAL' and len(sel_idx) != 4:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            elif cv_type == 'LCOD' and len(sel_idx)%2 != 0:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return 
        except:
            fc.pop_error("Error!!!", "Could not read atoms")
            return

        try:
            fpath = float(self.edit_path.text().strip())
        except:
            if cv_type == 'LCOD':
                fpath = None
                fpath_flag = False
            else:
                fc.pop_error("Error!!!", "Could not read final value of the path")
                return

        try:
            coeff = np.asarray(self.edit_coeff.text().strip().split(),
                               dtype=float)
            if cv_type == 'LCOD' and len(coeff) != len(atoms)//2:
                fc.pop_error("Error!!!", "Wrong number of coefficients")
            elif cv_type == 'LCOD' and len(coeff) == len(atoms)//2:
                coeff_flag = True
        except:
            if cv_type != 'LCOD':
                coeff = None
            else:
                fc.pop_error("Error!!!", "Could not read the coefficient")
                return

        try:
            harm = float(self.edit_harm.text().strip())
        except:
            harm = 600
            

        self.cv_dict['type'].append(cv_type)
        self.cv_dict['atoms'].append(sel_idx + 1)
        self.cv_dict['fpath'].append(fpath)
        self.cv_dict['coeff'].append(coeff)
        self.cv_dict['harm'].append(harm)

        # Print collective variables
        self.edit_show1.insertPlainText(f'type = {cv_type}\n')
        self.edit_show1.insertPlainText(f'atoms = {atoms}\n')
        if fpath_flag:
            self.edit_show1.insertPlainText(f'final path = {fpath}\n')
        if coeff_flag:
            self.edit_show1.insertPlainText(f'Coefficients = {coeff}\n')
        self.edit_show1.insertPlainText(f'Harmonic = {harm}\n/\n')
        
        # Clean edits
        self.edit_path.setText("")
        self.edit_coeff.setText("")
        self.edit_harm.setText("")

    def add_mask(self):
        mask_str = self.combo_mask.currentText().strip()
        
        if self.checkBox_all.isChecked():
            qmmm_mask = fc.get_mask_sel(mask_str, 'all',
                                        len(self.mask_list) + 1)
        elif self.checkBox_backbone.isChecked():
            qmmm_mask = fc.get_mask_sel(mask_str, 'backbone',
                                        len(self.mask_list) + 1)
        elif self.checkBox_sidechain.isChecked():
            qmmm_mask = fc.get_mask_sel(mask_str, 'sidechain',
                                        len(self.mask_list) + 1)
        
        self.mask_list.extend(qmmm_mask)
        self.edit_show2.insertPlainText(f'Mask {len(self.mask_list)}:')
        self.edit_show2.insertPlainText(f'{qmmm_mask}\n')

    def clear_cv(self):
        self.cv_dict = {'type':[], 'atoms':[], 'coeff':[], 'fpath':[],
                        'harm':[]}
        self.edit_show1.clear()
    
    def clear_mask(self):
        self.mask_list = []
        self.edit_show2.clear()

    # Tab Frames
    def onStateChanged(self):
        if self.checkBox_automatic.isChecked():
            # Disable manual
            self.edit_frames.setEnabled(False)
            # Enable automatic
            self.combo_filter.setEnabled(True)
            self.combo_atoms2.setEnabled(True)
            self.edit_cond.setEnabled(True)
            self.combo_cond.setEnabled(True)
            self.button_add2.setEnabled(True)
            self.button_clear2.setEnabled(True)
        elif self.checkBox_manual.isChecked():
            # Disable automatic
            self.button_add2.setEnabled(False)
            self.button_clear2.setEnabled(False)
            self.combo_filter.setEnabled(False)
            self.combo_atoms2.setEnabled(False)
            self.edit_cond.setEnabled(False)
            self.combo_cond.setEnabled(False)
            # Enable manual
            self.edit_frames.setEnabled(True)
    
    def sel_output(self):
        self.output_dir = str(QFileDialog.getExistingDirectory(self,
                                                               'Select output dir'))
        self.edit_output.setText(self.output_dir)

    def load_file(self):
        '''
        Load topology and trajectory file
        '''
        button_clicked = self.sender().objectName()
        if button_clicked == 'button_originalTraj':
            self.traj_file = fc.get_traj_top('traj', self.edit_trajectory)
        elif button_clicked == 'button_originalTop':
            self.top_file = fc.get_traj_top('top', self.edit_topology)
    
    def add_filter(self):
        filter_type = self.combo_filter.currentText()
        condition = self.combo_cond.currentText()
        condition_value = self.edit_cond.text()

        selection_raw = self.combo_atoms2.currentText()
        try:
            selection_str = 'index '
            if selection_raw in cmd.get_names('selections', 0):
                atoms = np.asarray([atom.id for atom in cmd.get_model(selection_raw).atom], dtype=int)
                for atom in atoms:
                    selection_str += f'{str(int(atom) - 1)} '
            else:
                for atom in selection_raw.strip().split():
                    selection_str += f'{str(int(atom) - 1)} '

            if filter_type == 'distance':
                sel_idx = fc.validate_sel(self.traj,
                                          f'{selection_str}', 2)
            elif filter_type == 'angle':
                sel_idx = fc.validate_sel(self.traj,
                                          f'{selection_str}', 3)
            elif filter_type == 'dihedral':
                sel_idx = fc.validate_sel(self.traj,
                                          f'{selection_str}', 4)
        except:
            fc.pop_error("Error!!!", "Could not read atoms")
        
        

        if condition_value == '':
            fc.pop_error("Error!!!", "No value for condition have been provided")
            return
        self.edit_cond.setText("")
        
        # Save conditions
        self.conditions_dict['type'].append(filter_type)
        self.conditions_dict['atoms'].append(sel_idx)
        self.conditions_dict['symbol'].append(condition)
        self.conditions_dict['value'].append(condition_value)

        self.edit_show3.insertPlainText(f"Filter: {filter_type}\n")
        self.edit_show3.insertPlainText(f"Atoms: {sel_idx + 1}\n")
        self.edit_show3.insertPlainText(f"Condition: {condition}")
        self.edit_show3.insertPlainText(f"{condition_value}\n/\n")
    
    def select_frames(self):    
        self.traj, need_top = fc.load_traj(self.traj_file, self.top_file)
        if not self.traj:
            fc.pop_error("Error!!!", "Load the topology and trajectory files")
            return
        
        if self.checkBox_manual.isChecked():
            if self.edit_frames.text() == '':
                fc.pop_error("Error!!!", "No frames were selected")
                return
            self.frames_sel = self.edit_frames.text().strip().split()
            self.frames_sel = np.asarray(self.frames_sel, dtype=int)

        elif self.checkBox_automatic.isChecked():
            self.frames_sel = np.arange(self.traj.n_frames)
            for condition in range(len(self.conditions_dict['type'])):
                atoms_ = self.conditions_dict['atoms'][condition]
                value_ = float(self.conditions_dict['value'][condition])
                measure = fc.measure(self.conditions_dict['type'][condition],
                                     self.traj, atoms_)
                if self.conditions_dict['symbol'][condition] == '>':
                    frames_ = np.where(measure > value_)[0]
                elif self.conditions_dict['symbol'][condition] == '<':
                    frames_ = np.where(measure < value_)[0]
                elif self.conditions_dict['symbol'][condition] == '>=':
                    frames_ = np.where(measure >= value_)[0]
                elif self.conditions_dict['symbol'][condition] == '<=':
                        frames_ = np.where(measure <= value_)[0]
                self.frames_sel = np.intersect1d(self.frames_sel, frames_)
                
            self.frames_sel += 1
        self.label_selected.setText(f"{len(self.frames_sel)} frames selected")
    
    def clear_filters(self):
        self.conditions_dict = {'type':[], 'atoms':[], 'symbol': [], 'value':[]}
        self.edit_show3.clear()
        
    def clear_mask(self):
        self.mask_list = []
        self.edit_show2.clear()
    
    def clear_frames(self):
        self.label_selected.setText("")
        self.fames_sel = None
        
    def generate(self):
        # Get qmmm.in params
        try:
            charge = int(self.edit_charge.text().strip())
        except:
            fc.pop_error('Error!!!', "Could not read the charge")
        
        try:
            steps = int(self.edit_steps.text().strip())
        except:
            steps = 15000
        
        try:
            timestep = float(self.edit_timestep.text().strip())
        except:
            timestep = 0.0002
                
        qmmm_dict = {'mask': self.mask_list,
                     'theory': self.combo_theory.currentText(),
                     'charge': charge,
                     'steps': steps,
                     'timestep': timestep}
        
        # Write cv.in file
        for frame_id, frame in enumerate(self.frames_sel):
            os.mkdir(f'{self.output_dir}/qmmm_{frame}')
            outdir = f'{self.output_dir}/qmmm_{frame}'
            fc.save_rst(self.traj, frame, outdir)
            fc.write_cv(self.traj, frame_id, self.cv_dict, outdir)
            fc.write_run(frame, outdir)
        fc.write_qmmm(qmmm_dict, self.output_dir)
        fc.write_sh(self.frames_sel, self.output_dir)
        fc.pop_error("Info", "Generation Terminated")
        return

    # Tab sMD results
    def get_qmmmlist(self):
        initial_dir = './'
        ext_filter = "All files (*)"

        file, _ = QFileDialog.getOpenFileName(None,
                                                  "Load File",
                                                  initial_dir,
                                                  ext_filter)
        self.qmmm_list = file
    
    def calc_jarz(self):
        kb = 0.001982923700 # boltzman constant
        if self.qmmm_list == '':
            fc.pop_error("Error!!!", "No qmmm list provided")
            return
        try:
            temp = int(self.edit_temp.text().strip().split()[0])
        except:
            temp = 300
        
        with open(self.qmmm_list, 'r') as inp:
            jobs = inp.read().splitlines()
        job_df = pd.read_csv(jobs[0], skiprows=3, skipfooter=3,
                             engine='python', header=None, sep=r'\s+')
        time_vect = np.asarray(job_df[job_df.columns[0]])

        pulling_works = np.zeros((len(jobs), len(time_vect)))

        for idx, job in enumerate(jobs):
            job_df = pd.read_csv(job, skiprows=3, skipfooter=3, engine='python',
                                 header=None, sep=r'\s+')

            pulling_works[idx] = np.asarray(job_df[job_df.columns[-1]])

        
        if self.check_cumulant.isChecked():
            energy, avg_work = fc.cumulant_expansion(pulling_works, time_vect,
                                                     len(jobs), temp, kb)
        elif self.check_exponential.isChecked():
            energy, avg_work = fc.exponential_average(pulling_works, time_vect,
                                                     len(jobs), temp, kb)
        m_energy, m_step = fc.find_maximum(energy, time_vect)
        if (m_energy is None):
            self.label_ts.setText('>> No maximum found :(')
        else:
            self.label_ts.setText(f'Maximum found at: {time_vect[m_step]:.2f} ps')
        
        self.ax1.clear()
        self.canvas1.draw()

        self.ax1.set_xlabel(r'Time $(ps)$')
        self.ax1.set_ylabel(r'Free Energy $(kcal\ mol^{-1})$')

        for work_vect in pulling_works:
            self.ax1.plot(time_vect, work_vect, alpha=0.2)
        self.ax1.plot(time_vect, energy, color='black')
        
        self.canvas1.draw()

