#!/bin/python

"""
@authors: Aliaa Abd Elhalim [aliaa.abdelhalim@edu.fh-joanneum.at]
          Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sanchez-Murcia [pedro.murcia@medunigraz.at]
"""
# Imports ======================================================================
# PyQt5
from PyQt5 import uic
from PyQt5.QtGui import QKeySequence
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from PyQt5.QtWidgets import (QMainWindow, QButtonGroup,QFileDialog, QShortcut)

# Plugin specific
from . import functions as fc

# Generals
import os
import numpy as np
import pandas as pd
from pymol import cmd

# for testing
# import functions as fc

# Simula Window ================================================================
class SimulaWindow(QMainWindow):
    def __init__(self, traj):
        super().__init__()
        uifile = os.path.join(os.path.dirname(__file__),
                              'ui_files/simula.ui')
        uic.loadUi(uifile, self) 

        # Init variables =======================================================
        self.conditions_dict = {'type':[], 'atoms':[], 'symbol': [], 'value':[]}
        self.cv_dict = {'type':[], 'atoms':[], 'coeff':[], 'fpath':[],
                        'harm':[]}
        self.traj = traj
        self.traj_file = ''
        self.top_file = ''
        self.frames_sel = []
        self.mask_list = []
        self.qmmm_list = ''

        # Init setup ===========================================================
        self.checkBox_all.setChecked(True)
        self.combo_cvtype.addItems(['DISTANCE', 'ANGLE', 'LCOD'])
        self.combo_theory.addItems(['DFTB3', 'EXTERN'])
        self.combo_filter.addItems(['distance', 'angle', 'dihedral'])
        self.combo_cond.addItems(['>', '<', '>=', '<='])
        self.combo_mask.addItems(cmd.get_names('selections', 0))
        self.edit_coeff.setEnabled(False)

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
        self.canvas1 = fc.MplCanvas(self, dpi=100)
        toolbar1 = NavigationToolbar2QT(self.canvas1, self)
        self.verticalLayout_14.addWidget(toolbar1)
        self.verticalLayout_14.addWidget(self.canvas1)
        self.ax1 = self.canvas1.fig.add_subplot(111)  
                
        # Connections ==========================================================
        # Tab sMD
        self.button_add1.clicked.connect(self.add_cv)
        self.button_clear1.clicked.connect(self.clear_cv)
        self.button_addmask.clicked.connect(self.add_mask)
        self.button_clearmask.clicked.connect(self.clear_mask)
        self.combo_cvtype.currentIndexChanged.connect(self.current_cv)

        # Tab Frames
        self.button_generate.clicked.connect(self.generate)
        self.checkBox_automatic.stateChanged.connect(self.onStateChanged)
        self.checkBox_manual.stateChanged.connect(self.onStateChanged)
        self.button_add2.clicked.connect(self.add_filter)
        self.button_clear2.clicked.connect(self.clear_filters)
        self.button_frames.clicked.connect(self.select_frames)
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
    
    # Tab sMD
    def current_cv(self):
        if self.combo_cvtype.currentText() != 'LCOD':
            self.edit_coeff.setEnabled(False)
            self.edit_path.setPlaceholderText("Final value of the CV")
        else:
            self.edit_coeff.setEnabled(True)
            self.edit_path.setPlaceholderText("Not required")

    def add_cv(self):
        cv_type = self.combo_cvtype.currentText()
        coeff_flag = False
        fpath_flag = True
        atoms_sel = []

        selection_raw = self.edit_atoms1.text()
        try:
            for atom in selection_raw.strip().split():
                atoms_sel.append(int(atom) - 1)
            atoms_sel = np.asarray(atoms_sel)

            if cv_type == 'DISTANCE' and len(atoms_sel) != 2:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            elif cv_type == 'ANGLE' and len(atoms_sel) != 3:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            elif cv_type == 'DIHEDRAL' and len(atoms_sel) != 4:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            elif cv_type == 'LCOD' and len(atoms_sel)%2 != 0:
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
                fc.pop_error("Error!!!", "Could not read final path")
                return

        try:
            coeff = np.asarray(self.edit_coeff.text().strip().split(),
                               dtype=float)
            if cv_type == 'LCOD' and len(coeff) != len(atoms_sel)//2:
                fc.pop_error("Error!!!", "Wrong number of coefficients")
            elif cv_type == 'LCOD' and len(coeff) == len(atoms_sel)//2:
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
        self.cv_dict['atoms'].append(atoms_sel + 1)
        self.cv_dict['fpath'].append(fpath)
        self.cv_dict['coeff'].append(coeff)
        self.cv_dict['harm'].append(harm)

        # Print collective variables
        self.edit_show1.insertPlainText(f'type = {cv_type}\n')
        self.edit_show1.insertPlainText(f'atoms = {atoms_sel + 1}\n')
        if fpath_flag:
            self.edit_show1.insertPlainText(f'final path = {fpath}\n')
        if coeff_flag:
            self.edit_show1.insertPlainText(f'Coefficients = {coeff}\n')
        self.edit_show1.insertPlainText(f'Harmonic = {harm}\n/\n')
        
        # Clean edits
        self.edit_path.setText("")
        self.edit_coeff.setText("")
        self.edit_harm.setText("")

        # Show in pymol
        if cv_type == 'DISTANCE':
            cmd.distance(f'cv_{len(self.cv_dict['type'])}',
                         f'index {atoms_sel[0] + 1}',
                         f'index {atoms_sel[1] + 1}')
        elif cv_type == 'ANGLE':
            cmd.angle(f'cv_{len(self.cv_dict['type'])}',
                         f'index {atoms_sel[0] + 1}',
                         f'index {atoms_sel[1] + 1}',
                         f'index {atoms_sel[2] + 1}')
        elif cv_type == 'DIHEDRAL':
            cmd.dihedral(f'cv_{len(self.cv_dict['type'])}',
                         f'index {atoms_sel[0] + 1}',
                         f'index {atoms_sel[1] + 1}',
                         f'index {atoms_sel[2] + 1}',
                         f'index {atoms_sel[3] + 1}')
        elif cv_type == 'LCOD':
            cmd.distance(f'cv_{len(self.cv_dict['type'])}_1',
                         f'index {atoms_sel[0] + 1}',
                         f'index {atoms_sel[1] + 1}')
            cmd.distance(f'cv_{len(self.cv_dict['type'])}_2',
                         f'index {atoms_sel[2] + 1}',
                         f'index {atoms_sel[3] + 1}')
    
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
        self.edit_show2.insertPlainText(f'Mask:')
        self.edit_show2.insertPlainText(f'{qmmm_mask}\n')

    def clear_cv(self):
        self.cv_dict = {'type':[], 'atoms':[], 'coeff':[], 'fpath':[],
                        'harm':[]}
        self.edit_show1.clear()
        cmd.delete('cv_*')
    
    def clear_mask(self):
        self.mask_list = []
        self.edit_show2.clear()
        cmd.delete('mask_*')

    # Tab Frames
    def onStateChanged(self):
        if self.checkBox_automatic.isChecked():
            # Disable manual
            self.edit_frames.setEnabled(False)
            # Enable automatic
            self.combo_filter.setEnabled(True)
            self.edit_atoms2.setEnabled(True)
            self.edit_cond.setEnabled(True)
            self.combo_cond.setEnabled(True)
            self.button_add2.setEnabled(True)
            self.button_clear2.setEnabled(True)
        elif self.checkBox_manual.isChecked():
            # Disable automatic
            self.button_add2.setEnabled(False)
            self.button_clear2.setEnabled(False)
            self.combo_filter.setEnabled(False)
            self.edit_atoms2.setEnabled(False)
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
        selection_raw = self.edit_atoms2.text()
        atoms_sel = []
        try:
            for atom in selection_raw.strip().split():
                atoms_sel.append(int(atom) - 1)
            atoms_sel = np.asarray(atoms_sel)

            if filter_type == 'distance' and len(atoms_sel) != 2:
                fc.pop_error('Error!!!', 'Wrong number of atoms')
                return
            elif filter_type == 'angle'and len(atoms_sel) != 3:
                fc.pop_error('Error!!!', 'Wrong number of atoms')
                return
            elif filter_type == 'dihedral'and len(atoms_sel) != 4:
                fc.pop_error('Error!!!', 'Wrong number of atoms')
                return
        except:
            fc.pop_error("Error!!!", "Could not read atoms")
            return
        self.edit_atoms2.setText("")

        if condition_value == '':
            fc.pop_error("Error!!!", "No condition have been provided")
            return
        self.edit_cond.setText("")
        
        # Save conditions
        self.conditions_dict['type'].append(filter_type)
        self.conditions_dict['atoms'].append(atoms_sel)
        self.conditions_dict['symbol'].append(condition)
        self.conditions_dict['value'].append(condition_value)

        self.edit_show3.insertPlainText(f"Filter: {filter_type}\n")
        self.edit_show3.insertPlainText(f"Atoms: {atoms_sel + 1}\n")
        self.edit_show3.insertPlainText(f"Condition: {condition}")
        self.edit_show3.insertPlainText(f"{condition_value}\n/\n")
    
    def select_frames(self):    
        self.traj, need_top = fc.load_traj(self.traj_file, self.top_file)
        self.frames_sel = []
        self.label_selected.setText("")
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
        if len(self.frames_sel) == 0:
            fc.pop_error("Error!!!", "No frames selected")
            return
        
        for frame in self.frames_sel:
            os.mkdir(f'{self.output_dir}/qmmm_{frame}')
            outdir = f'{self.output_dir}/qmmm_{frame}'
            fc.save_rst(self.traj, frame, outdir)
            fc.write_cv(self.traj, frame, self.cv_dict, outdir)
        fc.write_run(self.output_dir)
        fc.write_qmmm(qmmm_dict, self.output_dir)
        fc.write_sh(self.frames_sel, self.output_dir)
        fc.pop_message("Info", "Files have been generated")
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

