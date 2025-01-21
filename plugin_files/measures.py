#!/bin/python

"""
@author: Aliaa []
         Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
         Pedro A. Sanchez-Murcia [pedro.murcia@medunigraz.at]
"""

# Imports ======================================================================
# PyQt5
from PyQt5 import uic
from PyQt5.QtWidgets import (QWidget, QShortcut)
from PyQt5.QtGui import QKeySequence
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT

# Plugin specific
from . import functions as fc

# General
import os
import numpy as np
import pandas as pd
import seaborn as sns
from pymol import cmd

# for testing
#import functions as fc

# Measures Windows =============================================================
class MeasureWindow(QWidget):
    def __init__(self, traj):
        super().__init__()
        uifile = os.path.join(os.path.dirname(__file__),
                              'ui_files/measures.ui')
        uic.loadUi(uifile, self)

        # Init variables
        self.traj = traj
        self.out_path = './'
        self.data_dict = {'distance':[[], []], 'angle':[[], []],
                          'dihedral':[[], []], 'rmsd':[[], []]}

        # Canvas for each measure ==============================================
        # Distance
        self.canvas1 = fc.MplCanvas(self, dpi=100)
        toolbar1 = NavigationToolbar2QT(self.canvas1, self)
        self.verticalLayout_1.addWidget(toolbar1)
        self.verticalLayout_1.addWidget(self.canvas1)
        self.ax1 = self.canvas1.fig.add_subplot(111)

        # Angle
        self.canvas2 = fc.MplCanvas(self, dpi=100)
        toolbar2 = NavigationToolbar2QT(self.canvas2, self)
        self.verticalLayout_2.addWidget(toolbar2)
        self.verticalLayout_2.addWidget(self.canvas2)
        self.ax2 = self.canvas2.fig.add_subplot(111)

        # Dihedral
        self.canvas3 = fc.MplCanvas(self, dpi=100)
        toolbar3 = NavigationToolbar2QT(self.canvas3, self)
        self.verticalLayout_3.addWidget(toolbar3)
        self.verticalLayout_3.addWidget(self.canvas3)
        self.ax3 = self.canvas3.fig.add_subplot(111)

        # RMSD
        self.canvas4 = fc.MplCanvas(self, dpi=100)
        toolbar4 = NavigationToolbar2QT(self.canvas4, self)
        self.verticalLayout_4.addWidget(toolbar4)
        self.verticalLayout_4.addWidget(self.canvas4)
        self.ax4 = self.canvas4.fig.add_subplot(111)

        # Connections ==========================================================
        self.button_clear1.clicked.connect(self.clear)
        self.button_clear2.clicked.connect(self.clear)
        self.button_clear3.clicked.connect(self.clear)
        self.button_clear4.clicked.connect(self.clear)

        self.button_plot1.clicked.connect(self.add_plot)
        self.button_plot2.clicked.connect(self.add_plot)
        self.button_plot3.clicked.connect(self.add_plot)
        self.button_plot4.clicked.connect(self.add_plot)

        self.check_all.stateChanged.connect(self.state_checkbox)
        self.check_back.stateChanged.connect(self.state_checkbox)
        self.check_noH.stateChanged.connect(self.state_checkbox)

        self.button_csv1.clicked.connect(self.to_csv)
        self.button_csv2.clicked.connect(self.to_csv)
        self.button_csv3.clicked.connect(self.to_csv)
        self.button_csv4.clicked.connect(self.to_csv)

        # Refresh shortcut
        self.refresh = QShortcut(QKeySequence('F5'), self)
        self.refresh.activated.connect(self.update_pymol_names)

    # Functions ================================================================
    def update_pymol_names(self):
        self.combo_atoms1.clear()
        self.combo_atoms1.addItems(cmd.get_names('selections', 0))

        self.combo_atoms2.clear()
        self.combo_atoms2.addItems(cmd.get_names('selections', 0))

        self.combo_atoms3.clear()
        self.combo_atoms3.addItems(cmd.get_names('selections', 0))

        self.combo_atoms4.clear()
        self.combo_atoms4.addItems(cmd.get_names('selections', 0))

    def state_checkbox(self):
        sender = self.sender().objectName()
        if sender == 'check_all':
            if self.check_all.isChecked():
                self.check_back.setChecked(False)
                self.check_noH.setChecked(False)
        elif sender == 'check_back':
            if self.check_back.isChecked():
                self.check_all.setChecked(False)
                self.check_noH.setChecked(False)
        elif sender == 'check_noH':
            if self.check_noH.isChecked():
                self.check_all.setChecked(False)
                self.check_back.setChecked(False)
        return

    def add_plot(self):

        button_clicked = self.sender().objectName()

        # Distances
        if button_clicked == 'button_plot1':
            selection_raw = self.combo_atoms1.currentText()
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
                                          f'{selection_str}', 2)
            if len(sel_idx) != 2:
                return 
            
            label = self.edit_label1.text()
            if label == '':
                label = f"sel{len(self.data_dict['distance'][1]) + 1}"

            self.data_dict['distance'][0].append(sel_idx)
            self.data_dict['distance'][1].append(label)
        
            self.edit_label1.setText('')
            self.plot_measure('distance')


        # Angles
        if button_clicked == 'button_plot2':
            selection_raw = self.combo_atoms2.currentText()
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
                                      f'{selection_str}', 3)
            if len(sel_idx) != 3:
                return 
            
            label = self.edit_label2.text()
            if label == '':
                label = f"sel{len(self.data_dict['angle'][1]) + 1}"

            self.data_dict['angle'][0].append(sel_idx)
            self.data_dict['angle'][1].append(label)
        
            self.edit_label2.setText('')
            self.plot_measure('angle') 

        # Dihedrals
        if button_clicked == 'button_plot3':
            selection_raw = self.combo_atoms3.currentText()
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
                                      f'{selection_str}', 4)
            if len(sel_idx) != 4:
                return
            
            label = self.edit_label3.text()
            if label == '':
                label = f"sel{len(self.data_dict['dihedral'][1]) + 1}"

            self.data_dict['dihedral'][0].append(sel_idx)
            self.data_dict['dihedral'][1].append(label)
        
            self.edit_label3.setText('')
            self.plot_measure('dihedral')

        # RMSD
        if button_clicked == 'button_plot4':
            # check boxes
            if self.check_all.isChecked():
                selection_str = 'all'
            elif self.check_back.isChecked():
                selection_str = 'backbone'
            elif self.check_noH.isChecked():
                selection_str = 'name = noH '
            else:
                selection_raw = self.combo_atoms4.currentText()
                if selection_raw in cmd.get_names('selections', 0):
                    selection_str = 'index '
                    atoms = np.asarray([atom.id for atom in
                                        cmd.get_model(selection_raw).atom],
                                        dtype=int)
                    for atom in atoms:
                        selection_str += f'{str(int(atom) - 1)} '
                else:
                    for keyword in selection_raw.strip():
                        try:
                            selection_str += f'{int(keyword) - 1} '
                        except: 
                            selection_str += f'{keyword} '
            
            sel_idx = fc.validate_sel(self.traj, selection_str)
            if len(sel_idx) == 1:
                return 
            
            label = self.edit_label4.text()
            if label == '':
                label = f"sel{len(self.data_dict['rmsd'][1]) + 1}"

            self.data_dict['rmsd'][0].append(sel_idx)
            self.data_dict['rmsd'][1].append(label)
        
            self.edit_label4.setText('')
            
            self.check_all.setChecked(False)
            self.check_back.setChecked(False)
            self.check_noH.setChecked(False)
            self.plot_measure('rmsd')

    def plot_measure(self, sender):

        if sender == 'distance':
            ax=self.ax1
            canvas = self.canvas1
            ax.clear()
            canvas.draw()
        
            for atoms_, label_ in zip(self.data_dict['distance'][0],
                                      self.data_dict['distance'][1]):
                atoms_ = np.asarray(atoms_)
                atoms_ = np.reshape(atoms_, (-1, 2))
                measure = fc.measure('distance', self.traj, atoms_)
                ax.set_xlabel(r'Distance $(\AA)$')
                sns.kdeplot(x=measure.T[0], ax=ax, label=label_, fill=True)
                cmd.distance(label_, f'index {atoms_[0][0] + 1}',
                             f'index {atoms_[0][1] + 1}')

        elif sender == 'angle':
            ax=self.ax2
            canvas = self.canvas2
            ax.clear()
            canvas.draw()
        
            for atoms_, label_ in zip(self.data_dict['angle'][0],
                                      self.data_dict['angle'][1]):
                atoms_ = np.asarray(atoms_)
                atoms_ = np.reshape(atoms_, (-1, 3))
                measure = fc.measure('angle', self.traj, atoms_)
                ax.set_xlabel(r'Angle $(\AA)$')
                sns.kdeplot(x=measure.T[0], ax=ax, label=label_, fill=True)
                cmd.angle(label_, f'index {atoms_[0][0] + 1}',
                          f'index {atoms_[0][1] + 1}',
                          f'index {atoms_[0][2] + 1}')
        
        elif sender == 'dihedral':
            ax=self.ax3
            canvas = self.canvas3
            ax.clear()
            canvas.draw()
        
            for atoms_, label_ in zip(self.data_dict['dihedral'][0],
                                      self.data_dict['dihedral'][1]):
                atoms_ = np.asarray(atoms_)
                atoms_ = np.reshape(atoms_, (-1, 4))
                measure = fc.measure('dihedral', self.traj, atoms_)
                ax.set_xlabel(r'Dihedral $(\AA)$')
                sns.kdeplot(x=measure.T[0], ax=ax, label=label_, fill=True)
                cmd.dihedral(label_, f'index {atoms_[0][0] + 1}',
                             f'index {atoms_[0][1] + 1}',
                             f'index {atoms_[0][2] + 1}',
                             f'index {atoms_[0][3] + 1}')
        
        elif sender == 'rmsd':
            ax=self.ax4
            canvas = self.canvas4
            ax.clear()
            canvas.draw()
        
            for atoms_, label_ in zip(self.data_dict['rmsd'][0],
                                      self.data_dict['rmsd'][1]):
                atoms_ = np.asarray(atoms_)
                measure = fc.measure('rmsd', self.traj, atoms_)
                ax.set_ylabel(r'RMSD $(\AA)$')
                ax.set_xlabel('Frames')
                sns.lineplot(x=np.arange(self.traj.n_frames),
                             y=measure, ax=ax, label=label_)
        
        ax.legend(loc='upper right')
        canvas.draw()
        return
                
    def clear(self):
        button_clicked = self.sender().objectName()

        if button_clicked == 'button_clear1':
            self.data_dict['distance'] = [[], []]
            canvas = self.canvas1
            ax=self.ax1
        elif button_clicked == 'button_clear2':
            self.data_dict['angle'] = [[], []]
            canvas = self.canvas2
            ax=self.ax2
        elif button_clicked == 'button_clear3':
            self.data_dict['dihedral'] = [[], []]
            canvas = self.canvas3
            ax=self.ax3
        elif button_clicked == 'button_clear4':
            self.data_dict['rmsd'] = [[], []]
            canvas = self.canvas4
            ax=self.ax4
        ax.clear()
        canvas.draw()

    def to_csv(self):
        sender = self.sender().objectName()
        all_measures = []
        csv_file = ''
        labels = []

        if sender == 'button_csv1': 
            csv_file = 'distance.csv'       
            for atoms_, label_ in zip(self.data_dict['distance'][0],
                                      self.data_dict['distance'][1]):
                labels.append(label_)
                atoms_ = np.asarray(atoms_)
                atoms_ = np.reshape(atoms_, (-1, 2))
                measure = fc.measure('distance', self.traj, atoms_)
                all_measures.append(measure.T[0])

        elif sender == 'button_csv2':
            csv_file = 'angle.csv' 
            for atoms_, label_ in zip(self.data_dict['angle'][0],
                                      self.data_dict['angle'][1]):
                labels.append(label_)
                atoms_ = np.asarray(atoms_)
                atoms_ = np.reshape(atoms_, (-1, 3))
                measure = fc.measure('angle', self.traj, atoms_)
                all_measures.append(measure.T[0])
        
        elif sender == 'button_csv3':
            csv_file = 'dihedral.csv' 
            for atoms_, label_ in zip(self.data_dict['dihedral'][0],
                                      self.data_dict['dihedral'][1]):
                labels.append(label_)
                atoms_ = np.asarray(atoms_)
                atoms_ = np.reshape(atoms_, (-1, 4))
                measure = fc.measure('dihedral', self.traj, atoms_)
                all_measures.append(measure.T[0])
        
        elif sender == 'button_csv4':
            csv_file = 'rmsd.csv' 
            for atoms_, label_ in zip(self.data_dict['rmsd'][0],
                                      self.data_dict['rmsd'][1]):
                labels.append(label_)
                atoms_ = np.asarray(atoms_)
                measure = fc.measure('rmsd', self.traj, atoms_)
                all_measures.append(measure)

        all_measures_df = pd.DataFrame(np.asarray(all_measures).T,
                                       columns=labels)
        all_measures_df.to_csv(csv_file, index=False)