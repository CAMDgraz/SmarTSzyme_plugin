#!/bin/python

# Imports
from PyQt5 import uic
from PyQt5 import QtCore
from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QButtonGroup,
                             QFileDialog, QLineEdit)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure

import functions as fc
import seaborn as sns
import numpy as np
import os
from pymol import cmd


# General variables and functions ==============================================
class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)

# Main Window ==================================================================
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        uifile = os.path.join(os.path.dirname(__file__),
                              'ui_files/mainWindow.ui')
        uic.loadUi(uifile, self)
        self.show()

        # Init variables
        self.traj = None
        self.traj_file = ""
        self.top_file = ""

        # Connections
        self.actionMeasures.triggered.connect(self.show_Measures)
        self.actionSimula_Input.triggered.connect(self.show_Simula)
        #self.button_load.clicked.connect(self.load_system)
        self.button_trajectory.clicked.connect(self.load_file)
        self.button_topology.clicked.connect(self.load_file)

    # Functions
    def show_Measures(self):
        '''
        Shows the Measures window
        '''

        if not self.traj:
            fc.pop_error("Error!!!", "No system loaded")
            pass
        elif self.traj:
            self.measureWindow = MeasureWindow(self.traj)
            self.measureWindow.show()
    
    def show_Simula(self):
        '''
        Shows the Simula input window
        '''
        self.simulaWindow = SimulaWindow(self.traj)
        self.simulaWindow.show()
       
    def load_file(self):
        '''
        Load topology and trajectory file
        '''
        button_clicked = self.sender().objectName()
        if button_clicked == 'button_trajectory':
            self.traj_file = fc.get_traj_top('traj', self.edit_trajectory)
            
        elif button_clicked == 'button_topology':
            self.top_file = fc.get_traj_top('top', self.edit_topology)

    def load_system(self):
        self.traj = fc.load_traj(self.traj_file, self.top_file)      
        if self.traj:
            self.label_system.setText(str(self.traj))
            #system_name = os.path.splitext(os.path.basename(self.top_file))
            # pymol_objects = cmd.get_names('objects', 0)
            # if system_name not in pymol_objects:
            #     if not need_top:
            #         cmd.load(self.traj_file)
            #     else:
            #         cmd.load(self.top_file)
            #         cmd.load_traj(self.traj_file)
                

# Measures Window ==============================================================
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
        self.canvas1 = MplCanvas(self, width=5, height=4, dpi=100)
        toolbar1 = NavigationToolbar2QT(self.canvas1, self)
        self.verticalLayout_1.addWidget(toolbar1)
        self.verticalLayout_1.addWidget(self.canvas1)
        self.ax1 = self.canvas1.fig.add_subplot(111)

        # Angle
        self.canvas2 = MplCanvas(self, width=5, height=4, dpi=100)
        toolbar2 = NavigationToolbar2QT(self.canvas2, self)
        self.verticalLayout_2.addWidget(toolbar2)
        self.verticalLayout_2.addWidget(self.canvas2)
        self.ax2 = self.canvas2.fig.add_subplot(111)

        # Dihedral
        self.canvas3 = MplCanvas(self, width=5, height=4, dpi=100)
        toolbar3 = NavigationToolbar2QT(self.canvas3, self)
        self.verticalLayout_3.addWidget(toolbar3)
        self.verticalLayout_3.addWidget(self.canvas3)
        self.ax3 = self.canvas3.fig.add_subplot(111)

        # RMSD
        self.canvas4 = MplCanvas(self, width=5, height=4, dpi=100)
        toolbar4 = NavigationToolbar2QT(self.canvas4, self)
        self.verticalLayout_4.addWidget(toolbar4)
        self.verticalLayout_4.addWidget(self.canvas4)
        self.ax4 = self.canvas4.fig.add_subplot(111)

        # Connections
        self.button_add1.clicked.connect(self.add)
        self.button_add2.clicked.connect(self.add)
        self.button_add3.clicked.connect(self.add)
        self.button_add4.clicked.connect(self.add)

        self.button_remove1.clicked.connect(self.remove)
        self.button_remove2.clicked.connect(self.remove)
        self.button_remove3.clicked.connect(self.remove)
        self.button_remove4.clicked.connect(self.remove)

        self.button_clear1.clicked.connect(self.clear)
        self.button_clear2.clicked.connect(self.clear)
        self.button_clear3.clicked.connect(self.clear)
        self.button_clear4.clicked.connect(self.clear)

        self.button_plot1.clicked.connect(self.plot_measure)
        self.button_plot2.clicked.connect(self.plot_measure)
        self.button_plot3.clicked.connect(self.plot_measure)
        self.button_plot4.clicked.connect(self.plot_measure)

        #self.button_csv1.clicked.connect(self.to_csv)
        #self.button_csv2.clicked.connect(self.to_csv)
        #self.button_csv3.clicked.connect(self.to_csv)
        #self.button_csv4.clicked.connect(self.to_csv)

    def add(self):
        button_clicked = self.sender().objectName()
        
        # Distances
        if button_clicked == 'button_add1':
            sel_idx = fc.validate_sel(self.traj,
                                      f'index {self.edit_atom1.text()}', 2)
            label = self.edit_label1.text()
            if label == '':
                label = f"sel {len(self.data_dict['distance'][1]) + 1}"
            self.edit_show1.insertPlainText(f'Atoms {sel_idx} -- {label}\n')

            self.data_dict['distance'][0].append(sel_idx)
            self.data_dict['distance'][1].append(label)
        
            self.edit_atom1.setText('')
            self.edit_label1.setText('')

        # Angles
        if button_clicked == 'button_add2':
            sel_idx = fc.validate_sel(self.traj,
                                      f'index {self.edit_atom2.text()}', 3)
            label = self.edit_label2.text()
            if label == '':
                label = f"sel {len(self.data_dict['angle'][1]) + 1}"
            self.edit_show2.insertPlainText(f'Atoms {sel_idx} -- {label}\n')

            self.data_dict['angle'][0].append(sel_idx)
            self.data_dict['angle'][1].append(label)
        
            self.edit_atom2.setText('')
            self.edit_label2.setText('') 

        # Dihedrals
        if button_clicked == 'button_add3':
            sel_idx = fc.validate_sel(self.traj,
                                      f'index {self.edit_atom3.text()}', 4)
            label = self.edit_label3.text()
            if label == '':
                label = f"sel {len(self.data_dict['dihedral'][1]) + 1}"
            self.edit_show3.insertPlainText(f'Atoms {sel_idx} -- {label}\n')

            self.data_dict['dihedral'][0].append(sel_idx)
            self.data_dict['dihedral'][1].append(label)
        
            self.edit_atom3.setText('')
            self.edit_label3.setText('')

        # RMSD
        if button_clicked == 'button_add4':
            sel_idx = fc.validate_sel(self.traj, self.edit_atom4.text())
            label = self.edit_label4.text()
            if label == '':
                label = f"sel {len(self.data_dict['rmsd'][1]) + 1}"
            self.edit_show4.insertPlainText(f'Atoms {sel_idx} -- {label}\n')

            self.data_dict['rmsd'][0].append(sel_idx)
            self.data_dict['rmsd'][1].append(label)
        
            self.edit_atom4.setText('')
            self.edit_label4.setText('') 

    def remove(self):
        button_clicked = self.sender().objectName()

        if button_clicked == 'button_remove1':
            self.data_dict['distance'][0].pop(-1)
            self.data_dict['distance'][1].pop(-1)
            self.edit_show1.clear()
            if len(self.data_dict['distance'][0]) != 0:
                for atoms_, label_ in zip(self.data_dict['distance'][0],
                                          self.data_dict['distance'][1]):
                    self.edit_show1.insertPlainText(f'Atoms {atoms_} --')
                    self.edit_show1.insertPlainText(f'{label_}\n')
            self.ax1.clear()
            self.canvas1.draw()
            self.button_plot1.clicked.emit()

        elif button_clicked == 'button_remove2':
            self.data_dict['angle'][0].pop(-1)
            self.data_dict['angle'][1].pop(-1)
            self.edit_show2.clear()
            if len(self.data_dict['angle'][0]) != 0:
                for atoms_, label_ in zip(self.data_dict['angle'][0],
                                          self.data_dict['angle'][1]):
                    self.edit_show2.insertPlainText(f'Atoms {atoms_} --')
                    self.edit_show2.insertPlainText(f'{label_}\n')
            self.ax2.clear()
            self.canvas2.draw()
            self.button_plot2.clicked.emit()

        elif button_clicked == 'button_remove3':
            self.data_dict['dihedral'][0].pop(-1)
            self.data_dict['dihedral'][1].pop(-1)
            self.edit_show3.clear()
            if len(self.data_dict['dihedral'][0]) != 0:
                for atoms_, label_ in zip(self.data_dict['dihedral'][0],
                                          self.data_dict['dihedral'][1]):
                    self.edit_show3.insertPlainText(f'Atoms {atoms_} --')
                    self.edit_show3.insertPlainText(f'{label_}\n')
            self.ax3.clear()
            self.canvas3.draw()
            self.button_plot3.clicked.emit()

        elif button_clicked == 'button_remove4':
            self.data_dict['rmsd'][0].pop(-1)
            self.data_dict['rmsd'][1].pop(-1)
            self.edit_show4.clear()
            if len(self.data_dict['rmsd'][0]) != 0:
                for atoms_, label_ in zip(self.data_dict['rmsd'][0],
                                          self.data_dict['rmsd'][1]):
                    self.edit_show4.insertPlainText(f'Atoms {atoms_} --')
                    self.edit_show4.insertPlainText(f'{label_}\n')
            self.ax4.clear()
            self.canvas4.draw()
            self.button_plot4.clicked.emit()

        return 

    def plot_measure(self):
        button_clicked = self.sender().objectName()

        if button_clicked == 'button_plot1':
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

        elif button_clicked == 'button_plot2':
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
        
        elif button_clicked == 'button_plot3':
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
        
        elif button_clicked == 'button_plot4':
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
            self.edit_show1.clear()
        elif button_clicked == 'button_clear2':
            self.data_dict['angle'] = [[], []]
            canvas = self.canvas2
            ax=self.ax2
            self.edit_show2.clear()
        elif button_clicked == 'button_clear3':
            self.data_dict['dihedral'] = [[], []]
            canvas = self.canvas3
            ax=self.ax3
            self.edit_show3.clear()
        elif button_clicked == 'button_clear4':
            self.data_dict['rmsd'] = [[], []]
            canvas = self.canvas4
            ax=self.ax4
            self.edit_show4.clear()
        ax.clear()
        canvas.draw()

# Simula Window ================================================================
class SimulaWindow(QWidget):
    def __init__(self, traj):
        super().__init__()
        uifile = os.path.join(os.path.dirname(__file__),
                              'ui_files/simula_input.ui')
        uic.loadUi(uifile, self)   

        # Init variables:
        self.conditions_dict = {'type':[], 'atoms':[], 'symbol': [], 'value':[]}
        self.cv_dict = {'type':[], 'atoms':[], 'coeff':[], 'fpath':[],
                        'harm':[]}
        self.traj_file = ''
        self.top_file = ''
        self.frames_sel = None
        self.edit_mask = QLineEdit()

        # Combo boxes
        self.combo_cvtype.addItems(['DISTANCE', 'ANGLE', 'LCOD'])
        self.combo_theory.addItems(['DFTB3', 'EXTERN'])
        self.combo_filter.addItems(['distance', 'angle', 'dihedral'])
        self.combo_cond.addItems(['>', '<', '>=', '<='])
        self.combo_mask.addItems(['test', 'test2'])
        self.combo_mask.setLineEdit(self.edit_mask)
                
        # Connections
        self.button_add1.clicked.connect(self.add_cv)
        self.button_clear1.clicked.connect(self.clear)
        self.button_generate.clicked.connect(self.generate)
        self.checkBox_automatic.stateChanged.connect(self.onStateChanged)
        self.checkBox_manual.stateChanged.connect(self.onStateChanged)
        self.button_add2.clicked.connect(self.add_filter)
        self.button_clear2.clicked.connect(self.clear_filters)
        self.button_frames.clicked.connect(self.select_frames)
        self.button_originalTop.clicked.connect(self.load_file)
        self.button_originalTraj.clicked.connect(self.load_file)
        self.button_output.clicked.connect(self.sel_output)
        self.button_addmask.clicked.connect(self.add_mask)
        
        # Group check_boxes
        self.group_checked = QButtonGroup()
        self.group_checked.addButton(self.checkBox_automatic)
        self.group_checked.addButton(self.checkBox_manual)
        self.group_checked.setExclusive(True)


    def add_mask(self):
        print(self.combo_mask.currentText())

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
        return
    
    def add_filter(self):
        filter_type = self.combo_filter.currentText()
        atoms = self.edit_atoms2.text().strip().split()
        condition = self.combo_cond.currentText()
        condition_value = self.edit_cond.text()

        if atoms == '':
            fc.pop_error("Error!!!", "No atoms have been selected")
            return
        elif atoms != '':
            if filter_type == 'distance' and len(atoms) != 2:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            elif filter_type == 'angle' and len(atoms) != 3:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            elif filter_type == 'dihedral' and len(atoms) != 4:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return

        if condition_value == '':
            fc.pop_error("Error!!!", "No value for condition have been provided")
            return
        self.edit_atoms2.setText("")
        self.edit_cond.setText("")
        
        # Save conditions
        self.conditions_dict['type'].append(filter_type)
        self.conditions_dict['atoms'].append(atoms)
        self.conditions_dict['symbol'].append(condition)
        self.conditions_dict['value'].append(condition_value)

        self.edit_show2.insertPlainText(f"Filter: {filter_type}\n")
        self.edit_show2.insertPlainText(f"Atoms: {atoms}\n")
        self.edit_show2.insertPlainText(f"Condition: {condition}")
        self.edit_show2.insertPlainText(f"{condition_value}\n/\n")
        return
    
    def clear_filters(self):
        self.conditions_dict = {'type':[], 'atoms':[], 'symbol': [], 'value':[]}
        self.edit_show2.clear()
        return
    
    def select_frames(self):    
        self.traj = fc.load_traj(self.traj_file, self.top_file)
        if not self.traj:
            fc.pop_error("Error!!!", "Load the topology and trajectory files")
            return
        
        if self.checkBox_manual.isChecked():
            if self.edit_frames.text() == '':
                fc.pop_error("Error!!!", "No frames were selected")
                return
            self.frames_sel = self.edit_frames.text().strip().split(' ')
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
                
        self.edit_show2.insertPlainText(f"Total of {len(self.frames_sel)}")
        self.edit_show2.insertPlainText(f" frames selected")
        return
    
    def add_cv(self):
        atoms = self.edit_atoms.text().strip().split()
        try:
            atoms = np.asarray(atoms, dtype=int)
        except:
            fc.pop_error("Error!!!", "Could not read atoms")
        try:
            fpath = float(self.edit_path.text().strip())
        except:
            fc.pop_error("Error!!!", "Could not read final path")
            return
        cv_type = self.combo_cvtype.currentText()
        coeff_flag = False
      
        # Type of CV
        if cv_type == 'DISTANCE':
            if len(atoms) != 2:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            else:
                self.cv_dict['coeff'].append(None)
                    
        elif cv_type == 'ANGLE':
            if len(atoms) != 3:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            else:
                self.cv_dict['coeff'].append(None)

        elif cv_type == 'LCOD':
            if len(atoms)%2 != 0:
                fc.pop_error("Error!!!", "Wrong number of atoms")
                return
            else:
                try:
                    coeff = self.edit_coeff.text().strip().split(' ')
                    coeff = np.asarray(coeff, dtype=float)
                    coeff_flag = True
                except ValueError:
                    fc.pop_error("Error!!!", "Coefficients must be supplied")
                    return
                if len(coeff) != len(atoms)//2:
                    fc.pop_error("Error!!!", "Wrong number of coefficients")
                    return
                self.cv_dict['coeff'].append(coeff)

        self.cv_dict['type'].append(cv_type)
        self.cv_dict['atoms'].append(atoms)
        self.cv_dict['fpath'].append(fpath)

        if self.edit_harm.text() != '':
            harm = float(self.edit_harm.text().strip())
        else:
            harm = 600
        self.cv_dict['harm'].append(harm)

        # Print collective variables
        self.edit_show1.insertPlainText(f'type = {cv_type}\n')
        self.edit_show1.insertPlainText(f'atoms = {atoms}\n')
        self.edit_show1.insertPlainText(f'final path = {fpath}\n')
        if coeff_flag:
            self.edit_show1.insertPlainText(f'Coefficients = {coeff}\n')
        self.edit_show1.insertPlainText(f'Harmonic = {harm}\n/\n')
        
        # Clean edits
        self.edit_atoms.setText("")
        self.edit_path.setText("")
        self.edit_coeff.setText("")
        self.edit_harm.setText("")

    def load_file(self):
        '''
        Load topology and trajectory file
        '''
        button_clicked = self.sender().objectName()
        if button_clicked == 'button_originalTraj':
            self.traj_file = fc.get_traj_top('traj', self.edit_trajectory)
        elif button_clicked == 'button_originalTop':
            self.top_file = fc.get_traj_top('top', self.edit_topology)
        
    def clear(self):
        self.edit_show1.clear()
        
    def sel_output(self):
        self.output_dir = str(QFileDialog.getExistingDirectory(self,
                                                               'Select output dir'))
        self.edit_output.setText(self.output_dir)

    def generate(self):
        # Get qmmm.in params
        if self.edit_mask == '':
            fc.pop_error("Error!!!", "1")
            return
        elif self.edit_charge == '':
            fc.pop_error('Error!!!', "2")
            return
        elif self.edit_steps == '':
            fc.pop_error('Error!!!', "3")
            return
        elif self.edit_timestep == '':
            fc.pop_error('Error!!!', "4")
            return
        
        # qmmm_dict = {'mask': self.edit_mask.text().strip(),
        #              'theory': self.combo_theory.currentText(),
        #              'charge': int(self.edit_charge.text().strip()),
        #              'steps': int(self.edit_steps.text().strip()),
        #              'timestep': float(self.edit_timestep.text().strip())}
        
        # Write cv.in file
        # for frame_id, frame in enumerate(self.frames_sel):
        #     os.mkdir(f'{self.output_dir}/qmmm_{frame}')
        #     outdir = f'{self.output_dir}/qmmm_{frame}'
        #     fc.write_cv(self.traj, frame_id, self.cv_dict, outdir)
        #     fc.write_qmmm(qmmm_dict, outdir)
        # return
    

        




if __name__ == '__main__':
    app = QApplication([])
    window = MainWindow()
    app.exec_()
