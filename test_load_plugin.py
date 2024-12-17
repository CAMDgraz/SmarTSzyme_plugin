from PyQt5 import uic
from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QFileDialog,
                             QMessageBox)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure

import functions as fc
import seaborn as sns
import mdtraj as md
import numpy as np
import math
import os

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__() # Call the inherited classes __init__ method
        uifile = os.path.join(os.path.dirname(__file__), 'mainWindow.ui')
        uic.loadUi(uifile, self) # Load the .ui file
        self.show() # Show the GUI
        self.traj = None
        self.traj_file = ''
        self.top_file = ''
    
        # Connections actions
        self.actionMeasures.triggered.connect(self.show_Measures)
        self.actionSimula_Input.triggered.connect(self.show_Simula)

        # Connections buttons
        self.mdtraj_button.clicked.connect(self.loadMDTraj)
        self.button_trajectory.clicked.connect(self.get_traj)
        self.button_topology.clicked.connect(self.get_top)

    def show_Measures(self):
        if not self.traj:
            fc.pop_error('Error!!!', 'No system loaded')
            pass
        elif self.traj:
            self.measureWindow = MeasureWindow(self.traj)
            self.measureWindow.show()
    
    def show_Simula(self):
        if not self.traj:
            fc.pop_error('Error!!!', 'No system loaded')
            pass
        elif self.traj:
            self.simulaWindow = SimulaWindow(self.traj)
            self.simulaWindow.show()
       
    def get_traj(self):
        initial_dir = './'
        ext_filter = ''

        for ext in fc.valid_trajs:
            ext_filter += f'*.{ext};;'
        ext_filter += 'All files (*)'

        # Get path or open file dialog
        if self.edit_trajectory.text() == "":
            traj_file, _ = QFileDialog.getOpenFileName(self,
                                                       "Open Trajectory File",
                                                       initial_dir,
                                                       ext_filter)
            self.edit_trajectory.setText(traj_file)
        else:
            traj_file = self.edit_trajectory.text()
        
        # Check valid trajectry extension
        if fc.is_valid_traj(traj_file):
            self.traj_file = traj_file
        else:
            self.edit_trajectory.setText('')

    def get_top(self):
        initial_dir = './'
        ext_filter = ''

        for ext in fc.valid_tops:
            ext_filter += f'*.{ext};;'
        ext_filter += 'All files (*)'

        # Get path or open file dialog
        if self.edit_topology.text() == "":
            top_file, _ = QFileDialog.getOpenFileName(self,
                                                       "Open Topology File",
                                                       initial_dir,
                                                       ext_filter)
            self.edit_topology.setText(top_file)
        else:
            top_file = self.edit_topology.text()
        
        # Check valid topology extension
        if fc.is_valid_top(top_file):
            self.top_file = top_file
        else:
            self.edit_topology.setText('')

    def loadMDTraj(self):
        self.traj = fc.load_traj(self.traj_file, self.top_file)
        if self.traj:
            self.label_system.setText(str(self.traj))


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)

class MeasureWindow(QWidget):
    def __init__(self, traj):
        super().__init__()
        uifile = os.path.join(os.path.dirname(__file__), 'measures.ui')
        uic.loadUi(uifile, self) # Load the .ui file
        self.traj = traj

        # Plot canvas Distance
        self.canvas1 = MplCanvas(self, width=5, height=4, dpi=100)
        toolbar1 = NavigationToolbar2QT(self.canvas1, self)
        self.verticalLayout_1.addWidget(toolbar1)
        self.verticalLayout_1.addWidget(self.canvas1)
        self.ax1 = self.canvas1.fig.add_subplot(111)

        # Plot canvas Angle
        self.canvas2 = MplCanvas(self, width=5, height=4, dpi=100)
        toolbar2 = NavigationToolbar2QT(self.canvas2, self)
        self.verticalLayout_2.addWidget(toolbar2)
        self.verticalLayout_2.addWidget(self.canvas2)
        self.ax2 = self.canvas2.fig.add_subplot(111)

        # Plot canvas Dihedral
        self.canvas3 = MplCanvas(self, width=5, height=4, dpi=100)
        toolbar3 = NavigationToolbar2QT(self.canvas3, self)
        self.verticalLayout_3.addWidget(toolbar3)
        self.verticalLayout_3.addWidget(self.canvas3)
        self.ax3 = self.canvas3.fig.add_subplot(111)

        # Plot canvas RMSD
        self.canvas4 = MplCanvas(self, width=5, height=4, dpi=100)
        toolbar4 = NavigationToolbar2QT(self.canvas4, self)
        self.verticalLayout_4.addWidget(toolbar4)
        self.verticalLayout_4.addWidget(self.canvas4)
        self.ax4 = self.canvas4.fig.add_subplot(111)

        # Variables 1-> distance 2-> angle 3-> dihedral 4-> rmsd
        self.atoms1 = []
        self.labels1 = []
        self.atoms2 = []
        self.labels2 = []
        self.atoms3 = []
        self.labels3 = []
        self.atoms4 = []
        self.labels4 = []

        # Connections
        self.button_add1.clicked.connect(lambda: self.add(1))
        self.button_add2.clicked.connect(lambda: self.add(2))
        self.button_add3.clicked.connect(lambda: self.add(3))
        self.button_add4.clicked.connect(lambda: self.add(4))

        self.button_clear1.clicked.connect(lambda: self.clear(1))
        self.button_clear2.clicked.connect(lambda: self.clear(2))
        self.button_clear3.clicked.connect(lambda: self.clear(3))
        self.button_clear4.clicked.connect(lambda: self.clear(4))

        self.button_plot1.clicked.connect(lambda: fc.plot_measure(1,
                                                                   self.traj,
                                                                   self.atoms1,
                                                                   self.labels1,
                                                                   self.canvas1,
                                                                   self.ax1))
        self.button_plot2.clicked.connect(lambda: fc.plot_measure(2,
                                                                   self.traj,
                                                                   self.atoms2,
                                                                   self.labels2,
                                                                   self.canvas2,
                                                                   self.ax2))
        self.button_plot3.clicked.connect(lambda: fc.plot_measure(3,
                                                                   self.traj,
                                                                   self.atoms3,
                                                                   self.labels3,
                                                                   self.canvas3,
                                                                   self.ax3))
        self.button_plot4.clicked.connect(lambda: fc.plot_measure(4,
                                                                   self.traj,
                                                                   self.atoms4,
                                                                   self.labels4,
                                                                   self.canvas4,
                                                                   self.ax4))
        self.button_csv1.clicked.connect(lambda: fc.to_csv(1, self.traj,
                                                           self.atoms1,
                                                           self.labels1,
                                                           'distance.csv'))
        self.button_csv2.clicked.connect(lambda: fc.to_csv(2, self.traj,
                                                           self.atoms2,
                                                           self.labels2,
                                                           'angle.csv'))
        self.button_csv3.clicked.connect(lambda: fc.to_csv(3, self.traj,
                                                           self.atoms3,
                                                           self.labels3,
                                                           'dihedral.csv'))
        self.button_csv4.clicked.connect(lambda: fc.to_csv(4, self.traj,
                                                           self.atoms4,
                                                           self.labels4,
                                                           'rmsd.csv'))
        
    def add(self, type):
        # Distances
        if type == 1:
            selection = self.edit_atom1.text()
            sel_idx = fc.atoms_sel(self.traj, selection)
            if sel_idx.any():
                if len(sel_idx) != 2:
                    fc.pop_error('Selection Error',
                                 f'The number of atoms is incorrect')
            label = self.edit_label1.text()
            if label != '':
                self.edit_show1.insertPlainText(f'Atoms {sel_idx} -- {label}\n')
            else:
                label = f'Atoms {sel_idx}'
                self.edit_show1.insertPlainText(f'Atoms {sel_idx} -- {label}\n')

            self.atoms1.append(sel_idx)
            self.labels1.append(label)
        
            self.edit_atom1.setText('')
            self.edit_label1.setText('')

        # Angles
        if type == 2:
            selection = self.edit_atom2.text()
            sel_idx = fc.atoms_sel(self.traj, selection)
            if sel_idx.any():
                if len(sel_idx) != 3:
                    fc.pop_error('Selection Error',
                                 f'The number of atoms is incorrect')
            label = self.edit_label2.text()
            if label != '':
                self.edit_show2.insertPlainText(f'Atoms {sel_idx} -- {label}\n')
            else:
                label = f'Atoms {sel_idx}'
                self.edit_show2.insertPlainText(f'Atoms {sel_idx} -- {label}\n')

            self.atoms2.append(sel_idx)
            self.labels2.append(label)
        
            self.edit_atom2.setText('')
            self.edit_label2.setText('') 

        # Dihedrals
        if type == 3:
            selection = self.edit_atom3.text()
            sel_idx = fc.atoms_sel(self.traj, selection)
            if sel_idx.any():
                if len(sel_idx) != 4:
                    fc.pop_error('Selection Error',
                                 f'The number of atoms is incorrect')
            label = self.edit_label3.text()
            if label != '':
                self.edit_show3.insertPlainText(f'Atoms {sel_idx} -- {label}\n')
            else:
                label = f'Atoms {sel_idx}'
                self.edit_show3.insertPlainText(f'Atoms {sel_idx} -- {label}\n')

            self.atoms3.append(sel_idx)
            self.labels3.append(label)
        
            self.edit_atom3.setText('')
            self.edit_label3.setText('')

        # RMSD
        if type == 4:
            selection = self.edit_atom4.text()
            sel_idx = fc.atoms_sel(self.traj, selection)
            if sel_idx.any():
                pass
            label = self.edit_label4.text()
            if label != '':
                self.edit_show4.insertPlainText(f'Atoms {sel_idx} -- {label}\n')
            else:
                label = f'Atoms {sel_idx}'
                self.edit_show4.insertPlainText(f'Atoms {sel_idx} -- {label}\n')

            self.atoms4.append(sel_idx)
            self.labels4.append(label)
        
            self.edit_atom4.setText('')
            self.edit_label4.setText('') 
        

    def clear(self, type):
        if type == 1:
            self.atoms1 = []
            self.labels1 = []
            self.edit_show1.clear()
        elif type == 2:
            self.atoms2 = []
            self.labels2 = []
            self.edit_show2.clear()
        elif type == 3:
            self.atoms3 = []
            self.labels3 = []
            self.edit_show3.clear()
        elif type == 4:
            self.atoms4 = []
            self.labels4 = []
            self.edit_show4.clear()


class SimulaWindow(QWidget):
    def __init__(self, traj):
        super().__init__() # Call the inherited classes __init__ method
        uifile = os.path.join(os.path.dirname(__file__), 'simula_input.ui')
        uic.loadUi(uifile, self) # Load the .ui file   
        self.traj = traj

        # Combo boxes
        self.combo_cvtype.addItems(['DISTANCE', 'ANGLE', 'LCOD'])
        self.combo_frames.addItems(['distance', 'angle', 'dihedral'])
        self.combo_theory.addItems(['DFTB3', 'EXTERN'])

if __name__ == '__main__':
    app = QApplication([])
    window = MainWindow()
    app.exec_()
