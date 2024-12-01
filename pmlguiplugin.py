"""
Created on Mon Nov 18 09:39:43 2024

@author: aliaa
"""

from PyQt5.QtCore import QRegExp
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtWidgets import (QMainWindow, QApplication, QPushButton, QWidget, 
                             QTabWidget, QVBoxLayout, QLineEdit, QFileDialog,
                             QLabel, QHBoxLayout, QMessageBox)
from . import vispml
from . import md_analysis as mda
from pymol import cmd
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import mdtraj as md

#_______________________________Main Window____________________________________

class Plugin(QMainWindow):
    
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("PyMOL GUI Plugin")
        self.setGeometry(0, 0, 448, 350)
        self.setCentralWidget(TabsWidget(self))
        
#_________________________________Tabs Widget__________________________________
        
class TabsWidget(QWidget):
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # initialize loadTraj in tabswidget and dir
        self.loadTraj = None
        self.initial_dir = "/home/user/Desktop"
        
        # Layout for the TabWidget
        self.tabs_layout = QVBoxLayout(self)
        
        # Initialize tab screen
        self.tabs = QTabWidget(self)
        self.tab1 = QWidget(self)
        self.tab2 = QWidget(self)
        self.tab3 = QWidget(self)
        
        # Add tabs
        self.tabs.addTab(self.tab1, "Input Analysis")
        self.tabs.addTab(self.tab2, "Simula Settings")
        self.tabs.addTab(self.tab3, "Visualization")
        
#________________________________Tab 1_________________________________________
       
        # main layout for tab 1
        self.mainlayout_1 = QVBoxLayout()
        self.tab1.setLayout(self.mainlayout_1)
        
        # traj file button
        self.trajFileButton = QPushButton("Choose Trajectory File", self)
        self.mainlayout_1.addWidget(self.trajFileButton) 
        # Label for Trajectory File (added directly below the button)
        self.labelTrajSelected = QLabel("No Trajectory file selected", self)
        self.mainlayout_1.addWidget(self.labelTrajSelected)
        self.trajFileButton.clicked.connect(self.open_traj_file_dialog)
        
        # topo file button
        self.topFileButton = QPushButton("Choose Topology File", self)
        self.mainlayout_1.addWidget(self.topFileButton)
        self.labelTopoSelected = QLabel("No Topology file selected", self)
        self.mainlayout_1.addWidget(self.labelTopoSelected)
        self.topFileButton.clicked.connect(self.open_topo_file_dialog)

#..............................................................................
        
        # button for loading files in mdtraj
        self.load_layout = QHBoxLayout()
        self.loadInMdtrajButton = QPushButton("Load in MDTraj", self)
        self.mainlayout_1.addWidget(self.loadInMdtrajButton)
        self.loadInMdtrajButton.clicked.connect(self.load_in_mdtraj)
        
        # button for loading files in pymol 
        self.loadInPymolButton = QPushButton("Load in PyMOL", self)
        self.mainlayout_1.addWidget(self.loadInPymolButton)
        self.loadInPymolButton.clicked.connect(self.load_in_pymol)

        # add space 
        self.mainlayout_1.addSpacing(10)  
        
#..............................................................................
        
        # main layout and horizontal        
        self.horizontal_layout = QHBoxLayout()
        self.mainlayout_1.addLayout(self.horizontal_layout)
        
        # measure button
        self.left_layout = QVBoxLayout()
        self.horizontal_layout.addLayout(self.left_layout)
        self.measureButton = QPushButton("Measure", self)
        self.left_layout.addWidget(self.measureButton)
        self.measureButton.clicked.connect(self.open_new_window_plot)

        # rmsd button
        self.right_layout = QVBoxLayout()
        self.rmsdButton = QPushButton("RMSD", self)
        self.right_layout.addWidget(self.rmsdButton)        
        self.horizontal_layout.addLayout(self.right_layout)
        self.rmsdButton.clicked.connect(self.open_new_window_rmsd)
                
#________________________________Tab 2_________________________________________
        
        # main layout for tab 2
        self.mainlayout_2 = QVBoxLayout()        
        self.tab2.setLayout(self.mainlayout_2)
        
#______________________________________________________________________________

        # main layout for tab 3
        self.mainlayout_3 = QVBoxLayout()
        self.tab3.setLayout(self.mainlayout_3)

        # pdb file button
        self.choosePdbButton = QPushButton("Choose PDB File", self)
        self.mainlayout_3.addWidget(self.choosePdbButton)
        # Label for PDB File (added directly below the button)
        self.labelPdbSelected = QLabel("No *.pdb file selected", self)
        self.mainlayout_3.addWidget(self.labelPdbSelected)
        self.choosePdbButton.clicked.connect(self.open_pdb_file_dialog)
        
        # csv file button
        self.chooseCsvButton = QPushButton("Choose CSV File", self)
        self.mainlayout_3.addWidget(self.chooseCsvButton)        
        # Label for CSV File (added directly below the button)
        self.labelCsvSelected = QLabel("No *.csv file selected", self)
        self.mainlayout_3.addWidget(self.labelCsvSelected)
        self.chooseCsvButton.clicked.connect(self.open_csv_file_dialog)        
        
        # execute button
        self.executeButton = QPushButton("Execute", self)
        self.mainlayout_3.addWidget(self.executeButton)
        self.executeButton.clicked.connect(self.execute_vispml)          
        
#..............................................................................
        
        # Add tabs to widget
        self.tabs_layout.addWidget(self.tabs)
        self.setLayout(self.tabs_layout)
        
#____________________________Methods for Tab 1_________________________________

    def open_traj_file_dialog(self):
        # Filter to allow only Trajectory Files
        file_filter = "Trajectory File (*.arc; *.dcd; *.binpos; *.xtc; *.trr; *.hdf5; *.h5; *.ncdf; *.netcdf; *.nc; *.pdb.gz; *.pdb; *.lh5; *.crd; *.mdcrd; *.inpcrd; *.restrt; *.rst7; *.ncrst; *.lammpstrj; *.dtr; *.stk; *.gro; *.xyz.gz; *.xyz; *.tng; *.xml; *.mol2; *.hoomdxml; *.gsd)"  # Update this filter as needed

        # Open the file dialog
        self.traj_file, _ = QFileDialog.getOpenFileName(self, "Open Trajectory File", self.initial_dir, file_filter)
        
        if mda.is_valid_traj(self.traj_file, mda.valid_trajs) == True and self.traj_file:
            self.labelTrajSelected.setText(f"Selected Trajectory file: {self.traj_file}")
        else:
            self.labelTrajSelected.setText("No Trajectory file selected")
            
        self.needForTop_Tab1 = mda.need_for_top(self.traj_file)
        
        if self.needForTop_Tab1 == False:
            self.topFileButton.setEnabled(False)   
            self.topFileButton.setText("No Seperate Topology File needed")
        # TODO: fix when pdb is selected and then a traj that needs a top is selected button is still greyed out    

#..............................................................................

    def open_topo_file_dialog(self):
        # Filter to allow only Topology Files (e.g., PDB files)
        file_filter = "Topology Files (*.pdb; *.pdb.g; *.h5; *.lh5; *.prmtop; *.parm7; *.prm7; *.psf; *.mol2; *.hoomdxml; *.gro; *.arc; *.hdf5; *.gsd)"  # Update this filter as needed
        
        # Open the file dialog
        self.top_file, _ = QFileDialog.getOpenFileName(self, "Open Topology File", self.initial_dir, file_filter)
    
        if self.top_file:
            # Update the label to display the selected file path
            self.labelTopoSelected.setText(f"Selected Topology file: {self.top_file}")
        else:
            self.labelTopoSelected.setText("No Topology file selected")

#..............................................................................            

    def load_in_pymol(self):
         
        cmd.load(self.traj_file)
        
        if mda.need_for_top(self.traj_file) == True:
            cmd.load(self.top_file)
            cmd.load_traj(self.traj_file)         
            
#..............................................................................
            
    def load_in_mdtraj(self):
        
        if mda.need_for_top(self.traj_file) == True:
            self.loadTraj = mda.load_traj(self.traj_file, self.top_file, mda.valid_trajs, mda.valid_tops)
        else:      
            self.loadTraj = mda.load_traj(self.traj_file, self.traj_file, mda.valid_trajs, mda.valid_tops)
        
        print(self.loadTraj) 

#..............................................................................
        
    def open_new_window_rmsd(self):
        # Create an instance of the new window and show it
        self.new_window = RMSDWindow(self.loadTraj)
        self.new_window.show()
        
#..............................................................................
        
    def open_new_window_plot(self):
        # Create an instance of the new window and show it
        self.new_window_plot = MeasuresWindow(self.loadTraj)
        self.new_window_plot.show()
            
#____________________________Methods for Tab 2_________________________________
    


#____________________________Methods for Tab 3_________________________________

    def open_pdb_file_dialog(self):
        # Filter to allow only PDB Files
        file_filter = "PDB Files (*.pdb)"
   
        # Open the file dialog
        self.file_pdb, _ = QFileDialog.getOpenFileName(self, "Choose PDB File", self.initial_dir, file_filter)
   
        if self.file_pdb:
            self.labelPdbSelected.setText(f"{self.file_pdb}")
        else:
            self.labelPdbSelected.setText("No PDB file selected")
         
#..............................................................................
         
    def open_csv_file_dialog(self):
        # filter to allow only CSV File
        file_filter = "CSV File (*.csv)"
   
        # open the file dialog
        self.file_csv, _ = QFileDialog.getOpenFileName(self, "Choose CSV File", self.initial_dir, file_filter)
        
        if self.file_csv:
            self.labelCsvSelected.setText(f"{self.file_csv}")
        else:
            self.labelCsvSelected.setText("No CSV file selected")
     
#..............................................................................
     
    def execute_vispml(self):

        # Generate results in the current directory
        self.visualize = vispml.results_pdb(vispml.read_pdb(self.file_pdb),vispml.normalize_flux(self.file_csv))
    
        # Load the generated PDB file in PyMOL
        cmd.load("results.pdb")
    
        # Apply visualization settings directly in PyMOL
        cmd.show("cartoon")
        cmd.bg_color("white")
        cmd.spectrum("b", "rainbow")
    
        # Select and color atoms with B-factor == 0
        cmd.select("bfactor_zero", "b < 0.0001")  # Select atoms with zero B-factor
        cmd.color("black", "bfactor_zero")
    
        # Color substrate (non-polymer residues) gray
        cmd.select("substrate", "not polymer")
        cmd.color("grey", "substrate")

#_____________________________Measures Window__________________________________    

class MeasuresWindow(QWidget):
    def __init__(self, loadTraj, parent=None):
        super().__init__(parent)
        
        # initialize variables
        self.atom_list = []
        self.loadTraj = loadTraj
        
        # main layout for measures window
        self.measures_layout = QHBoxLayout()
        self.setLayout(self.measures_layout)
        
        self.setWindowTitle("Measures Calculation")
        self.setGeometry(500,150,1000,700)
        
        # initialize plot canvas and toolbar
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
#..............................................................................        
        
        # layout for the left side of measures window
        self.left_layout = QVBoxLayout()
        self.measures_layout.addLayout(self.left_layout)
        
        # toolbar layout
        self.left_first_layout = QHBoxLayout()
        self.left_layout.addLayout(self.left_first_layout)
        self.left_first_layout.addWidget(self.toolbar)

        # layout for plot
        self.left_second_layout = QHBoxLayout()
        self.left_layout.addLayout(self.left_second_layout)
        self.left_second_layout.addWidget(self.canvas)
        
#..............................................................................        
        
        # layout for the left side of measures window
        self.right_layout = QVBoxLayout()
        self.measures_layout.addLayout(self.right_layout)
        
        # layout for first line of right side of measures window
        self.right_first_layout = QHBoxLayout()
        self.right_layout.addLayout(self.right_first_layout)
        
        # select atoms in text field
        self.select_atoms = QLineEdit(self)
        self.select_atoms.setPlaceholderText("Select Atoms")
        self.right_first_layout.addWidget(self.select_atoms)
        
        # define the regex pattern for frames and atoms
        pattern = r"^\d+(,\d+)*$"
        regex = QRegExp(pattern)
        # apply QRegExpValidator for both inputs (frames and atoms)
        validator_select_atoms = QRegExpValidator(regex, self.select_atoms)
        self.select_atoms.setValidator(validator_select_atoms)

        # add button to save atom selection
        self.addButton = QPushButton("Add", self)
        self.right_first_layout.addWidget(self.addButton)
        self.addButton.clicked.connect(self.add_atoms)

        # clear button to delete atoms selection
        self.clearButton = QPushButton("Clear", self)
        self.right_first_layout.addWidget(self.clearButton)
        self.clearButton.clicked.connect(self.clear_atoms)
        
#..............................................................................                
        
        # layout for first line of right side of  measures window
        self.right_second_layout = QHBoxLayout()
        self.right_layout.addLayout(self.right_second_layout)
        
        # plot button to plot atom selection
        self.plotButton = QPushButton("Plot", self)
        self.right_second_layout.addWidget(self.plotButton)
        self.plotButton.clicked.connect(self.plot_measures)

        # clear button to delete plot
        self.clearplotButton = QPushButton("Clear Plot", self)
        self.right_second_layout.addWidget(self.clearplotButton)
        self.clearplotButton.clicked.connect(self.clear_plot)

#_________________________Methods for Measures Window__________________________

    def add_atoms(self):
        atom_text = self.select_atoms.text()
        splitted = atom_text.split(',')
        
        if len(splitted) > 4 or len(splitted) < 2:
            self.clear_atoms()
            top_popup = QMessageBox(self)
            top_popup.setWindowTitle("Error")
            top_popup.setText("Wrong number of atoms!")
            top_popup.setIcon(QMessageBox.Information)
            top_popup.setStandardButtons(QMessageBox.Ok)
            top_popup.exec_()
        else:
            self.atom_list.append(splitted)
        print(self.atom_list)

#..............................................................................                

    def clear_atoms(self):
        self.atom_list = []
        print("Atoms list has been cleared.")

#..............................................................................                

    def plot_measures(self):
        """plot based on the atom list"""
        
        if not self.atom_list:
            print("No data to plot.")
            return
        
        self.figure.clear()
        ax = self.figure.add_subplot()
        
        # Table with top info
        table, _ = self.loadTraj.topology.to_dataframe()
        
        for atoms in self.atom_list:
            atoms = np.asarray(atoms, dtype=int)
            atoms -= 1
            
            label = ''
            for atom in atoms:
                label += (f':{table.iloc[atom, 3]}@{table.iloc[atom, 1]} ')
            # Plot dihedral
            if len(atoms) == 4:
                atoms = np.reshape(atoms, (-1, 4))
                dihedral = md.compute_dihedrals(self.loadTraj, atoms)
                sns.kdeplot(x=dihedral.T[0], ax=ax, label=label)
                ax.set_title("Density Plot")
                ax.set_xlabel("Dihedral")
                ax.set_ylabel("Density")
            # Plot angle
            elif len(atoms) == 3:
                atoms = np.reshape(atoms, (-1, 3))
                angle = md.compute_angles(self.loadTraj, atoms)
                sns.kdeplot(x=angle.T[0], ax=ax, label=label)
                ax.set_title("Density Plot")
                ax.set_xlabel("Angle")
                ax.set_ylabel("Density")
            # Plot distance
            elif len(atoms) == 2:
                atoms = np.reshape(atoms, (-1, 2))
                distance = md.compute_distances(self.loadTraj, atoms)
                sns.kdeplot(x=distance.T[0], ax=ax, label=label)
                ax.set_title("Density Plot")
                ax.set_xlabel("Distance")
                ax.set_ylabel("Density")
                
        ax.legend()  # Add legend to differentiate between plots
        self.canvas.draw()  
 
#..............................................................................                        
        
    def clear_plot(self):
        self.figure.clear()
        ax = self.figure.add_subplot()
        ax.set_title("Density Plot")
        ax.set_xlabel("x")
        ax.set_ylabel("Density") 
                      
        print("Measures plot has been cleared.")
        
        self.canvas.draw()  
        
#_________________________________RMSD Window__________________________________
        
class RMSDWindow(QWidget):
    def __init__(self, loadTraj, parent=None):
        super().__init__(parent)
        
        # initialize loadTraj and lists
        self.loadTraj = loadTraj
        self.selection_list = []
        self.labels = []
        
        # layout for rmsd window
        self.rmsd_layout = QHBoxLayout()
        self.setLayout(self.rmsd_layout)
            
        self.setWindowTitle("RMSD Calculation")
        self.setGeometry(500,150,1000,700)
        
        # matplotlib figure and canvas instance to plot on
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        # layout for left side of rmsd window
        self.left_layout = QVBoxLayout()
        self.rmsd_layout.addLayout(self.left_layout)
        
        # layout for first line of left side of rmsd window
        self.left_first_layout = QHBoxLayout()
        self.left_layout.addLayout(self.left_first_layout)
        
        # layout for second line of right side of rmsd window
        self.left_second_layout = QHBoxLayout()
        self.left_layout.addLayout(self.left_second_layout)
        self.left_second_layout.addWidget(self.canvas)
        self.left_first_layout.addWidget(self.toolbar)
        
        # layout for right side of rmsd window
        self.right_layout = QVBoxLayout()
        self.rmsd_layout.addLayout(self.right_layout)
        
        # layout for first line of right side of rmsd window
        self.right_first_layout = QHBoxLayout()
        self.right_layout.addLayout(self.right_first_layout)
        
        self.select_selection = QLineEdit(self)
        self.select_selection.setPlaceholderText("Selection")
        
        self.right_first_layout.addWidget(self.select_selection)
        self.addButton_selection = QPushButton("Add", self)
        self.right_first_layout.addWidget(self.addButton_selection)
        self.addButton_selection.clicked.connect(self.add_selection)

        self.clearButton_selection = QPushButton("Clear", self)
        self.right_first_layout.addWidget(self.clearButton_selection)
        self.clearButton_selection.clicked.connect(self.clear_selection)
        
        # layout for second line of right side of rmsd window
        self.right_second_layout = QHBoxLayout()
        self.right_layout.addLayout(self.right_second_layout)

        self.select_frames = QLineEdit(self)
        self.select_frames.setPlaceholderText("Select Frames: first:last:stride")
        self.right_second_layout.addWidget(self.select_frames)
        
        self.right_second_layout.addWidget(self.select_frames)
        self.addButton_frames = QPushButton("Add", self)
        self.right_second_layout.addWidget(self.addButton_frames)
        self.addButton_frames.clicked.connect(self.select_frames_function)

        # layout for third line of right side of rmsd window
        self.right_third_layout = QHBoxLayout()
        self.right_layout.addLayout(self.right_third_layout)
        
        self.plotButton = QPushButton("Plot", self)
        self.right_third_layout.addWidget(self.plotButton)
        self.plotButton.clicked.connect(self.plot_rmsd)
        
        self.clearplotButton = QPushButton("Clear Plot", self)
        self.right_third_layout.addWidget(self.clearplotButton)
        self.clearplotButton.clicked.connect(self.clear_plot)
        
#___________________________Methods for RMSD Window____________________________        
    
    def add_selection(self):
        selection_text = self.select_selection.text()
        sel_idx = mda.atoms_sel(self.loadTraj, selection_text)
        if sel_idx.any():
            self.selection_list.append(sel_idx)
            self.labels.append(f'{self.select_selection.text()}')
            print(self.selection_list)
            # print(self.labels)
        
#..............................................................................                        

    def select_frames_function(self):
        self.first = 0 #index of frame not frame
        self.last = None
        self.stride = 1
        
        self.frames_text = self.select_frames.text().split(':')
        if self.frames_text[0] != '':
            self.first = int(self.frames_text[0]) - 1
        if self.frames_text[1] != '':
            self.last = int(self.frames_text[1]) - 1
        if self.frames_text[2] != '':
            self.stride = self.frames_text[2]

        self.traj = mda.range_traj(self.loadTraj, self.first, self.last,
                                   self.stride)
    
#..............................................................................                        

    def plot_rmsd(self):
        """Plot rmsd based on the selection and the frames."""
        
        if not self.selection_list:
            print("No data to plot.")
            return
        
        self.figure.clear()
        ax = self.figure.add_subplot()
        ax.set_ylabel(r'RMSD ($\AA$)')
        ax.set_xlabel(r'Frames')
        ax.set_title(r'RMSD')
        
        if not self.last:
            x_axis = np.arange(self.first + 1, self.loadTraj.n_frames + 1,
                               self.stride)
        else:
            x_axis = np.arange(self.first + 1, self.last+1, self.stride)
    
        rmsd_per_sel = np.zeros((len(self.selection_list), self.traj.n_frames))
        
        for idx, selection in enumerate(self.selection_list):
            rmsd = md.rmsd(self.traj, self.traj, frame=0, atom_indices=selection)
            rmsd *= 10
            rmsd_per_sel[idx] = rmsd
            
            ax.plot(x_axis, rmsd, label = self.labels[idx])
            
            ax.legend()
        self.canvas.draw() 
        
#..............................................................................                        

    def clear_selection(self):
        self.selection_list = []
        print("selection list cleared.")
        
#..............................................................................                        

    def clear_plot(self):
        self.figure.clear()
        ax = self.figure.add_subplot()
        ax.set_ylabel(r'RMSD ($\AA$)')
        ax.set_xlabel(r'Frames')
        ax.set_title(r'RMSD')
        
        print("RMSD plot has been cleared.")
        self.canvas.draw()  
        
#______________________________________________________________________________

if __name__ == '__main__':
    app = QApplication([])
    p = Plugin()
    p.show()
    app.exec_()
