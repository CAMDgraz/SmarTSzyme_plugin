# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 09:39:43 2024

@author: aliaa
"""

from PyQt5.QtCore import QRegExp
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QTabWidget, QVBoxLayout, QLineEdit, QFileDialog, QLabel, QHBoxLayout
import vispml
import subprocess



initial_dir = "/home/user/Documents"

class Plugin(QMainWindow):
    
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("PyMOL GUI Plugin")
        self.setGeometry(0, 0, 448, 350)

        self.setCentralWidget(TabsWidget(self))
#______________________________________________________________________________

class TabsWidget(QWidget):
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Layout for the TabWidget
        self.layout = QVBoxLayout(self)
        
        # Initialize tab screen
        self.tabs = QTabWidget(self)
        self.tab1 = QWidget(self)
        self.tab2 = QWidget(self)
        self.tab3 = QWidget(self)
        
        # Add tabs
        self.tabs.addTab(self.tab1, "Molecular Dynamics")
        self.tabs.addTab(self.tab2, "Quantum Mechanics")
        self.tabs.addTab(self.tab3, "Visualization")
#______________________________________________________________________________     
       
        # Create Tab 1
        self.tab1.layout = QVBoxLayout()

        # topo file button
        self.pushButton1 = QPushButton("Choose Topology File", self)
        self.tab1.layout.addWidget(self.pushButton1)
        
        # Label for Topology File (added directly below the button)
        self.labelTopoSelected = QLabel("No Topology file selected", self)
        self.tab1.layout.addWidget(self.labelTopoSelected)
        
        # traj file button
        self.pushButton2 = QPushButton("Choose Trajectory File", self)
        self.tab1.layout.addWidget(self.pushButton2)
        
        # Label for Trajectory File (added directly below the button)
        self.labelTrajSelected = QLabel("No Trajectory file selected", self)
        self.tab1.layout.addWidget(self.labelTrajSelected)
        
        # Create the QLineEdits for frames and atoms with validation
        self.select_frames_tab1 = QLineEdit(self)
        self.select_frames_tab1.setPlaceholderText("Select Frames: e.g. 1-5 OR 1,2,3 OR 1-3,5")
        self.select_atoms_tab1 = QLineEdit(self)
        self.select_atoms_tab1.setPlaceholderText("Select Atoms: e.g. 1-5 OR 1,2,3 OR 1-3,5 ")
        
        # Add the input fields to the layout
        self.tab1.layout.addWidget(self.select_frames_tab1)
        self.tab1.layout.addWidget(self.select_atoms_tab1)
        
        # Define the regex pattern for frames and atoms
        pattern = r"^(\d+(-\d+)?(,\d+(-\d+)?)*,?)$"
        regex = QRegExp(pattern)
        
        # Apply QRegExpValidator for both inputs (frames and atoms)
        validator_frames_tab1 = QRegExpValidator(regex, self.select_frames_tab1)
        validator_atoms_tab1 = QRegExpValidator(regex, self.select_atoms_tab1)
        
        self.select_frames_tab1.setValidator(validator_frames_tab1)
        self.select_atoms_tab1.setValidator(validator_atoms_tab1)
        
        # Create Save button
        self.save_button_tab1 = QPushButton("Execute", self)
        self.save_button_tab1.clicked.connect(self.on_save_button_clicked_tab1)
        
        # Create a layout for the button and input fields
        input_layout_tab1 = QHBoxLayout()
        input_layout_tab1.addWidget(self.select_frames_tab1)
        input_layout_tab1.addWidget(self.select_atoms_tab1)
        input_layout_tab1.addWidget(self.save_button_tab1)
        
        # Add the layout to the tab
        self.tab1.layout.addLayout(input_layout_tab1)

        # Output labels (initially hidden)
        self.output_label_frames_tab1 = QLabel("Selected Frames: ", self)
        self.output_label_atoms_tab1 = QLabel("Selected Atoms: ", self)
        self.output_label_frames_tab1.setVisible(False)
        self.output_label_atoms_tab1.setVisible(False)
        
        self.tab1.layout.addWidget(self.output_label_frames_tab1)
        self.tab1.layout.addWidget(self.output_label_atoms_tab1)
        
        self.tab1.setLayout(self.tab1.layout)
#______________________________________________________________________________
        
        # Create Tab 2
        self.tab2.layout = QVBoxLayout()
        
        # topo file button
        self.pushButton3 = QPushButton("Choose Topology File", self)
        self.tab2.layout.addWidget(self.pushButton3)
        
        # Label for Topology File (added directly below the button)
        self.labelTopo2Selected = QLabel("No Topology file selected", self)
        self.tab2.layout.addWidget(self.labelTopo2Selected)
        
        # traj file button
        self.pushButton4 = QPushButton("Choose Trajectory File", self)
        self.tab2.layout.addWidget(self.pushButton4)
        
        # Label for Trajectory File (added directly below the button)
        self.labelTraj2Selected = QLabel("No Trajectory file selected", self)
        self.tab2.layout.addWidget(self.labelTraj2Selected)
        
        # Create the QLineEdits for frames and atoms with validation
        self.select_frames_tab2 = QLineEdit(self)
        self.select_frames_tab2.setPlaceholderText("Select Frames: e.g. 1-5 OR 1,2,3 OR 1-3,5")
        self.select_atoms_tab2 = QLineEdit(self)
        self.select_atoms_tab2.setPlaceholderText("Select Atoms: e.g. 1-5 OR 1,2,3 OR 1-3,5 ")
        
        # Add the input fields to the layout
        self.tab2.layout.addWidget(self.select_frames_tab2)
        self.tab2.layout.addWidget(self.select_atoms_tab2)
        
        # Define the regex pattern for frames and atoms
        pattern = r"^(\d+(-\d+)?(,\d+(-\d+)?)*,?)$"
        regex = QRegExp(pattern)
        
        # Apply QRegExpValidator for both inputs (frames and atoms)
        validator_frames_tab2 = QRegExpValidator(regex, self.select_frames_tab2)
        validator_atoms_tab2 = QRegExpValidator(regex, self.select_atoms_tab2)
        
        self.select_frames_tab2.setValidator(validator_frames_tab2)
        self.select_atoms_tab2.setValidator(validator_atoms_tab2)
        
        # Create Save button
        self.save_button_tab2 = QPushButton("Execute", self)
        self.save_button_tab2.clicked.connect(self.on_save_button_clicked_tab2)
        
        # Create a layout for the button and input fields
        input_layout_tab2 = QHBoxLayout()
        input_layout_tab2.addWidget(self.select_frames_tab2)
        input_layout_tab2.addWidget(self.select_atoms_tab2)
        input_layout_tab2.addWidget(self.save_button_tab2)
        
        # Add the layout to the tab
        self.tab2.layout.addLayout(input_layout_tab2)
        
        # Output labels (initially hidden)
        self.output_label_frames_tab2 = QLabel("Selected Frames: ", self)
        self.output_label_atoms_tab2 = QLabel("Selected Atoms: ", self)
        self.output_label_frames_tab2.setVisible(False)
        self.output_label_atoms_tab2.setVisible(False)
        
        self.tab2.layout.addWidget(self.output_label_frames_tab2)
        self.tab2.layout.addWidget(self.output_label_atoms_tab2)
        
        self.tab2.setLayout(self.tab2.layout)
#______________________________________________________________________________

        # Create Tab 3
        self.tab3.layout = QVBoxLayout()
        
        # pdb file button
        self.pushButton5 = QPushButton("Choose PDB File", self)
        self.tab3.layout.addWidget(self.pushButton5)
        
        # Label for PDB File (added directly below the button)
        self.labelPdbSelected = QLabel("No *.pdb file selected", self)
        self.tab3.layout.addWidget(self.labelPdbSelected)
        
        # csv file button
        self.pushButton6 = QPushButton("Choose CSV File", self)
        self.tab3.layout.addWidget(self.pushButton6)
        
        # Label for CSV File (added directly below the button)
        self.labelCsvSelected = QLabel("No *.csv file selected", self)
        self.tab3.layout.addWidget(self.labelCsvSelected)
        
        self.executeButton3 = QPushButton("Execute", self)
        self.tab3.layout.addWidget(self.executeButton3)
        
        
        self.tab3.setLayout(self.tab3.layout)
#______________________________________________________________________________
        
        # Connect buttons to file open methods
        self.pushButton1.clicked.connect(self.open_topo_file_dialog_tab1)
        self.pushButton2.clicked.connect(self.open_traj_file_dialog_tab1)
        self.pushButton3.clicked.connect(self.open_topo_file_dialog_tab2)
        self.pushButton4.clicked.connect(self.open_traj_file_dialog_tab2)
        self.pushButton5.clicked.connect(self.open_pdb_file_dialog)
        self.pushButton6.clicked.connect(self.open_csv_file_dialog)
        self.executeButton3.clicked.connect(self.execute_vispml)
        
        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
#______________________________________________________________________________

    def open_topo_file_dialog_tab1(self):
        # Filter to allow only Topology Files (e.g., PDB files)
        file_filter = "Topology Files (*.pdb; *.pdb.g; *.h5; *.lh5; *.prmtop; *.parm7; *.prm7; *.psf; *.mol2; *.hoomdxml; *.gro; *.arc; *.hdf5; *.gsd)"  # Update this filter as needed
        
        # Open the file dialog
        file, _ = QFileDialog.getOpenFileName(self, "Open Topology File", initial_dir, file_filter)
    
        if file:
            # Update the label to display the selected file path
            self.labelTopoSelected.setText(f"Selected Topology file: {file}")
        else:
            self.labelTopoSelected.setText("No Topology file selected")
            

    def open_traj_file_dialog_tab1(self):
        # Filter to allow only Trajectory Files
        file_filter = "Trajectory File (*.arc; *.dcd; *.binpos; *.xtc; *.trr; *.hdf5; *.h5; *.ncdf; *.netcdf; *.nc; *.pdb.gz; *.pdb; *.lh5; *.crd; *.mdcrd; *.inpcrd; *.restrt; *.rst7; *.ncrst; *.lammpstrj; *.dtr; *.stk; *.gro; *.xyz.gz; *.xyz; *.tng; *.xml; *.mol2; *.hoomdxml; *.gsd)"  # Update this filter as needed

        # Open the file dialog
        file, _ = QFileDialog.getOpenFileName(self, "Open Trajectory File", initial_dir, file_filter)

        if file:
            self.labelTrajSelected.setText(f"Selected Trajectory file: {file}")
        else:
            self.labelTrajSelected.setText("No Trajectory file selected") 
            
            
    def on_save_button_clicked_tab1(self):
        # Get the text from the QLineEdits and update the output labels
        entered_frames_tab1 = self.select_frames_tab1.text()
        entered_atoms_tab1 = self.select_atoms_tab1.text()
        
        if entered_frames_tab1:
            self.output_label_frames_tab1.setText(f"Selected Frames: {entered_frames_tab1}")
            self.output_label_frames_tab1.setVisible(True)
        
        if entered_atoms_tab1:
            self.output_label_atoms_tab1.setText(f"Selected Atoms: {entered_atoms_tab1}")
            self.output_label_atoms_tab1.setVisible(True)    
#______________________________________________________________________________
    
    def open_topo_file_dialog_tab2(self):
        # Filter to allow only Topology Files (e.g., PDB files)
        file_filter = "Topology Files (*.pdb; *.pdb.g; *.h5; *.lh5; *.prmtop; *.parm7; *.prm7; *.psf; *.mol2; *.hoomdxml; *.gro; *.arc; *.hdf5; *.gsd)"  # Update this filter as needed
        
        # Open the file dialog
        file, _ = QFileDialog.getOpenFileName(self, "Open Topology File", initial_dir, file_filter)
    
        if file:
            # Update the label to display the selected file path
            self.labelTopo2Selected.setText(f"Selected Topology file: {file}")
        else:
            self.labelTopo2Selected.setText("No Topology file selected")
            

    def open_traj_file_dialog_tab2(self):
        # Filter to allow only Trajectory Files
        file_filter = "Trajectory File (*.arc; *.dcd; *.binpos; *.xtc; *.trr; *.hdf5; *.h5; *.ncdf; *.netcdf; *.nc; *.pdb.gz; *.pdb; *.lh5; *.crd; *.mdcrd; *.inpcrd; *.restrt; *.rst7; *.ncrst; *.lammpstrj; *.dtr; *.stk; *.gro; *.xyz.gz; *.xyz; *.tng; *.xml; *.mol2; *.hoomdxml; *.gsd)"  # Update this filter as needed

        # Open the file dialog
        file, _ = QFileDialog.getOpenFileName(self, "Open Trajectory File", initial_dir, file_filter)

        if file:
            self.labelTraj2Selected.setText(f"Selected Trajectory file: {file}")
        else:
            self.labelTraj2Selected.setText("No Trajectory file selected") 
    
            
    def on_save_button_clicked_tab2(self):
        # Get the text from the QLineEdits and update the output labels
        entered_frames_tab2 = self.select_frames_tab2.text()
        entered_atoms_tab2 = self.select_atoms_tab2.text()
        
        if entered_frames_tab2:
            self.output_label_frames_tab2.setText(f"Selected Frames: {entered_frames_tab2}")
            self.output_label_frames_tab2.setVisible(True)
        
        if entered_atoms_tab2:
            self.output_label_atoms_tab2.setText(f"Selected Atoms: {entered_atoms_tab2}")
            self.output_label_atoms_tab2.setVisible(True)
            
#______________________________________________________________________________

    def open_pdb_file_dialog(self):
        # Filter to allow only PDB Files
        file_filter = "PDB Files (*.pdb)"
   
        # Open the file dialog
        self.file_pdb, _ = QFileDialog.getOpenFileName(self, "Choose PDB File", initial_dir, file_filter)
   
        if self.file_pdb:
            self.labelPdbSelected.setText(f"{self.file_pdb}")
        else:
            self.labelPdbSelected.setText("No PDB file selected")
         
            
    def open_csv_file_dialog(self):
        # Filter to allow only CSV File
        file_filter = "CSV File (*.csv)"
   
        # Open the file dialog
        self.file_csv, _ = QFileDialog.getOpenFileName(self, "Choose CSV File", initial_dir, file_filter)
        
   
        if self.file_csv:
            self.labelCsvSelected.setText(f"{self.file_csv}")
        else:
            self.labelCsvSelected.setText("No CSV file selected")

#______________________________________________________________________________

    def execute_simula(self):
        return

#______________________________________________________________________________
    
    def execute_reduce(self):
        return


#______________________________________________________________________________
     
    def execute_vispml(self):
        
        self.visualize = vispml.results_pdb(vispml.read_pdb(self.file_pdb), vispml.normalize_flux(self.file_csv))
        # print(self.file_pdb)
        # print(vispml.read_pdb(self.file_pdb))
        # print(self.file_csv)
        # print(vispml.normalize_flux(self.file_csv))
        # print(self.visualize)
        
        subprocess.run(["pymol", "results.pdb", "b-factor.pml"])
        
#______________________________________________________________________________
if __name__ == '__main__':
    app = QApplication([])
    p = Plugin()
    p.show()
    app.exec_()


 


# C:\Users\Medizinische Chemie\.conda\envs\aliaa_env\Lib\site-packages\pymol\plugins\aliaa_gui