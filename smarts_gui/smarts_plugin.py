#!/bin/python

"""
@authors: Aliaa Abd Elhalim [aliaa.abdelhalim@edu.fh-joanneum.at]
          Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sanchez-Murcia [pedro.murcia@medunigraz.at]
"""
# Imports ======================================================================
# PyQt5
from PyQt5 import uic
from PyQt5.QtWidgets import (QApplication, QMainWindow)

# Plugin specific
from . import simula
from . import reduce
from . import measures
from . import functions as fc

# Generals
import os
from pymol import cmd

# for testing
# import functions as fc
# import measures
# import reduce
# import simula

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
        self.actionSimula.triggered.connect(self.show_Simula)
        self.actionReduce.triggered.connect(self.show_Result)
        self.button_trajectory.clicked.connect(self.load_file)
        self.button_topology.clicked.connect(self.load_file)
        self.button_load.clicked.connect(self.load_system)

    # Windows
    def show_Measures(self):
        '''
        Shows the Measures window
        '''

        if not self.traj:
            fc.pop_error("Error!!!", "No system loaded")
            pass
        elif self.traj:
            self.measureWindow = measures.MeasureWindow(self.traj)
            self.measureWindow.show()
    
    def show_Simula(self):
        '''
        Shows the Simula input window
        '''
        self.simulaWindow = simula.SimulaWindow(self.traj)
        self.simulaWindow.show()
    
    def show_Result(self):
        '''
        Show the Reduce window
        '''
        self.reduceWindow = reduce.ReduceWindow(self.traj)
        self.reduceWindow.show()
       
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
        self.traj, need_top = fc.load_traj(self.traj_file, self.top_file)      
        if self.traj:
            self.label_system.setText(str(self.traj))
            system_name = os.path.splitext(os.path.basename(self.top_file))[0]
            pymol_objects = cmd.get_names('objects', 0)
            if system_name not in pymol_objects:
                if not need_top:
                    cmd.load(self.traj_file)
                else:
                    cmd.load(self.top_file)
                    cmd.load_traj(self.traj_file)
        return

if __name__ == '__main__':
    app = QApplication([])
    window = MainWindow()
    app.exec_()
