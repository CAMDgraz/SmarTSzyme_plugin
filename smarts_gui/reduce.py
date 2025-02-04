#!/bin/python

"""
@authors: Aliaa Abd Elhalim [aliaa.abdelhalim@edu.fh-joanneum.at]
          Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sanchez-Murcia [pedro.murcia@medunigraz.at]
"""
# Imports ======================================================================
# PyQt5
from PyQt5 import uic
from PyQt5.QtWidgets import (QMainWindow, QFileDialog)
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT

# Plugin specific
from . import functions as fc

# Generals
import os
import heapq
import numpy as np
import pandas as pd
import seaborn as sns
from pymol import cmd

# for testing
# import functions as fc

# Reduce Window ================================================================
class ReduceWindow(QMainWindow):
    def __init__(self, traj):
        super().__init__()
        uifile = os.path.join(os.path.dirname(__file__),
                              'ui_files/reduce.ui')
        uic.loadUi(uifile, self)

        # Init variables
        self.pdb = ''
        self.csv = ''

        self.canvas1 = fc.MplCanvas(self, dpi=100)
        toolbar1 = NavigationToolbar2QT(self.canvas1, self)
        self.verticalLayout_2.addWidget(toolbar1)
        self.verticalLayout_2.addWidget(self.canvas1)
        self.ax1 = self.canvas1.fig.add_subplot(111)


        # Conections ===========================================================
        self.button_pdb.clicked.connect(self.load_file)
        self.button_csv.clicked.connect(self.load_file)
        self.button_pymol.clicked.connect(self.show_pymol)
        self.button_plot.clicked.connect(self.plot)

    def load_file(self):
        '''
        Load pdb and csv files
        '''
        sender = self.sender().objectName()

        initial_dir = './'
        ext_filter = "All files (*)"

        file, _ = QFileDialog.getOpenFileName(None,
                                                  "Load File",
                                                  initial_dir,
                                                  ext_filter)
        ext = os.path.splitext(os.path.basename(file))[1]
        if sender == 'button_pdb':
            if ext != '.pdb':
                fc.pop_error("Error!!!", "Not a pdb file")
            else:
                self.pdb = file
        elif sender == 'button_csv':
            if ext != '.csv':
                fc.pop_error("Error!!!", "Not a csv file")
            else:
                self.csv = file

    def show_pymol(self):
        norm_csv = fc.normalize_flux(self.csv)
        nest_pdb = fc.read_pdb(self.pdb)
        result_pdb = fc.results_pdb(nest_pdb, norm_csv)
        cmd.load(result_pdb)
        cmd.spectrum('b', selection='results')
    
    def plot(self):
        for_plotting = []
        heapq.heapify(for_plotting)
        try:
            csv_df = pd.read_csv(self.csv)
            residues = np.asarray(csv_df.residue)
            fluxes = np.asarray(csv_df.flux)
            stab_index = np.where(fluxes < 0)[0]
            destab_index = np.where(fluxes > 0)[0]
        except:
            fc.pop_error("Error!!!", "Could not read csv file")
            return
        
        ordered_destab = []
        ordered_stab = []
        for_plotting_destab = []
        for_plotting_stab = []
        heapq.heapify(ordered_destab)
        heapq.heapify(ordered_stab)
        heapq.heapify(for_plotting_destab)
        heapq.heapify(for_plotting_stab)

        for res, flux in zip(residues[stab_index], fluxes[stab_index]):
                heapq.heappush(ordered_stab, (flux, res))
        for res, flux in zip(residues[destab_index], fluxes[destab_index]):
                heapq.heappush(ordered_destab, (flux, res))

        # Get top stabilizing 
        try:
            top_lower = int(self.edit_lowest.text().strip().split()[0])
            lowest = heapq.nsmallest(top_lower, ordered_stab)
        except:
            lowest = ordered_stab

        self.edit_show1.clear()    
        for res_flux in lowest:
            self.edit_show1.insertPlainText(f'{res_flux[1]} -> {res_flux[0]}\n')
            heapq.heappush(for_plotting_stab, (res_flux[1], res_flux[0]))

        # Get top destabilizing            
        try:
            top_higher =int(self.edit_higher.text().strip().split()[0])
            highest = heapq.nlargest(top_higher, ordered_destab)
        except:
            highest = ordered_destab
        
        self.edit_show2.clear()
        for res_flux in highest:
            self.edit_show2.insertPlainText(f'{res_flux[1]} -> {res_flux[0]}\n')
            heapq.heappush(for_plotting_destab, (res_flux[1], res_flux[0]))

        # plot destabilizing and stabilizing
        self.ax1.clear()
        self.canvas1.draw()

        x_ax_stab = [res_flux[0] for res_flux in for_plotting_stab]
        y_ax_stab = [res_flux[1] for res_flux in for_plotting_stab]

        sns.barplot(x=x_ax_stab, y=y_ax_stab, ax=self.ax1, color='blue',
                     label='Stabilizing')
        
        x_ax_destab = [res_flux[0] for res_flux in for_plotting_destab]
        y_ax_destab = [res_flux[1] for res_flux in for_plotting_destab]
        sns.barplot(x=x_ax_destab, y=y_ax_destab, ax=self.ax1, color='red',
                     label='Destabilizing')
        
        xticks = self.ax1.get_xticks()
        xticklabels = self.ax1.get_xticklabels()
        self.ax1.set_xticks(xticks)
        self.ax1.set_xticklabels(xticklabels, rotation=45, fontsize=7)

        yticks = self.ax1.get_xticks()
        yticklabels = self.ax1.get_xticklabels()
        self.ax1.set_xticks(yticks)
        self.ax1.set_xticklabels(yticklabels, fontsize=7)


        self.ax1.set_xlabel(r'Residues')
        self.ax1.set_ylabel(r'Flux')
        self.ax1.legend(loc='upper right')
        self.canvas1.draw()


        

        
        
        
            
