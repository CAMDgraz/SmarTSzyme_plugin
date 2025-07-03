#!/bin/python

"""
@authors: Aliaa Abd Elhalim [aliaa.abdelhalim@edu.fh-joanneum.at]
          Daniel Platero-Rochart [daniel.platero-rochart@medunigraz.at]
          Pedro A. Sanchez-Murcia [pedro.murcia@medunigraz.at]
"""
# Imports ======================================================================
from PyQt5 import uic
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout,
                             QTableWidgetItem, QMessageBox, QButtonGroup,
                             QFileDialog, QProgressDialog)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure
import os
import numpy as np
import pandas as pd
from pymol import cmd
import seaborn as sns
from . import smd
#import prepare_smd as smd
import shutil
import glob
from . import smarTS

# Plots class ==================================================================
class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, placeholder, parent=None, dpi=75):
        fig_width = placeholder.width()/dpi
        fig_height = placeholder.height()/dpi

        self.fig = Figure(figsize=(fig_width, fig_height), dpi=dpi)
        super().__init__(self.fig)
    
    def update_plot(self, nrows, ncols):
        self.fig.clear()
        self.axes = []
        for i in range(nrows):
            ax = self.fig.add_subplot(nrows, ncols, i+1)
            self.axes.append(ax)
        self.fig.subplots_adjust(left=0.15, right=0.99, bottom=0.15, top=0.99,
                                 hspace=0.5) 
        self.draw()

# Main Window ==================================================================
class SmarTSWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        uifile = os.path.join(os.path.dirname(__file__),
                              'smarTS.ui')
        uic.loadUi(uifile, self)
        self.show()

        # Internal variables ===================================================
        self.ax_labels = {'DISTANCE':r'Distance $(\AA)$',
                          'ANGLE':r'Angle $(degrees)$',
                          'DIHEDRAL':r'Dihedral $(degrees)$'}
        
        self.measures_pd = pd.DataFrame()

        self.lineEdit_atoms = [self.lineEdit_atom1, self.lineEdit_atom2,
                               self.lineEdit_atom3, self.lineEdit_atom4]
        
        self.score = None

        # Radio buttons ========================================================
        self.radioChain = QButtonGroup(self)
        self.radioChain.addButton(self.check_bck)
        self.radioChain.addButton(self.check_sidechain)
        self.radioChain.addButton(self.check_all)

        self.radioFrames = QButtonGroup(self)
        self.radioFrames.addButton(self.check_automatic)
        self.radioFrames.addButton(self.check_manual)

        self.radioSMD = QButtonGroup(self)
        self.radioSMD.addButton(self.check_cumulative)
        self.radioSMD.addButton(self.check_exponential)

        self.radioTop = QButtonGroup(self)
        self.radioTop.addButton(self.check_topmax)
        self.radioTop.addButton(self.check_topmin)

        # Canvas for plots =====================================================
        self.layout_canvas = QVBoxLayout(self.widget)
        self.canvas_measures = MplCanvas(self.widget)
        toolbar = NavigationToolbar2QT(self.canvas_measures, self)
        self.layout_canvas.addWidget(toolbar)
        self.layout_canvas.addWidget(self.canvas_measures)

        self.layout_canvas2 = QVBoxLayout(self.widget_2)
        self.canvas_pmf = MplCanvas(self.widget)
        toolbar2 = NavigationToolbar2QT(self.canvas_pmf, self)
        self.layout_canvas2.addWidget(toolbar2)
        self.layout_canvas2.addWidget(self.canvas_pmf)

        self.layout_canvas3 = QVBoxLayout(self.widget_3)
        self.canvas_smart = MplCanvas(self.widget)
        toolbar3 = NavigationToolbar2QT(self.canvas_smart, self)
        self.layout_canvas3.addWidget(toolbar3)
        self.layout_canvas3.addWidget(self.canvas_smart)

        # Connections ==========================================================
        self.combo_measures.currentIndexChanged.connect(self.measure_changed)
        self.combo_plottype.currentIndexChanged.connect(self.plot_measure)
        self.button_plot.clicked.connect(self.add_measure)
        self.button_delete.clicked.connect(lambda: self.delete_entry(
                                                   self.tableWidget_measures))
        self.combo_cvtype.currentIndexChanged.connect(self.update_combo_cv)
        self.button_addcv.clicked.connect(self.add_cv)
        self.button_delcv.clicked.connect(lambda: self.delete_entry(
                                                  self.tableWidget_cv))
        self.button_addmask.clicked.connect(self.add_mask)
        self.button_delmask.clicked.connect(lambda: self.delete_entry(
                                                    self.tableWidget_mask))
        self.button_addfilter.clicked.connect(self.add_filter)
        self.button_delfilter.clicked.connect(lambda: self.delete_entry(
                                                      self.tableWidget_filter))
        self.button_delframe.clicked.connect(lambda: self.delete_entry(
                                                     self.tableWidget_frames))
        self.button_selframes.clicked.connect(self.select_frames)
        self.check_allmeasures.stateChanged.connect(lambda: self.sel_all(
                                                     self.tableWidget_measures))
        self.check_allcv.stateChanged.connect(lambda: self.sel_all(
                                                     self.tableWidget_cv))
        self.check_allmask.stateChanged.connect(lambda: self.sel_all(
                                                     self.tableWidget_mask))
        self.check_allfilters.stateChanged.connect(lambda: self.sel_all(
                                                     self.tableWidget_filter))
        self.check_allframes.stateChanged.connect(lambda: self.sel_all(
                                                     self.tableWidget_frames))
        self.button_updatemask.clicked.connect(self.update_mask)
        self.button_export.clicked.connect(self.export_csv)
        self.check_manual.toggled.connect(self.sele_mode)
        self.check_automatic.toggled.connect(self.sele_mode)
        self.button_mdpath.clicked.connect(self.load_file)
        self.button_topology.clicked.connect(self.load_file)
        self.button_smdlist.clicked.connect(self.load_file)
        self.button_generatesmd.clicked.connect(self.generate_smd)
        self.button_plotpmf.clicked.connect(self.plot_pmf)
        self.button_matrices.clicked.connect(self.load_file)
        self.button_coupling.clicked.connect(self.load_file)
        self.button_score.clicked.connect(self.mutability_score)
        self.button_plottop.clicked.connect(self.plot_score)
        self.button_plotexp.clicked.connect(self.plot_experimental)
        self.button_pymoltop.clicked.connect(self.show_pymol)

    def measure_changed(self, index):
        # Distance by default
        self.label_1.setEnabled(True)
        self.label_2.setEnabled(True)
        self.label_3.setEnabled(False)
        self.label_4.setEnabled(False)
       
        self.lineEdit_atom1.setEnabled(True)
        self.lineEdit_atom2.setEnabled(True)
        self.lineEdit_atom3.setEnabled(False)
        self.lineEdit_atom4.setEnabled(False)

        if index >= 1: # Angle selected
            self.label_3.setEnabled(True)       
            self.lineEdit_atom3.setEnabled(True)
        if index >= 2: # Dihedral selected
            self.label_4.setEnabled(True)
            self.lineEdit_atom4.setEnabled(True)

        return 
    
    def add_measure(self):
        # Type of measure
        index_comboBox = self.combo_measures.currentIndex()
        if index_comboBox == 0:
            measure_type, natoms =  'DISTANCE', 2
        elif index_comboBox == 1:
            measure_type, natoms = 'ANGLE', 3
        elif index_comboBox == 2:
            measure_type, natoms = 'DIHEDRAL', 4
        
        # Check and store atoms
        measure_atoms = []
        for edit_atom in self.lineEdit_atoms[:natoms]:
            atom_str = edit_atom.text().strip()
            try:
                atom = int(atom_str)
                measure_atoms.append(atom)
            except ValueError:
                QMessageBox.critical(self, "Error", "Wrong format of atoms")
                return

        # Get label
        measure_label = self.lineEdit_label.text().strip()
        if measure_label == "":
            QMessageBox.critical(self, "Error", "No label provided")
            return

        if self.duplicated_label(measure_label):
            return
        
        # Add to table
        table_row = self.tableWidget_measures.rowCount()
        self.tableWidget_measures.setRowCount(table_row + 1)
        check_sel = QTableWidgetItem()
        check_sel.setFlags(Qt.ItemFlag.ItemIsUserCheckable | 
                           Qt.ItemFlag.ItemIsEnabled)
        check_sel.setCheckState(Qt.CheckState.Unchecked)
        self.tableWidget_measures.setItem(table_row, 0,
                                            check_sel)
        self.tableWidget_measures.setItem(table_row, 1,
                                           QTableWidgetItem(measure_type))
        self.tableWidget_measures.setItem(table_row, 2,
                                           QTableWidgetItem(measure_label))
        atoms_cell = ",".join(str(atom) for atom in measure_atoms)
        self.tableWidget_measures.setItem(table_row, 3,
                                           QTableWidgetItem(atoms_cell))
        self.tableWidget_measures.resizeColumnsToContents()
                
        # Calculate and plot measure 
        measure = self.compute_measure(measure_type, measure_atoms,
                                       measure_label)
        self.measures_pd[measure_label] = measure
        #self.measures_pd[measure_label] = np.random.randint(1, 20, size=100) # testing
        self.plot_measure()

        # Clean line edits and update menus
        for edit_atom in self.lineEdit_atoms:
            edit_atom.setText('')
        self.lineEdit_label.setText('')
        self.update_combo_cv()
        self.update_combo_filter()
        return
    
    def duplicated_label(self, label_new):
        rows = self.tableWidget_measures.rowCount()
        for row in range(rows):
            label = self.tableWidget_measures.item(row, 2).text()
            if label == label_new:
                QMessageBox.critical(self, "Error", "Duplicated label")
                return True
        return False
    
    def plot_measure(self):
        # Create axes for each plot and get ax info
        ax_counter = []
        plot_info = [] # type, label, ax id and ax label
        for row in range(self.tableWidget_measures.rowCount()):
            measure_type = self.tableWidget_measures.item(row, 1).text()
            if measure_type not in ax_counter:
                ax_counter.append(measure_type)
            plot_info.append((measure_type, ax_counter.index(measure_type)))
        
        self.canvas_measures.update_plot(len(ax_counter), 1)

        # Get measures info
        plot_type = self.combo_plottype.currentIndex()
        for column, info in zip(self.measures_pd.columns, plot_info):
            ax = self.canvas_measures.axes[info[1]]
            if plot_type == 0:
                sns.kdeplot(x=self.measures_pd[column], ax=ax,
                            fill=True, label=column)
                ax.set(xlabel=self.ax_labels[info[0]], ylabel='Density')
            elif plot_type == 1:
                x = range(1, len(self.measures_pd[column]) + 1)
                sns.lineplot(x=x, y=self.measures_pd[column],
                             ax=ax, label=column)
                ax.set(xlabel=r'Frames', ylabel=self.ax_labels[info[0]])
            ax.legend(loc='upper right')
        self.canvas_measures.draw()
        return

    def export_csv(self):
        file, _ = QFileDialog.getSaveFileName(self, "Export csv")
        # Get measures
        if file:
            self.measures_pd.to_csv(file, index=False)
        else:
            QMessageBox.critical(self, "Error", "No output file selected")
        return
    
    def sel_all(self, table):
        sender = self.sender()
        state = sender.checkState()
        rows = table.rowCount()

        if state == Qt.CheckState.Checked:
            for row in range(rows):
                check_item = table.item(row, 0)
                check_item.setCheckState(Qt.Checked)

        else:
            for row in range(rows):
                check_item = table.item(row, 0)
                check_item.setCheckState(Qt.Unchecked)
        return
    
    def delete_entry(self, table):
        rows = table.rowCount()
        for row in range(rows-1, -1, -1):
            check_item = table.item(row, 0)
            if check_item.checkState() == Qt.CheckState.Checked:
                if table == self.tableWidget_mask:
                    cmd.hide('sticks', f'mask_{row}')
                    cmd.delete(f'mask_{row}')
                    
                    # Update mask names
                    try:
                        for maskid in range(row, rows):
                            cmd.set_name(f"mask_{maskid + 1}", f"mask_{maskid}")
                        for maskid in range(0, rows):
                            cmd.show('sticks', f'mask_{maskid}')
                    except:
                        pass

                if table == self.tableWidget_measures:
                    label = self.tableWidget_measures.item(row, 2).text()
                    self.measures_pd.drop(label, axis=1, inplace=True)
                    # Clean filters
                    rows_filter = self.tableWidget_filter.rowCount()
                    for row_filter in range(rows_filter-1, -1, -1):
                        filter_label = self.tableWidget_filter.item(row_filter,
                                                                     1)
                        if filter_label not in self.measures_pd.columns():
                            self.tableWidget_filter.removeRow(row_filter)
                            self.tableWidget_frames.setRowCount(0)
                    # Clean CVs
                    rows_cv = self.tableWidget_cv.rowCount()
                    for row_cv in range(rows_cv-1, -1, -1):
                        cv_label = self.tableWidget_cv.item(row_cv, 2).text()
                        cv_label = cv_label.split(',')
                        for cv_label_ in cv_label:
                            if cv_label_ not in self.measures_pd.columns:
                                self.tableWidget_cv.removeRow(row_cv)
                    # Delete pymol object
                    cmd.delete(label)

                table.removeRow(row)

        # Redraw plots
        if table == self.tableWidget_measures:
            self.plot_measure()
        return

    def compute_measure(self, measure_type, atoms, label):
        n_frames = cmd.count_states()
        measure_values = np.zeros(n_frames)

        if measure_type == "DISTANCE":
            for frame in range(n_frames):
                measure_values[frame] = cmd.get_distance(f'id {atoms[0]}',
                                                         f'id {atoms[1]}',
                                                         state=frame+1)
                cmd.distance(label, f"id {atoms[0]}", f"id {atoms[1]}")
                
        elif measure_type == "ANGLE":
            for frame in range(n_frames):
                measure_values[frame] = cmd.get_angle(f'id {atoms[0]}',
                                                      f'id {atoms[1]}',
                                                      f'id {atoms[2]}',
                                                      state=frame + 1)
                cmd.angle(label, f"id {atoms[0]}", f"id {atoms[1]}",
                          f"id {atoms[2]}")
        elif measure_type == "DIHEDRAL":
            for frame in range(n_frames):
                measure_values[frame] = cmd.get_dihedral(f'id {atoms[0]}',
                                                         f'id {atoms[1]}',
                                                         f'id {atoms[2]}',
                                                         f'id {atoms[3]}',
                                                         state=frame + 1)
                cmd.dihedral(label, f"id {atoms[0]}", f"id {atoms[1]}",
                             f"id {atoms[2]}", f"id {atoms[3]}")
        return measure_values
    
    def update_combo_cv(self):
        index = self.combo_cvtype.currentIndex()

        # Clean cv boxes
        self.combo_cv1.clear()
        self.combo_cv2.clear()

        # Populate with the measures
        for row in range(self.tableWidget_measures.rowCount()):
            measure_type = self.tableWidget_measures.item(row, 1).text()
            measure_label = self.tableWidget_measures.item(row, 2).text()
            if index == 0 and measure_type == 'DISTANCE':
                self.combo_cv1.addItem(measure_label)
            elif index == 1 and measure_type == 'ANGLE':
                self.combo_cv1.addItem(measure_label)
            elif index == 2 and measure_type == 'DIHEDRAL':
                self.combo_cv1.addItem(measure_label)
            elif index == 3 and measure_type == 'DISTANCE':
                self.combo_cv1.addItem(measure_label)
                # Activate combo cv2
                self.combo_cv2.setEnabled(True)
                self.lineEdit_coeff1.setEnabled(True)
                self.lineEdit_coeff2.setEnabled(True)
                self.combo_cv2.addItem(measure_label)
        return

    def add_cv(self):
        cv_type = self.combo_cvtype.currentText()
        if cv_type != "LCOD":
            cv_label = self.combo_cv1.currentText()
            cv_coeff = ""
        else:
            d1 = self.combo_cv1.currentText()
            d2 = self.combo_cv2.currentText()
            r1 = self.lineEdit_coeff1.text().strip()
            r2 = self.lineEdit_coeff2.text().strip()
            try:
                r1 = int(r1)
                r2 = int(r2)
            except ValueError:
                QMessageBox.critical(self, "Error",
                                     "Wrong format of Ceofficients") 
                return
            cv_coeff = f"{r1},{r2}"
            cv_label = f"{d1},{d2}"
        cv_final_value = self.lineEdit_value.text().strip()
        try:
            cv_final_value = float(cv_final_value)
        except ValueError:
            QMessageBox.critical(self, "Error", "Wrong format of Final Value")
            return
        cv_harm = self.lineEdit_harm.text().strip()
        if cv_harm == "":
            cv_harm = "600"
        try:
            int(cv_harm)
        except ValueError:
            QMessageBox.critical(self, "Error",
                                 "Wrong format of Harmonic constant")

        # Add to table
        table_row = self.tableWidget_cv.rowCount()
        self.tableWidget_cv.setRowCount(table_row + 1)

        check_sel = QTableWidgetItem()
        check_sel.setFlags(Qt.ItemFlag.ItemIsUserCheckable | 
                           Qt.ItemFlag.ItemIsEnabled)
        check_sel.setCheckState(Qt.CheckState.Unchecked)
        self.tableWidget_cv.setItem(table_row, 0, check_sel)
        self.tableWidget_cv.setItem(table_row, 1,
                                      QTableWidgetItem(cv_type))
        self.tableWidget_cv.setItem(table_row, 2,
                                      QTableWidgetItem(cv_label))
        self.tableWidget_cv.setItem(table_row, 3,
                                      QTableWidgetItem(f"{cv_coeff}"))
        self.tableWidget_cv.setItem(table_row, 4,
                                      QTableWidgetItem(f"{cv_final_value}"))
        self.tableWidget_cv.setItem(table_row, 5,
                                      QTableWidgetItem(f"{cv_harm}"))
        self.tableWidget_cv.resizeColumnsToContents()
        return

    def update_mask(self):
        selections_pymol = cmd.get_names('selections', 0)
        self.combo_mask.addItems(selections_pymol)
        return
    
    def add_mask(self):
        qmmm_mask = []
        selection = self.radioChain.checkedButton().text()
        mask_str = self.combo_mask.currentText()
        pymol_selection = cmd.get_model(f"{mask_str}").atom
        
        residues = set([atom.resi for atom in pymol_selection])

        if selection == 'All':
            mask_selection = 'resid '
            for resid in residues:
                qmmm_mask_ = f':{str(resid)}'
                qmmm_mask.append(qmmm_mask_)
                mask_selection += f"{resid}+"
        elif selection == 'Backbone':
            mask_selection = 'id '
            for resid in residues:
                qmmm_mask_ = f':{str(resid)}@'
                for atom in cmd.get_model(f'resid {resid}').atom:
                    if atom.name in ['N', 'H', 'CA', 'HA', 'C', 'O']:
                        mask_selection += f"{str(atom.id)}+"
                        qmmm_mask_ += f'{atom.name},'
                qmmm_mask.append(qmmm_mask_)
            
        elif selection == 'Side Chain':
            mask_selection = 'id '
            for resid in residues:
                qmmm_mask_ = f':{str(resid)}@'
                for atom in cmd.get_model(f'resid {resid}').atom:
                    if atom.name not in ['N', 'H', 'CA', 'HA', 'C', 'O']:
                        mask_selection += f'{str(atom.id)}+'
                        qmmm_mask_ += f'{atom.name},'
                qmmm_mask.append(qmmm_mask_)

        # Check for duplicates
        qmmm_mask_str = " " .join(str(atoms) for atoms in qmmm_mask)
        rows_mask = self.tableWidget_mask.rowCount()
        for row in range(rows_mask):
            mask = self.tableWidget_mask.item(row, 1).text()
            if mask == qmmm_mask_str:
                QMessageBox.critical(self, "Info", "Duplicated mask")
                return
        
        # Add to table
        table_row = self.tableWidget_mask.rowCount()
        self.tableWidget_mask.setRowCount(table_row + 1)

        check_sel = QTableWidgetItem()
        check_sel.setFlags(Qt.ItemFlag.ItemIsUserCheckable | 
                           Qt.ItemFlag.ItemIsEnabled)
        check_sel.setCheckState(Qt.CheckState.Unchecked)
        self.tableWidget_mask.setItem(table_row, 0, check_sel)
        self.tableWidget_mask.setItem(table_row, 1,
                                      QTableWidgetItem(qmmm_mask_str))
        self.tableWidget_mask.resizeColumnsToContents()
        
        # Show in pymol
        cmd.select(f"mask_{table_row}", f'{mask_selection}')
        cmd.show('sticks', f'mask_{table_row}')
        return

    def update_combo_filter(self):
        self.combo_filter.clear()
        for row in range(self.tableWidget_measures.rowCount()):
            measure_label = self.tableWidget_measures.item(row, 2).text()
            self.combo_filter.addItem(measure_label)
        return 
    
    def add_filter(self):
        filter_label = self.combo_filter.currentText().strip()
        filter_condition = self.combo_condition.currentText().strip()
        filter_threshold = self.lineEdit_filter.text().strip()
        try:
            float(filter_threshold)
        except ValueError:
            QMessageBox.critical(self, "Error",
                                 "Wrong format of the filter threshold")
            return
        # Add to table
        table_row = self.tableWidget_filter.rowCount()
        self.tableWidget_filter.setRowCount(table_row + 1)
        check_sel = QTableWidgetItem()
        check_sel.setFlags(Qt.ItemFlag.ItemIsUserCheckable | 
                           Qt.ItemFlag.ItemIsEnabled)
        check_sel.setCheckState(Qt.CheckState.Unchecked)
        self.tableWidget_filter.setItem(table_row, 0, check_sel)
        self.tableWidget_filter.setItem(table_row, 1,
                                          QTableWidgetItem(filter_label))
        self.tableWidget_filter.setItem(table_row, 2,
                                          QTableWidgetItem(filter_condition))
        self.tableWidget_filter.setItem(table_row, 3,
                                          QTableWidgetItem(str(filter_threshold)))
        self.tableWidget_filter.resizeColumnsToContents()
        return

    def select_frames(self):
        frames_label = self.radioFrames.checkedButton().text()
        if frames_label == "Manual":
            self.frames_list = self.lineEdit_frames.text().strip().split(" ")
            self.frames_list = np.asarray(self.frames_list, dtype=int)
            self.frames_list -= 1
        elif frames_label == "Automatic":
            self.frames_list = np.arange(0, cmd.count_states(), 1)
            # Read filters
            for row in range(self.tableWidget_filter.rowCount()):
                filter_label = self.tableWidget_filter.item(row, 1).text()
                filter_condition = self.tableWidget_filter.item(row, 2).text()
                filter_threshold = self.tableWidget_filter.item(row, 3).text()
                filter_threshold = float(filter_threshold)

                if filter_condition == ">":
                    frames_true = np.where(self.measures_pd[filter_label]
                                           > filter_threshold)[0]
                elif filter_condition == "<":
                    frames_true = np.where(self.measures_pd[filter_label]
                                           < filter_threshold)[0]
                elif filter_condition == "=":
                    frames_true = np.where(self.measures_pd[filter_label]
                                           == filter_threshold)[0]
                self.frames_list = np.intersect1d(self.frames_list, frames_true)
            
        if len(self.frames_list) == 0:
            QMessageBox.critical(self, "Error",
                                       "No frames have been selected")
        # Add to table
        self.tableWidget_frames.setRowCount(0)
        for frame in self.frames_list:
            table_row = self.tableWidget_frames.rowCount()
            self.tableWidget_frames.setRowCount(table_row + 1)
            check_sel = QTableWidgetItem()
            check_sel.setFlags(Qt.ItemFlag.ItemIsUserCheckable | 
                               Qt.ItemFlag.ItemIsEnabled)
            check_sel.setCheckState(Qt.CheckState.Unchecked)
            self.tableWidget_frames.setItem(table_row, 0, check_sel)
            frame_ = int(frame) + 1
            self.tableWidget_frames.setItem(table_row, 1,
                                            QTableWidgetItem(str(frame_)))
        return
    
    def sele_mode(self):
        frames_label = self.radioFrames.checkedButton().text()
        if frames_label == "Manual":
            # Disable the Automatic options
            self.combo_filter.setEnabled(False)
            self.combo_condition.setEnabled(False)
            self.lineEdit_filter.setEnabled(False)
            self.button_addfilter.setEnabled(False)
            self.button_delfilter.setEnabled(False)
            # Enable the Manual options
            self.lineEdit_frames.setEnabled(True)
        else:
            # Enable the Automatic options
            self.combo_filter.setEnabled(True)
            self.combo_condition.setEnabled(True)
            self.lineEdit_filter.setEnabled(True)
            self.button_addfilter.setEnabled(True)
            self.button_delfilter.setEnabled(True)
            # Disbale the Manual options
            self.lineEdit_frames.setEnabled(False)

    def load_file(self):
        sender = self.sender().objectName()
        initial_dir = "./"
        ext_filter = "All files (*)"
        if sender == "button_mdpath":
            self.rst_path = QFileDialog.getExistingDirectory(self,
                                                                "MD rst path")
        elif sender == "button_topology":
            self.topo, _ = QFileDialog.getOpenFileName(self,
                                                       "Load topology file",
                                                       initial_dir, ext_filter)
        elif sender == "button_smdlist":
            self.smd_list, _ = QFileDialog.getOpenFileName(self,
                                                           "sMD jobs path",
                                                           initial_dir,
                                                           ext_filter)
        elif sender == "button_matrices":
            self.matrices_path = QFileDialog.getExistingDirectory(self,
                                                             "Matrices Path")
        elif sender == "button_coupling":
            self.coupling_path = QFileDialog.getExistingDirectory(self,
                                                                  "Coupling Path")
        return

    def data_smd_files(self):
        # Check data completness
        self.charge = self.lineEdit_charge.text()
        if self.charge == "":
            QMessageBox.critical(self, "Error", "No charge provided")
            return
        else:
            try:
                _ = int(self.charge)
            except ValueError:
                QMessageBox.critical(self, "Error", "Wrong format of charge")
                return
        
        self.steps = self.lineEdit_steps.text()
        if self.steps == "":
            self.steps = 15000 # default number of steps
        else:
            try:
                _ = int(self.steps)
            except ValueError:
                QMessageBox.critical(self, "Error",
                                     "Wrong format of the number of steps")
                return
        
        self.time_step = self.lineEdit_timestep.text()
        if self.time_step == "":
            self.time_step = 0.0002 # default time step
        else:
            try:
                _ = int(self.time_step)
            except ValueError:
                QMessageBox.critical(self, "Error",
                                     "Wrong format of the time step")
                return
        
        self.theory = self.combo_theory.currentText()
        
        # Get data from tables Widgets
        table_measures = self.tableWidget_measures
        table_cv = self.tableWidget_cv
        table_masks = self.tableWidget_mask
        table_frames = self.tableWidget_frames

        # CV data
        cv_rows = table_cv.rowCount()
        if cv_rows == 0:
            QMessageBox.critical(self, "Error", "No CV defined")
            return
        self.cv_dict = dict()
        for row in range(cv_rows):
            cv_type = table_cv.item(row, 1).text()
            cv_key = table_cv.item(row, 2).text()
            cv_label = cv_key.split(",")
            cv_coeff = table_cv.item(row, 3).text().split(",")
            cv_final = table_cv.item(row, 4).text()
            cv_harm = table_cv.item(row, 5).text()
            self.cv_dict[cv_key] = (cv_type, cv_label, cv_coeff, cv_final,
                                    cv_harm)
        
        measure_rows = table_measures.rowCount()
        if measure_rows == 0:
            QMessageBox.critical(self, "Error", "No measures defined")
            return
        self.measure_dict = dict()
        for row in range(measure_rows):
            measure_type = table_measures.item(row, 1).text()
            measure_label = table_measures.item(row, 2).text()
            measure_atoms = table_measures.item(row, 3).text().split(",")
            self.measure_dict[measure_label] = (measure_type, measure_label,
                                                measure_atoms)
        
        frames_row = table_frames.rowCount()
        if frames_row == 0:
            QMessageBox.critical(self, "Error", "No frames selected")
            return
        self.frames_list = []
        for row in range(frames_row):
            frame = table_frames.item(row, 1).text()
            self.frames_list.append(int(frame))
        
        mask_rows = table_masks.rowCount()
        if mask_rows == 0:
            QMessageBox.critical(self, "Error", "No mask selected")
            return
        self.masks_list = []
        for row in range(mask_rows):
            mask = table_masks.item(row, 1).text().split(" ")
            self.masks_list.extend(mask)
        return
    
    def generate_smd(self):
        self.data_smd_files()
        output = QFileDialog.getExistingDirectory(self, "Output dir")
        if not output:
            QMessageBox.critical(self, "Error", "No output path selected")
            return
        progress = QProgressDialog("Creating files ...", None, 0, 0, self)
        progress.setWindowModality(Qt.ApplicationModal)
        progress.setCancelButton(None)
        progress.setMinimumDuration(0)
        progress.show()
        QApplication.processEvents()
        for frame in self.frames_list:
            try:
                os.mkdir(f"{output}/frame_{frame}")
            except:
                pass

            shutil.copyfile(self.topo, f"{output}/frame_{frame}/top.top")
            rst_files = glob.glob(f"{self.rst_path}/*{frame}*.rst")
            if len(rst_files) == 0:
                QMessageBox.critical(self, "Error",
                                     f"No rst file found for frame {frame}")
                continue
            shutil.copyfile(f"{rst_files[0]}",
                            f"{output}/frame_{frame}/frame.rst")

            frame_ = frame - 1
            smd.write_cv(f"{output}/frame_{frame}/cv.in", self.cv_dict,
                         self.measure_dict, self.measures_pd, frame_)
            smd.write_qmmm(f"{output}/frame_{frame}/qmmm.in",
                           self.masks_list, self.theory,
                           self.steps, self.time_step, self.charge,
                           self.lineEdit_amber.text().strip(),
                           300)
            smd.write_jobrun(f"{output}/frame_{frame}/runjob.sh",
                             self.lineEdit_amber.text().strip(),
                             self.lineEdit_openmpi.text().strip())
        with open(f"{output}/smd_list.txt", "w") as f:
                for frame in self.frames_list:
                    f.write(f"{output}/frame_{frame}\n")
        with open(f"{output}/run.sh", "w") as f:
            f.write('#!/bin/bash\n\n')
            f.write('while IFS= read -r line; do\n')
            f.write('cd ${line}\n')
            f.write('sbatch runjob.sh\n')
            f.write('done < smd_jobs.txt\n')
        progress.close()
        return

    def plot_pmf(self):
        if self.smd_list == "":
            QMessageBox.critical(self, "Error",
                                 "No List of smd jobs provided")
            return
        
        smd_time = self.lineEdit_smdtime.text().strip()
        if smd_time == "":
            QMessageBox.critical(self, "Error", "No sMD time provided")
            return
        try:
            smd_time = float(smd_time)
        except ValueError:
            QMessageBox.critical(self, "Error", "Wrong format of the sMD time")
            return
        
        smd_step = self.lineEdit_smdstep.text().strip()
        if smd_step == "":
            QMessageBox.critical(self, "Error",
                                 "No sMD time step provided")
            return
        try:
            smd_step = float(smd_step)
        except ValueError:
            QMessageBox.critical(self, "Error", "Wrong format of sMD time step")
        
        smd_temp = self.lineEdit_smdT.text().strip()
        if smd_temp == "":
            smd_temp = 300
        else:
            try:
                smd_temp = float(smd_temp)
            except ValueError:
                QMessageBox.critical(self, "Error",
                                     "Wrong format of Temperature")
                return
        
        # Get PMF and plot
        smd_type = self.radioSMD.checkedButton().text()
        self.canvas_pmf.update_plot(1, 1)
        m_energy, m_time = smd.smd_results(self.smd_list, smd_time, smd_step,
                                           smd_type, 300,
                                           self.canvas_pmf.axes[0])
        if m_energy and m_time:
            self.label_maximum.setText(f"{m_energy:.2f}")
            self.label_time.setText(f"{m_time:.2f}")
        else:
            QMessageBox.critical(self, "Info", "No maximum found")
        self.canvas_pmf.draw()
        return

    def mutability_score(self):
        if self.matrices_path == "":
            QMessageBox.critical(self, "Error", "No path to matrices found")
            return
        if self.coupling_path == "":
            QMessageBox.critical(self, "Error", "No path to coupling found")
            return
        
        self.batch = self.lineEdit_batch.text().strip()
        if self.batch == "":
            self.batch = 1
        else:
            try:
                self.batch = int(self.batch)
            except ValueError:
                QMessageBox.critical(self, "Error",
                                     "Wrong format of Batch size")
                return

        self.nres = self.lineEdit_nres.text().strip()
        if self.nres == "":
            QMessageBox.critical(self, "Error", "No. of residues not found")
            return
        else:
            try:
                self.nres = int(self.nres)
            except ValueError:
                QMessageBox.critical(self, "Error",
                                     "Wrong format of the No. of residues")
                return
        
        progress = QProgressDialog("Calculating Score ...", None, 0, 0, self)
        progress.setWindowModality(Qt.ApplicationModal)
        progress.setCancelButton(None)
        progress.setMinimumDuration(0)
        progress.show()
        QApplication.processEvents()
        self.score, self.batches = smarTS.score(self.coupling_path,
                                                self.matrices_path,
                                                self.nres, self.batch)
        progress.close()
        return 
    
    def plot_score(self):
        self.topn = self.lineEdit_topn.text().strip()
        if self.topn == "":
            self.topn = 10
        else:
            try:
                self.topn = int(self.topn)
            except ValueError:
                QMessageBox.critical(self, "Error",
                                     "Wrong format of Top N")
                return
        
        # Plot top residues
        self.canvas_smart.update_plot(1, 1)
        top_plot = self.radioTop.checkedButton().text()
        try:
            smarTS.plot_score(self.score, self.batches, self.topn, self.nres,
                              self.canvas_smart.axes[0], top_plot)
        except AttributeError:
            QMessageBox.critical(self, "Error",
                                 "Calculate Mutability Score first")
            return
        
        self.canvas_smart.draw()
        return

    def plot_experimental(self):
        expres = self.lineEdit_expRes.text().strip().split()
        try:
            expres = [int(res) for res in expres]
        except ValueError:
            QMessageBox.critical(self, "Error",
                                 "Wrong format of experimental residues")
        # Plot experimental residues
        self.canvas_smart.update_plot(1, 1)
        smarTS.plot_experimental(expres, self.score, self.batches,
                                 self.canvas_smart.axes[0])
        self.canvas_smart.draw()
        return

    def show_pymol(self):
        self.topn = self.lineEdit_topn.text().strip()
        if self.topn == "":
            self.topn = 10
        else:
            try:
                self.topn = int(self.topn)
            except ValueError:
                QMessageBox.critical(self, "Error",
                                     "Wrong format of Top N")
                return
            
        top_stab, top_dest = smarTS.get_top(self.score, self.topn, self.nres)
        top_stab_id = np.where(top_stab[-1] == 1)[0]
        top_dest_id = np.where(top_dest[-1] == 1)[0]

        top_stab_pymol = 'resid '
        for resid in top_stab_id:
            resid += 1
            top_stab_pymol += f"{resid}+"

        top_dest_pymol = 'resid '
        for resid in top_dest_id:
            resid += 1
            top_dest_pymol += f"{resid}+"
        # Create pymol selection and show as sticks
        cmd.select(f"topStab", f"{top_stab_pymol}")
        cmd.select(f"topDest", f"{top_dest_pymol}")

        cmd.show("sticks", "topStab")
        cmd.show("sticks", "topDest")

        # Change the beta factor
        for resid, score in zip(range(1, self.nres + 1), self.score[-1]):
            cmd.alter(f"resid {resid}", f"b={score}")
        cmd.spectrum("b", "blue_white_red", "all", maximum=self.score.max(),
                     minimum=self.score.min())
        return
        
if __name__ == '__main__':
    app = QApplication([])
    window = SmarTSWindow()
    app.exec_()
