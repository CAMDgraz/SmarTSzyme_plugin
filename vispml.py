#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 10:37:25 2024
@author: aliaa
""" 

# whether it's normalized, not normalized, range -1 to 1 or 0 to 1 or 0 to 100
# the values themselves are very close to each other and always show the same colors in pymol

import pandas as pd
import pymol
from pymol import cmd

range_min = 0
range_max = 1

path_to_pdb = r'C:\Users\Medizinische Chemie\Desktop\aliaa\ndt.pdb'
path_to_csv = r'C:\Users\Medizinische Chemie\Desktop\aliaa\results.csv'

def normalize_flux(csv_file):
    
    """
    Normalizes residue dependent flux values from csv file
    
    Parameters
    ----------
    csv file: str
        Path to csv file 
        contains residues (int), fluxes (float)
                    
    Return
    ------
    df['flux_norm']: DataFrame
        flux normalized with min-max normalization
        rounded to 4 decimal places
    """

    df = pd.read_csv(csv_file)
    df.loc[df['flux']!= 0, 'flux_norm'] = ((range_min) + ((df['flux'] - df['flux'].min()) * (range_max - (range_min))/ (df['flux'].max() - df['flux'].min()))).round(4)
    df.loc[df['flux_norm'].isna(), 'flux_norm'] = (range_max - (- range_min)) / 2 #middle value

    
    
    return df['flux_norm']
        

def read_pdb(pdb_file):
    
    """
    Reads pdb file as a nested list
    
    Parameters
    ----------
    pdb file: str
        Path to pdb file 
        
    Return
    ------
    pdb_list: nested list with specified coloumns
    """ 
    
    pdb_list = []
    with open(pdb_file) as pdbfile:
        for line in pdbfile:
            if line[:4] == 'ATOM' or line[:6] == "HETATM":
            
                # Split the line
                split_line = [line[:6],             # ATOM / HETATM (str)
                              int(line[6:11]),      # Atom serial number (int)
                              line[12:16],          # Atom name (str)
                              line[17:20],          # Residue name (str)
                              int(line[22:26]),     # Residue sequence number (int)
                              float(line[30:38]),   # X orthogonal Å coordinate (float)
                              float(line[38:46]),   # Y orthogonal Å coordinate (float)
                              float(line[46:54]),   # Z orthogonal Å coordinate (float)
                              float(line[54:60]),   # Occupancy (float)
                              float(line[60:66]),   # Temperature factor (float)
                              line[76:78]]          # Element symbol (str)                          
                
                pdb_list.append(split_line)
                
    return pdb_list             
    

def results_pdb(pdb_list,csv_flux):

    """
    Replaces beta-factors from original pdb 
    with previously created pdb nested list and calculated normalized flux values
    
    Parameters
    ----------
    pdb_list: nested list
        
    flux_norm: list of str
                
    Return
    ------
    mod_pdb: str 
        generates pdb file with modified beta-factors 
        depending on residue sequence number
    """
     
    for i in range(len(pdb_list)): 
        for j in range(len(csv_flux)):
            if j+1 == pdb_list[i][4]:
                pdb_list[i][9] = csv_flux[j]
    
    with open('results.pdb', 'w') as mod_pdb: 
        for i in range(len(pdb_list)):
            pdb_list[i][0] = pdb_list[i][0].ljust(6)
            pdb_list[i][1] = str(pdb_list[i][1]).rjust(5)
            pdb_list[i][2] = pdb_list[i][2].ljust(4)
            pdb_list[i][3] = pdb_list[i][3].rjust(3)
            pdb_list[i][4] = str(pdb_list[i][4]).rjust(6)
            pdb_list[i][5] = f"{pdb_list[i][5]:8.3f}".rjust(8)
            pdb_list[i][6] = f"{pdb_list[i][6]:8.3f}".rjust(8)
            pdb_list[i][7] = f"{pdb_list[i][7]:8.3f}".rjust(8)
            pdb_list[i][8] = f"{pdb_list[i][8]:6.2f}".rjust(6)
            pdb_list[i][9] = str(pdb_list[i][9]).rjust(6)
            pdb_list[i][10] = pdb_list[i][10].ljust(2)
            mod_pdb.write("%s%s %s %s%s    %s%s%s%s%s          %s\n"
                          % (pdb_list[i][0],pdb_list[i][1],pdb_list[i][2],
                             pdb_list[i][3],pdb_list[i][4],pdb_list[i][5],
                             pdb_list[i][6],pdb_list[i][7],pdb_list[i][8],
                             pdb_list[i][9],pdb_list[i][10]))
    
        
    return mod_pdb

# a = read_pdb(path_to_pdb)
# b = normalize_flux(path_to_csv)
# c = results_pdb(a,b)
