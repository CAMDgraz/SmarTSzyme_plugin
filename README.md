# SmarTSzyme-GUI
<div align="center">
<img src="images/header.png" alt="load" width="1000"/>
</div>

*SmarTSzyme-GUI* is a Graphical User Interface to aid in the selection and setup of steered QM/MM MD simulations. At the current release it only supports Amber simulation software.

## Requirements
For using the plugin you need to create a new conda environment as follows (PyMOL will also be installed in the environment):
```bash
$ conda env create -f environment.yml
$ conda activate smarts_pymol 
```
Once the environment has been created execute the following command to zip the plugin:
```bash
$ sh pack.sh
```

## Installation
The installation of the plugin, is done via the plugin manager as follows:\
1. Launch the PyMOL installed with the created environment (/path/to/env/bin/pymol).\
2. Open the plugin manager *Plugin -> Plugin Manager*.\
3. In the *Install New Plugin* tab click the *Choose File...* button. Select the **plugin.zip** file created before.

## How to use
The plugin has been designed to aid the user in the first step of the workflow followed in [SmarTSzyme](https://github.com/CAMDGraz/SmarTSzyme) along with the visualization of the results.

The initial window of the SmarTSzyme-GUI consists of four options: Measurements, sMD setup, sMD results, and SmarTSzyme. 



**Measurements**. With this menu, geometrical features like distances, angles and dihedrals can be monitores on the GUI. After loading a trajectory as an object into PyMOL, the user can select the atom IDs to be traced along the trajectory. Each of the measurements need to be labeled. The different measurements can be plotted simultaneously on the GUI. As output, the use can select to plot these data as a distribution of the value or as a evolution along the different frames of the trajectory. These data can be exported as images. 

**sMD setup**. SmarTSzyme analyzes several chemical dynamics trajectories of the simulation of the catalyzed reaction to identify the residues that explain the (de)stabilization of the ETS with respect to the ES complex. In the current implementation, SmarTSzyme reads steered QM/MM (sMD) simulations. Thus, the menu sMD setup helps the users to analyze previous trajectories of the substrate-enzyme complex to select those snapshots to be used as initial geometries for the sMD calculations and to prepare the input files that are needed for running several sMD's using the QM/MM interface of Amber (https://ambermd.org). 

It consists of 4 submenus: Collective Variables, QM region and simulation parameters, Frames selection, and MD files.

COLLECTIVE VARIABLES
QM REGION AND SIMULATION PARAMETERS
FRAMES SELECTION
MD FILES

**sMD results**. The user needs to collect the files with the work values for the explored collective variables. The name files need to be listed in a file and ping pointed with the 'List' option. It needs also to be specified the total ammount of time of the sMD (i.e., 2.0 ps), the separation between the printed values of work (i.e., 0.01 ps), and the temperature. By clicking 'PMF', the GUI prints all profiles and it calculates the average using the Jarzynski equality. There are two ways of averaging: exponential or cumulative. The GUI identifies the energy maximum and at which time point is found. 

**SmarTSzyme**. This is the way to plot the results after one SmarTSzyme analysis. Before running the GUI, it is needed to calculate the matrices needed by SmarTSzyme with the command-line code. Please, check [SmarTSzyme](https://github.com/CAMDGraz/SmarTSzyme) to know how to run SmarTSzyme. 

## License
*SmarTSzyme-GUI* is licensed under GNU General Public License v3.0.

## Citation
The corresponding publication is under preparation

## Contact
**Laboratory of Computed-Aided Molecular Design Graz (Sánchez-Murcia's group)**

Division of Medicinal Chemistry\
Otto-Loewi Research Center\
Medical University of Graz\
Neue Stiftingstalstraße 6/III\
A-8010 Graz, Austria
 
In case of questions and/or suggestions you can contact us at: pedro.murcia@medunigraz.at and daniel.platero-rochart@medunigraz.at

<img src="images/logo.png" alt="load" width="200"/>

