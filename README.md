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

After openning the plugin the user is asked to load the topology and trajectory files. All formats accepted by [MDTraj](https://mdtraj.org/1.9.4/index.html) are also available in the plugin. Once the trajectory is loaded, the user can select between three different windows under the *Tools* menu:

1. *Tools -> Measures:* Using atoms id we offer the possibility to plot and save the distribution of different measures (i.e. distances, angles, dihedrals and RMSD).
2. *Tools -> Simula:* Contains three tabs dedicated to the setup and results of sMD simulations. In the first one the user can define the collective variables (cv) and the QM/MM settings. In the second tab a second topology and trajectory are requested (system including solvent, ions etc.) and the user can select the starting frames for the sMD simulations. The frame selection is donen either by manually specifing the desires snapshots or using measures like distances and angles to select a subset of frames that comply with the conditions. The plugin will output each frame in a separate folder including already the all necesary files to run the calculations. Once the sMD calculations were runned, the user can make use of the third tab to visualize and identify the maximum in the mean free energy.
3. *Tools -> Reduce:* Using the output from [SmarTSzyme](https://github.com/CAMDGraz/SmarTSzyme), the user can visualize the Top N residues identified and display a color coded structure of the stabilizing (blue) and destabilizing (red) residues.

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

