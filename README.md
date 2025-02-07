# Simula: A PyMOL plugin for steered Molecular Dynamics and the SmarTSzyme.

*Simula* is a Graphical User Interface to aid in the selection and setup of steered QM/MM MD simulations. At the current release it only supports Amber simulation software.

## How to use
The plugin has been designed to aid the user in the first step of the workflow followed in [SmarTSzyme](https://github.com/CAMDGraz/SmarTSzyme) along with the visualization of the results.

After openning the plugin the user is asked to load the topology and trajectory files. All formats accepted by [MDTraj](https://mdtraj.org/1.9.4/index.html) are also available in **Simula**. Once the trajectory is loaded, the user can select between three different windows under the *Tools* menu:

1. *Measures:* Using atoms id we offer the possibility to plot and save the distribution of different measures (i.e. distances, angles, dihedrals and RMSD).
2. *Simula:* Contains three tabs dedicated to the setup and results of sMD simulations. In the first one the user can define the collective variables (cv) and the QM/MM settings. In the second tab a second topology and trajectory are requested (system including solvent, ions etc.) and the user can select the starting frames for the sMD simulations. The frame selection is donen either by manually specifing the desires snapshots or using measures like distances and angles to select a subset of frames that comply with the conditions. The plugin will output each frame in a separate folder including already the all necesary files to run the calculations. Once the sMD calculations were runned, the user can make use of the third tab to visualize and identify the maximum in the mean free energy.
3. *Reduce:* Using the output from [SmarTSzyme](https://github.com/CAMDGraz/SmarTSzyme), the user can visualize the Top N residues identified and display a color coded structure of the most prominent residues.



