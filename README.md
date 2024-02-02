# single-cell-simulator

This repository contains the code necessary to reproduce supporting simulations and related figures in:
> All-optical reporting of inhibitory receptor driving force in the nervous system, by
> Selfe, Steyn, Shorer, Burman, Düsterwald, Abdelfattah, Schreiter, Newey, Akerman and Raimondo

## :scroll: Structure
```
├── Single-Cell-Simulator     <- This repository
|   ├── Source_code          <- Source files for analysis and figures
|   |   └── ...
|   ├── hdf5_files          <- Simulation data (checkpoints)
|   |   └── ...
|   ├── Jupyter_notebooks_for_analysis     <- Notebooks for analysis and figures
|   |   └── ...
|   └── Jupyter_outputs                <- Figures from notebooks
```

## :clipboard: Requirements

**OS**
Tested on Windows 11 Build 22631.3085
Tested on macOS 14 Sonoma

**Python**
Python Version 3.12.1.

**Python Packages**
h5py==3.10.0
ipykernel==6.26.0
ipywidgets==8.1.1
matplotlib==3.8.1
numpy==1.26.1
seaborn==0.13.0

## :minidisc: Installation

The repository is coded in python with Jupyter notebooks for ease of interaction. Requirements are minimal, but included in requirements.txt.

The Jupyter Notebooks can be run by via the web browser using the jupyter notebook site: [Jupyter Notebook](https://jupyter.org/).
Alternatively, it can be installed locally using Python's package manager, pip, or with Anaconda which is a free an open source Python distribution that includes Jupyter. Installation and setup using this method would take approximately 10 to 15 minutes. 

## :computer: Instructions and Demo

**New Simulations**

To run a new simulation open the SCS_Controller.py file located in the Source_code folder. This python script has step by step instructions to initialize the parameters to run a biophysically realistic single cell neuronal simulation. 

The run time of the simulation varies depending on the total duration of the simulation, the time step, and the presence or absence of additional features such as synapses. 

We recommend using a time step of between 10^-4 and 10^-6, and simulation lengths of between 1s to 10000s. Simuation lengths that are longer than these may result in simulation times exceeding 24 hours. For standard simulations within the parameters mentioned, the simulation time is typically on order of minutes to hours.

All simulations are saved as an HDF file with the file name specified by the user when running the SCS_Controller.py script. 
Simulations can be run, beginning afresh, by specifying the initialization parameters and setting the sim_type variable = "New". 
Alternatively, simulations can be initialized from the final values from an previously run simulation. This requires the HDF file from the previously run simulation to be placed in the working directory, and the sim_type variable = "Extend".

An example of an already run simulation is placed under the folder hdf5_files.


**Plotting results**

To plot the results including ion dynamics of a simulation use the Jupyter notebook file ORCHID_Dynamics.ipynb. If the simulation includes additional synapses use the ORCHID_Synapses.ipynb file. These files will required the input of the name of the HDF5 file where the simulation results are stored. The HDF5 file to be analyzed need to be in the folder "HDF5_files". The figures take a few seconds to minutes render based on the size of the HDF5 file. 

Analysis of simulation data are reproducible through the sharing of HDF5 files produced by the software.

<!--
## :loudspeaker: Citation

If you are using this code, please cite [the paper](https://www.biorxiv.org/content/10.1101/2023.08.30.555464v2):

    @article{selfe2023,
      title={All-optical reporting of inhibitory receptor driving force in the nervous system},
      author={Jakob, AMV and Gershman, SJ},
      author={Selfe, JS and Steyn, TJS and Shorer, EF and Burman, RJ and Düsterwald, KM and Abdelfattah, AS and Schreiter, ER and Newey, SE and Akerman, CJ and Raimondo, JV}
      journal={},
      year={2023}
    }
-->

## :envelope: Contact

If you have any questions, please get in touch with the corresponding author
[Joseph Raimondo](mailto:joseph.raimondo@gmail.com).