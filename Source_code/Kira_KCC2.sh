#!/bin/sh
#SBATCH --account=humbio
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --job-name=ORCHID_KCC2
#SBATCH --mail-user=joseph.raimondo@uct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL


#Load the software you have installed using the following 
module load python/miniconda3-py39
source activate single-cell-simulator

python SCS_Controller_v1_KCC2.py