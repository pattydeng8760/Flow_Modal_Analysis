####################################################################################
# Applying Dynamic Mode Analysis to transient data
# Input Parameters 
# Author: Patrick Deng
####################################################################################
import os
import numpy as np

Data_type = ''
# # data type corresponds to the source file type
# data type corresponds to the source file type
# If the source file type is FWH data Data_type == 'FWH'
# If the source file type is from OUTPUT_POSTPROC Data_type == 'OUTPUT'
# If the source file type is from ISOSURFACE Data_type == 'EXTRACT'

#### THE DIRECTORY INFORMATION
# 1) The source transient solution directory
sol_dir = '/home/p/plavoie/denggua1/scratch/Bombardier_LES/SuperC_Small/RUN_Speedup/FWH_Airfoil_All/FWH_Data/'
# 2) The source transient solution name preamble string
sol_file = 'FWH_Airfoil_All'
# 2) The destination directory
dest_dir = '/home/p/plavoie/denggua1/scratch/Bombardier_LES/SuperC_Small/PostProc/DMD_Airfoil/'
# 3) The extracted and filtered raw data destination subdirectory
dest_file = 'FWH_Data'                          # Can personally Change
# 4) The mesh directory
mesh_Dir = '/home/p/plavoie/denggua1/scratch/Bombardier_LES/SuperC_Small/MESH/'
# 5) The source mesh file name in .h5 format
mesh_fileName = 'Bombardier_5AOA_Small.mesh.h5' 
# 6) The destination file name prefix
mesh_srcname = 'Airfoil_Surface_Mesh'           # Can personally Change
# 7) The destination DMD source data file
DMD_datafile = 'DMD_Data'                       # Can personally Change
# 8) The destination DMD output data directory
DMD_output_dir = 'DMD_Airfoil_Output'                   # Can personally Change
# 9) The destination DMD reconstructed field directory
DMD_reconstruct_dir = 'DMD_Reconstruction_Airfoil_Modes'        # Can personally Change

#### THE FILE SAMPLING OPTION
Option = 1
# Option 1 = Only take n files (maxfile), Option 2 = Skip every n files either until maxfiles or end of directory
nskip = 1               # Number of files to skip
max_file = 3600        # Maximum number of files


## THE DATA INFORMATION
# The DMD variable of interest
var = 'pressure'
# The timestep size from run.params
tstep = 1.578e-8
# The data sampling interval from run.params
samp_intl = 200
# Data recording interval
if Option == 1:
    intl = samp_intl*1
elif Option == 2:
    intl = samp_intl*nskip
# The DMD timestep size
dt = tstep*float(intl)
# The DMD frequency of interest
freq_start = 500
freq_end = 15000
df = 250
N = (freq_end-freq_start)/df+1
freq_interest = np.linspace(freq_start, freq_end, int(N), endpoint=True)
# The sweep range above and below the above frequency of interest

    