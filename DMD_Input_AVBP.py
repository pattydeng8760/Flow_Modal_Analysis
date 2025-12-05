####################################################################################
# Applying Dynamic Mode Analysis to transient data
# Input Parameters from AVBP OUTPUT_POSTPROC
# Author: Patrick Deng
####################################################################################
import os
import numpy as np

Data_type = 'OUTPUT'
# data type corresponds to the source file type
# If the source file type is FWH data Data_type == 'FWH'
# If the source file type is from OUTPUT_POSTPROC Data_type == 'OUTPUT'


#### THE DIRECTORY INFORMATION
# 1) The source transient solution directory
sol_dir = '/home/p/plavoie/denggua1/scratch/Bombardier_LES/SuperC_Small/RUN_Speedup/OUTPUT_POSTPROC/Z_cut_combine/'
# 2) The source transient solution name preamble string
#sol_file = 'x_095chord_cut'
sol_file = 'z_tip_onein_cut'
# 2) The destination directory
dest_dir = '/home/p/plavoie/denggua1/scratch/Bombardier_LES/SuperC_Small/PostProc/DMD_Airfoil/'
# 3) The extracted and filtered raw data destination subdirectory
#dest_file = 'X_095_onein'                          # Can personally Change
dest_file = 'x_tip_onein' 
# 4) The mesh directory
mesh_Dir = sol_dir
# 5) The source mesh file name in .h5 format
#mesh_fileName = 'x_095chord_mesh-2400.h5' 
mesh_fileName = 'z_tip_onein_mesh-2400.h5' 
# 6) The destination file name prefix
mesh_srcname = 'z_tip_onein_Mesh'           # Can personally Change
# 7) The destination DMD source data file
DMD_datafile = 'DMD_Data'                       # Can personally Change
# 8) The destination DMD output data directory
DMD_output_dir = 'DMD_Output_z_tip_onein'                   # Can personally Change
# 9) The destination DMD reconstructed field directory
DMD_reconstruct_dir = 'DMD_Reconstruction_z_tip_onein_Modes'        # Can personally Change

#### THE FILE SAMPLING OPTION
Option = 1
# Option 1 = Only take n files (maxfile), Option 2 = Skip every n files either until maxfiles or end of directory
nskip = 1               # Number of files to skip
max_file = 150         # Maximum number of files


## THE DATA INFORMATION
# The DMD variable of interest
var = 'dilatation'
# The timestep size from run.params
tstep = 1.578e-8
# The data sampling interval from run.params
samp_intl = 2000
# Data recording interval
if Option == 1:
    intl = samp_intl*1
elif Option == 2:
    intl = samp_intl*nskip
# The DMD timestep size
dt = tstep*float(intl)
# The DMD frequency of interest
freq_start = 500
freq_end = 7000
df = 250
N = (freq_end-freq_start)/(2*df)+1
freq_interest = np.linspace(freq_start, freq_end, int(N), endpoint=True)
#freq_interest = [10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000]
# The sweep range above and below the above frequency of interest
