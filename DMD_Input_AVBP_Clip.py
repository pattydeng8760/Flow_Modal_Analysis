####################################################################################
# Applying Dynamic Mode Analysis to transient data
# Input Parameters from AVBP Post Extract dataset from Full solutions (CutPlane_Extact.py)
# Author: Patrick Deng
####################################################################################
import os
import numpy as np
#### THE DIRECTORY INFORMATION
# data type corresponds to the source file type
# If the source file type is FWH data Data_type == 'FWH'
# If the source file type is from OUTPUT_POSTPROC Data_type == 'OUTPUT'
# If the source file type is from ISOSURFACE Data_type == 'EXTRACT'
# If the source file type is from Clip Data_type == 'CLIP'
Data_type = 'CLIP'
# The DMD variable of interest
var = 'vort_x'
# 1) The source transient solution directory
sol_dir = '/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/Cut_Clip_New/'
# 2) The source transient solution name preamble string
sol_file = 'B_10AOA'
# 2) The destination directory
dest_dir = '/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/PostProc/DMD_Airfoil/'
# 3) The extracted and filtered raw data destination subdirectory                       
dest_file = 'Cut_Clip'
# 4) The mesh directory
mesh_Dir = '/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/MESH_ZONE_Apr24/'
# 5) The source mesh file name in .h5 format
mesh_fileName = 'Bombardier_10AOA_Combine_Apr24.mesh.h5'
# 6) The destination file name prefix
mesh_srcname =  dest_file +'_Mesh'          # Can personally Change
# 7) The destination DMD source data file
DMD_datafile = 'DMD_Data_' + dest_file + '_' + var                        # Can personally Change
# 8) The destination DMD output data directory
DMD_output_dir = 'DMD_Output_'+ var + '_'+ dest_file                # Can personally Change
# 9) The destination DMD reconstructed field directory
DMD_reconstruct_dir = 'DMD_Reconstruction_' + dest_file +  '_Modes'      # Can personally Change


#### THE FILE SAMPLING OPTION
option = 2
# option 1 = Only take n files (maxfile), option 2 = Skip every n files either until maxfiles or end of directory
nskip = 2               # Number of files to skip
max_file = 150         # Maximum number of files
reload = False
triple_decomp = False   # Perform the triple decomposition of the flow field
vortex = 'PV'                   # Adjust the vortex identification method to center on the specified vortex, not applicable for not triple decomposition
## THE DATA INFORMATION
# The timestep size from run.params
tstep = 1.878e-8
# The data sampling interval from run.params
samp_intl = 2000
# Data recording interval
if option == 1:
    intl = samp_intl*1
elif option == 2:
    intl = samp_intl*nskip
# The DMD timestep size
dt = tstep*float(intl)
# The DMD frequency of interest
freq_start = 500
freq_end = 6000
df = 250
N = (freq_end-freq_start)/(2*df)+1
freq_interest = np.linspace(freq_start, freq_end, int(N), endpoint=True)
#freq_interest = [10000,11000,10500,13000,14000,15000,16000,17000,18000,19000,20000]
# The sweep range above and below the above frequency of interest
