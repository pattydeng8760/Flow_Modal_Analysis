####################################################################################
# Applying Dynamic Mode Analysis to transient data
# Input Parameters from AVBP Post Extract dataset from Full solutions (CutPlane_Extact.py)
# Author: Patrick Deng
####################################################################################
import os
import numpy as np
import sys
#### THE DIRECTORY INFORMATION
# data type corresponds to the source file type
# If the source file type is FWH data Data_type == 'FWH'
# If the source file type is from OUTPUT_POSTPROC Data_type == 'OUTPUT'
# If the source file type is from ISOSURFACE Data_type == 'EXTRACT'
# If the source file type is from Clip Data_type == 'CLIP'
Data_type = 'EXTRACT'
# 1) The DMD variable of interest
var = 'pressure'
# 2) The source transient solution directory
sol_dir = sys.argv[1]
# 3) The source transient solution name preamble string
sol_file = 'B_10AOA'
# 4) The destination directory
dest_dir = os.getcwd()
# 5) The extracted and filtered raw data destination subdirectory, this is also the cut location                      
dest_file = os.path.basename(os.path.normpath(sol_dir))
# 6) The mesh directory
mesh_Dir = '/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/MESH_ZONE_Apr24/'
#mesh_Dir = '/project/p/plavoie/denggua1/BBDB_10AOA/MESH_ZONE_NEW/'
# 7) The source mesh file name in .h5 format
mesh_fileName = 'Bombardier_10AOA_Combine_Apr24.mesh.h5'
# 8) The destination DMD output data directory
DMD_output_dir = 'DMD_Output_'+ var + '_'+ dest_file                # Can personally Change

#### THE FILE SAMPLING OPTION
option = 2
extract_only = True if option == 3 else False
# option 1 = Only take n files (maxfile), option 2 = Skip every n files either until maxfiles or end of directory, option 3 = Take all files
nskip = 1               # Number of files to skip
nstart = 100
maxfile = 400        # Maximum number of files
reload_file=True         # Re-copy the files from the source directory
skip=False         # Re-load the raw data from the copied source files/skip the copying process
noise_lvl = 1e-4
triple_decomp = False   # Perform the triple decomposition of the flow field
vortex = 'PV'
## THE DATA INFORMATION
# The timestep size from run.params
tstep = 1.67e-8
# The data sampling interval from run.params
samp_intl = 2000
# Data recording interval
if option != 2:
    intl = samp_intl*1
elif option == 2:
    intl = samp_intl*nskip
# The DMD timestep size
dt = tstep*float(intl)
# The DMD frequency of interest
freq_start = 500
freq_end = 7000
df = 250
N = (freq_end-freq_start)/(2*df)+1
freq_interest = np.linspace(freq_start, freq_end, int(N), endpoint=True)
#freq_interest = [10000,11000,10500,13000,14000,15000,16000,17000,18000,19000,20000]
# The sweep range above and below the above frequency of interest
