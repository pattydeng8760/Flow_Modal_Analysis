####################################################################################
# Applying Dynamic Mode Analysis to transient data
# The Main DMD Function to be called
# Author: Patrick Deng
####################################################################################
# Importing the required modules
import os
import sys
# the subroutines
from DMD_Functions import print, Extract_Files, Extract_Surface, Extract_Data, DMD_Compute, Reconstruct_DMD_Field, timer
# the input parameters file
#from DMD_Input_AVBP_FWH import *
#from DMD_Input_Pressure_Hessian import *
from DMD_Input_AVBP_Extract import *
#from DMD_Input_AVBP_Extract_TripleDecomp import *
#from DMD_Input_AVBP_Clip import *
#### THE OUTPUT LOG FILE
sys.stdout = open(os.path.join('log_'+DMD_output_dir.replace('_Output', '')+'.txt'), "w", buffering=1)

text = 'Starting the DMD Program'
print(f'\n{text:.^120}\n') 

@timer
def main():
    output, DMD_output = Extract_Files(sol_dir,sol_file,dest_dir,dest_file,DMD_output_dir,option=option,maxfile=maxfile,nskip=nskip,reload=reload_file,skip=skip)
    DMD_mesh, mesh_node = Extract_Surface(mesh_Dir, mesh_fileName, DMD_output, Data_type, dest_file)
    ntime, DMD_source = Extract_Data(sol_dir, DMD_output, mesh_node, dest_file, Data_type, var, option=option,nstart=nstart,maxfile=maxfile,nskip=nskip,fluc=False, load_existing=skip)
    if extract_only: 
        return
    else:
        peak_modes, DMD_result_file, DMD_spectra_file, DMD_modes_file = DMD_Compute(DMD_mesh, DMD_source, DMD_output, ntime, var, dt, df, freq_interest, Data_type,\
            noise_lvl=noise_lvl, triple_decomp=triple_decomp, vortex=vortex)
        Reconstruct_DMD_Field(DMD_modes_file, DMD_spectra_file, peak_modes, DMD_output, var, dest_file, dtype='float32', nperiod=3)

if __name__ ==  '__main__':
    main()

text = 'DMD Program Complete!'
print(f'\n{text:.^120}\n') 