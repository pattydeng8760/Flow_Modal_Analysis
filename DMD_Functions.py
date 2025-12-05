####################################################################################
# Applying Dynamic Mode Analysis to transient data
# The subfile with all required functions
# Functions to be called in DMD.py
# Author: Patrick Deng
# DO NOT MODIFY THIS FILE
####################################################################################
# The required modules
import os
import numpy as np
import sys
import time
from antares import *
import h5py
import copy
import glob
import builtins
import re
import shutil
import matplotlib.pyplot as plt

# Printing any on-screen print functions into a log file
def print(text:str,**kwargs):
    """ print function to print to the screen and to a log file
    """
    builtins.print(text,**kwargs)
    os.fsync(sys.stdout)

def timer(func):
    """ Decorator to time the function func to track the time taken for the function to run"""
    def inner(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        elapsed = end - start
        print('The total compute time is: {0:1.0f} min, {1:1.0f} s'.format(np.floor(elapsed/60), np.mod(elapsed,60))) 
        return elapsed
    return inner

# The cut location class
class constants:
    """ Constants for the cut location"""
    def __init__(self, tip_gap, span):
        # Midspan
        self.z_mid_span = (tip_gap + span/2)
        # 2 inch from tip
        self.z_2inch_tip = (tip_gap - 0.0508)
        # 1 inch from tip
        self.z_1inch_tip = (tip_gap - 0.0254)
        # 0.25 inch from tip
        self.z_025inch_tip = (tip_gap - 0.00635)
        self.z_5mm_tip = (tip_gap - 0.005)
        self.z_25mm_tip = (tip_gap - 0.025)
        self.z_tip_gap = tip_gap
    
# Extracting files in the source directory into the present DMD directory to be conditioned.
def Extract_Files(sol_dirName:str, sol_file:str, dest_dir:str, dest_File:str, DMD_output_dir:str\
    ,option:int=1, maxfile:int=500, nskip:int=1, reload:bool=True, skip:bool=False):
    """Extracting the Transient Files (post cut-extraction) located in the source directory into the present directorry
    This is to avoid overwriting the original files in the source directory
    Also initiating the output directory to store intermediate files.
    Args:
        sol_dirName (str): The source directory where the transient files are located
        sol_file (str): The base name of the transient files
        dest_dir (str): The destination directory where the transient files are to be stored
        dest_File (str): The base name of the transient files to be stored
        DMD_output_dir (str): The destination directory where the DMD output data is to be stored
        option (int): The option to extract the files, 1 for n files, 2 for skip every n files
        maxfile (int): The maximum number of files to extract
        nskip (int): The number of files to skip
        reload (bool): The option to reload the files if they are already extracted
        skip (bool): The option to skip the extraction if the data is already extracted
    Returns:
        ntime (int): The total number of timesteps extracted
        output (str): The destination directory where the transient files are stored
    """
    text = 'Extracting the Transient Files'
    print(f'\n{text:.^80}\n')  
    # The Output directory for the transient data
    output = os.path.join(dest_dir,dest_File)
    
    # Making the output directory for the DMD data
    DMD_output = os.path.join(dest_dir, DMD_output_dir)
    if os.path.exists(DMD_output) and reload==True:
        shutil.rmtree(DMD_output)
    os.makedirs(DMD_output, exist_ok = True)
    if skip:
        print('The extraction is skipped...')
        arr = sorted(glob.glob('{0:s}/{1:s}*.h5'.format(sol_dirName,sol_file)))
        n_files = len(arr)
    else:
        # An array with all the files in the source directory
        arr,arr_dest = [], []
        arr += sorted(glob.glob('{0:s}/{1:s}*.h5'.format(sol_dirName,sol_file)))
        n_files = len(arr)
        if os.path.exists(output):
            print('There are {0:d} files in the source directory'.format(len(arr)))
            arr_dest += sorted(glob.glob('{0:s}/*.h5'.format(output)))
            print('There are {0:d} files in the destination directory'.format(len(arr_dest)))
            print('After extraction there should be {0:d} files in the destination directory'.format(int(np.floor(len(arr))/nskip)))
        # Checking if all the files from the source directory are already extracted
        if ((int(len(arr_dest)) == int(np.floor(len(arr))/nskip) and os.path.exists(output)) or (int(len(arr_dest)) >= maxfile)) or (reload == False and os.path.exists(output)):
            text = 'All the files are already extracted to destination directory'
            print(f'{text}')
        else:
            print('----> Extracting the transient files from the source directory...')
            # Removing the directory if exists
            if os.path.exists(output):
                shutil.rmtree(output)
            os.makedirs(output, exist_ok = False)
            if option != 2:         # Only take n files
                for i in range(0,np.int(len(arr))):
                    source = os.path.join(sol_dirName,arr[i])
                    destination = output
                    shutil.copy(source,destination)
                    if np.mod(i,10)==0:
                        print(source)
                arr2 = os.listdir(output)
                arr2 = list(arr2)
                arr2.sort()
            elif option == 2:         # Skip every n files
                arr = []
                arr += sorted(glob.glob('{0:s}/{1:s}*.h5'.format(sol_dirName,sol_file)))
                for i in range(0,np.min([np.int(len(arr)/nskip),np.int(maxfile)])):
                    source = os.path.join(sol_dirName,arr[np.int(i*nskip)])
                    destination = output
                    shutil.copy(source,destination)
                    if np.mod(i,10)==0:
                        print(source)
        list_files=[]
        list_files+=sorted(glob.glob('{0:s}/*.h5'.format(output)))
        ntime = np.shape(list_files)[0]
        text = 'The total number of timesteps extracted is: {0:s}'.format(str(ntime))
        print(f'{text}')  
        text = 'The post sampled transient source data is stored in: {0:s}'.format(output)
        print(f'\n{text}')  
    text = 'File Extraction Complete'
    print(f'\n{text:.^80}\n')  
    return output, DMD_output

# Mappign the cut selection location
def map_cut(data_type:str,cut_style:str,tip_gap:float,span:float,AoA:int):
    """Function to map the cut selection to the location of the cut plane given the tip gap and span constnats
    Args:
        data_type (str): The style of the cut, either 'plane' or 'cylinder'
        cut_style (str): The location of the cut plane, defined either explicitly or with a %_TE location designating the distance from the trailing edge
        tip_gap (float): the tip gap size
        span (float): the span size
        AoA (int): the angle of attack
    Returns:
        origin (list): The origin of the cut plane
        normal (list): The normal of the cut plane
    """
    z_loc = constants(tip_gap,span)
    if data_type == 'EXTRACT':
        if cut_style.find("midspan") != -1:
            origin = [1.225,0.,z_loc.z_mid_span]
            normal = [0.,0.,1.]
        elif cut_style.find("2inch_tip") != -1:
            origin = [1.225,0.,z_loc.z_2inch_tip]
            normal = [0.,0.,1.]
        elif cut_style.find("1inch_tip") != -1:  
            origin= [1.225,0.,z_loc.z_1inch_tip]
            normal = [0.,0.,1.]
        elif cut_style.find("025inch_tip") != -1:  
            origin = [1.225,0.,z_loc.z_025inch_tip]
            normal = [0.,0.,1.]
        elif cut_style.find("25mm_tip") != -1:  
            origin= [1.225,0.,z_loc.z_25mm_tip]
            normal = [0.,0.,1.]
        elif cut_style.find("5mm_tip") != -1:  
            origin= [1.225,0.,z_loc.z_5mm_tip]
            normal = [0.,0.,1.]
        elif cut_style.find('PIV1') != -1:
            x,y,z = 1.42222035, 0, z_loc.z_mid_span
            origin = [x,y,z]
            normal  = [1,0,0]
        elif cut_style.find('PIV2') != -1:
            x,y,z = 1.48172998, 0, z_loc.z_mid_span
            origin = [x,y,z]
            normal  = [1,0,0]
        elif cut_style.find('PIV3') != -1:
            x,y,z = 1.5641908, 0, z_loc.z_mid_span
            origin= [x,y,z]
            normal = [1,0,0]
        elif cut_style.find("TE") != -1:
            Loc = float(re.findall(r"\d+", cut_style)[0])/100
            PIV = 1.25 + np.array(Loc)*0.3048*np.cos(AoA*np.pi/180)
            origin =  [PIV,0.,z_loc.z_mid_span]
            normal = [1.,0.,0.]
    elif data_type == 'CLIP':
        origin = [1.42222035,0.,z_loc.z_tip_gap]
        normal = [1.,0.,0.]
    print('The selected cut is of style: {0}'.format(data_type))
    print('The selected cut origin is: {0}'.format(origin))
    print('The selected cut normal is: {0}'.format(normal))
    return origin,normal


# Extracting the airfoil surface mesh for the base cut_locationof DMD
def Extract_Surface(mesh_fileDir:str,mesh_fileName:str,DMD_output:str, data_type:str,cut_location:str,\
    span:float=-0.2286,tip_gap:float=-0.1034,AoA:int=10,reload:bool=False):
    """ Extracting the surface mesh for the DMD computation
    Args:
        mesh_fileDir (str): The directory where the avbp mesh file is located (post zone merge)
        mesh_fileName (str): The name of the AVBP mesh file of the whole domain
        DMD_output(str): The destination directory where the extracted mesh is to be stored
        data_type (str): The type of data to be extracted, either 'FWH', 'OUTPUT', 'EXTRACT', or 'CLIP'
        cut_location (str): The location of the cut plane, defined either explicitly via PIV planes
            or with a %_TE location designating the distance from the trailing edge
        span (float): the span size
        tip_gap (float): the tip gap size
        AoA (int): the angle of attack
        reload (bool): The option to reload the files if they are already extracted
    Returns:
        dest_mesh (str): The destination directory where the extracted mesh is stored
        nodes (int): The number of nodes on the surface
    """
    text = 'Extracting Mesh'
    print(f'\n{text:.^80}\n')
    # The name of the mesh file is the cut location name
    mesh_name = cut_location + '_Mesh'  
    DMD_mesh = os.path.join(DMD_output,mesh_name+'.h5')
    if os.path.exists(DMD_mesh) == True and reload == False:
        print('----> LES Mesh already extracted at {0:s}'.format(DMD_mesh))
        # Loading the mesh
        r = Reader('hdf_antares')
        r['filename'] = DMD_mesh
        mesh = r.read()
        mesh.show()
        nodes = mesh[0][0]['x'].shape[0]
    elif data_type == 'FWH':
        ## Loading the mesh
        text = '----> Extracting the Surface Mesh from FWH Surface'
        print(f'{text}')  
        ## Loading the Main LES Mesh File
        r = Reader('hdf_avbp')
        r['filename'] = os.path.join(mesh_fileDir,mesh_fileName)
        base  = r.read() # b is the Base object of the Antares API
        airfoil_base = Family()
        airfoil_base['Airfoil_Surface'] = base.families['Patches']['Airfoil_Surface']
        airfoil_base['Airfoil_Trailing_Edge'] = base.families['Patches']['Airfoil_Trailing_Edge']
        airfoil_base['Airfoil_Side_LE'] = base.families['Patches']['Airfoil_Side_LE']
        airfoil_base['Airfoil_Side_Mid'] = base.families['Patches']['Airfoil_Side_Mid']
        airfoil_base['Airfoil_Side_TE']  = base.families['Patches']['Airfoil_Side_TE']
        base.families['SKIN'] = airfoil_base
        skin_base = base[base.families['SKIN']]
        text = '----> The Extracted Base Objects'
        print(f'\n{text}')   
        print(skin_base)
        ## Merging the extracted base objects to the same zone
        text = 'Merging the Base Objects'
        print(f'\n{text}')   
        myt = Treatment('merge')
        myt['base'] = skin_base
        myt['duplicates_detection'] = False
        myt['tolerance_decimals'] = 13
        # Writing the extraced mesh
        text = '----> Writing the Mesh File'
        print(f'\n{text}')   
        merged = myt.execute()
        writer = Writer('hdf_antares')
        writer['base'] = merged
        writer['filename'] = os.path.join(DMD_output,mesh_name)
        writer.dump()
        # The data for the original mesh
        text = '----> The Original Mesh Surface'
        print(f'\n{text}')   
        skin_base.show()
        # The data for the extracted and merged mesh
        text = '----> The Post Extraced Mesh Surface'
        print(f'\n{text}')   
        merged.show()
        nodes = merged['0000'][0].shape[0]
        print('\nThe number of nodes on the surface is: {0:d}'.format(nodes))
        DMD_mesh = os.path.join(DMD_output,mesh_name+'.h5')
        print('\nThe Extracted surface mesh is saved in {0:s}'.format(DMD_mesh))
    elif data_type == 'OUTPUT':
        # This is not suggested as requires extraction at run time, better to save the full soltuion and run Cutplanes post processing
        text = '----> Extracting the output mesh from AVBP OUTPUT POSTPROC'
        print(f'{text}')   
        r = Reader('hdf_avbp')
        mesh = os.path.join(mesh_fileDir,mesh_fileName)
        DMD_mesh = os.path.join(DMD_output,mesh_name+'.h5')
        r['filename'] = mesh
        base  = r.read() # b is the Base object of the Antares API
        shutil.copyfile(mesh, DMD_mesh)
        nodes = base['0000'][0].shape[0]
        text = '----> The AVBP OUTPUT POSTPROC Database Mesh Surface'
        base.show()
        print('\nThe number of nodes on the surface is: {0:d}'.format(nodes))
        print('\nThe copied mesh is saved in {0:s}'.format(DMD_mesh))
    elif data_type == 'EXTRACT' or data_type == 'CLIP':
        # Loading the mesh
        text = '----> Extracting the Cut Mesh'
        print(f'{text}')  
        # Loading the Main LES Mesh File
        r = Reader('hdf_avbp')
        r['filename'] = os.path.join(mesh_fileDir,mesh_fileName)
        base  = r.read() # b is the Base object of the Antares API
        # Compute the origin and normal of the cut plane location
        if data_type == 'CLIP':
            t= Treatment('clip')
            t['base'] = base
            cut_style = 'cylinder'
            t['type'] = cut_style
            t['axis'] = 'x'
            t['radius'] = 0.1
        elif data_type == 'EXTRACT':
            t= Treatment('cut')
            t['base'] = base
            t['type'] = 'plane'
            cut_style = cut_location
        origin, normal = map_cut(data_type,cut_style,tip_gap,span,AoA)
        t['origin'] = origin
        t['normal'] = normal
        inter = t.execute()     
        # Writing the extraced mesh
        text = '----> Writing the Mesh File'
        print(f'\n{text}')   
        writer = Writer('hdf_antares')
        writer['base'] = inter
        writer['filename'] = os.path.join(DMD_output,mesh_name)
        writer.dump()
        # The data for the extracted and merged mesh
        text = '----> The Post Extraced Mesh Surface'
        print(f'\n{text}')   
        inter.show()
        nodes = inter['0000'][0].shape[0]
        print('\nThe number of nodes on the surface is: {0:d}'.format(nodes))
        DMD_mesh = os.path.join(DMD_output,mesh_name+'.h5')
        print('\nThe Extracted surface mesh is saved in: {0:s}'.format(DMD_mesh))
    text = 'Mesh Extraction Complete!'
    print(f'\n{text:.^80}\n')  
    return DMD_mesh, nodes

# Extracting data in the post sampled source directory into a unified file to be read by the DMD script
def Extract_Data(sol_dir:str,DMD_output:str,nodes:float,cut_location:str,data_type:str,var:str, \
    option:int=1,nstart:int=0,nskip:int=1,maxfile:int=1000,fluc:bool=False, load_existing:bool=False):
    """ Extracting the source data into a unified hdf5 file to be read by the DMD script and numpy file for backup
    Args:
        sol_dir (str): The source directory where the transient files are located 
        DMD_output (str): The destination directory where the extracted data is to be stored
        nodes (int): The number of nodes on the surface from the meah file (to compare with the data to ensure consistency)
        data_type (str): The type of data to be extracted, either 'FWH', 'OUTPUT', 'EXTRACT', or 'CLIP'
        cut_location (str): The location of the cut plane, defined either explicitly via PIV planes
        var (str): The variable of interest to be extracted, default is 'pressure'
        fluc (bool): The option to extract the fluctuating part of the data for analysis, default is False
        load_existing (bool): The option to load the existing data from the copied directory from Extract_Files, default is False, if True, the data is reloaded from the saved hdf5 file
    Returns:
        ntime (int): The total number of timesteps extracted
        DMD_source (str): The destination path where the extracted data is stored as hdf5 file to be read by the DMD
    """
    text = 'Beginning Data Extraction'
    print(f'\n{text:.^80}\n')  
    # The name of the output data file
    data_files=sorted(glob.glob('{0:s}/DMD_Data*.hdf5'.format(DMD_output)))
    'DMD_Data_' + var + '_'+ cut_location                      # Can personally Change
    if not load_existing:
        print('---->loading the raw data from the copied directory')
        # The directory information
        l_orig=sorted(glob.glob('{0:s}/*.h5'.format(sol_dir)))
        l_count = np.linspace(0,len(l_orig)-1,len(l_orig),dtype=int)
        nb_source = len(l_orig)
        maxfile = np.min([maxfile,nb_source])
        if option == 1:
            print('Option 1: Extracting the first {0:d} files sequentially:'.format(maxfile))
            l,l_count=l_orig[nstart:nstart+maxfile],l_count[nstart:nstart+maxfile]
        elif option == 2:
            print('Option 2: Extracting every {0:d} files sequentially until {1:d}:'.format(nskip, maxfile))
            l_orig,l_count=l_orig[nstart::nskip],l_count[nstart::nskip]
            l,l_count=l_orig[:maxfile],l_count[:maxfile]
        elif option == 3:
            print('Option 3: Extracting all files:')
            l_orig,l_count=l_orig,l_count
            l,l_count=l_orig,l_count
        nb_files=len(l)
        print('The number of extracted files is: {0:d}.'.format(nb_files))
        print('The variable of interest is: {0:s}.'.format(var))
        # Number of nodal points from the surface
        nb_points = nodes
        # Pre allocating space
        data=np.zeros((nb_points,nb_files),dtype=np.float64)
        # Running the loop to extract the pressure data
        for it,filename in enumerate(l):
            if data_type == 'FWH':
                f=h5py.File(filename,'r')
                press=f['frame_data/pressure'][()]
                #print(len(press))
                data[:,it]=press
                f.close()
            elif data_type == 'OUTPUT':
                r = Reader('hdf_avbp')
                r['filename'] = filename
                base  = r.read()
                press = base[0][0][(var, 'node')]
                data[:,it]=press
            elif data_type == 'EXTRACT' or data_type == 'CLIP':
                r = Reader('hdf_antares')
                r['filename'] = filename
                base  = r.read()
                press = base[0][0][(var, 'node')]
                data[:,it]=press
            # Check for the dimensions
            if it == 0:
                if nb_points != press.shape[0]:
                    print('The number of nodes in the datafile is: {0:d}'.format(press.shape[0]))
                    print('The number of nodes in the meshfile is: {0:d}'.format(nb_points))
                    raise ValueError('The number of nodes on the surface is not consistent with the mesh file')
                else:
                    print('The number of nodes in the datafile is: {0:d} equal to that of the mesh file: {1:d}\n'.format(press.shape[0],nb_points))
            if np.mod(it,10) == 0 or it==nb_files-1:
                print('Iteration {0:03d} of {1:03d}: Extracting file {2:03d} of {3:03d} from the original'.format(it, nb_files, l_count[it], nb_source))
        data_mean = np.mean(data,axis=1)
        data_mean = data_mean.reshape((nb_points,1))
        if fluc :
            print('Extracting the fluctuating part of the data')
            data = data - data_mean  # Subtracting the mean value to get the fluctuating part
        print('\nThe maximum and mimumum {0:s} in the instantaneous field is: {1:.4f}, {2:.4f}'.format(var,np.max(data),np.min(data)))
        # The number of timesteps
        ntime = data.shape[1]
        print('The number of timesteps extracted is: {0:d}'.format(ntime))
        # Saving the output as hdf5 file 
        if option != 3: 
            DMD_source = os.path.join(DMD_output,'DMD_Data_{0:s}_{1:s}_{2:04d}inst.hdf5'.format(var,cut_location,ntime))
        else: 
            DMD_source = os.path.join(DMD_output,'DMD_Data_{0:s}_{1:s}_max_inst.hdf5'.format(var,cut_location))
        f=h5py.File(DMD_source,'w')
        data_float32 = data.astype(np.float64)
        f[var]=data_float32
        f.close()
    elif load_existing and os.path.exists(os.path.join(DMD_output,'DMD_Data_{0:s}_{1:s}_max_inst.hdf5'.format(var,cut_location))):
        print('---->Loading the data from the saved .hdf5 file with the full data')
        with h5py.File(os.path.join(DMD_output,'DMD_Data_{0:s}_{1:s}_max_inst.hdf5'.format(var,cut_location)),'r') as f:
            data = f[var][()]
            if option == 1:
                data = data[:,nstart:nstart+maxfile]
            elif option == 2:
                data = data[:,nstart::nskip]
                data = data[:,:maxfile]
            ntime = data.shape[1]
        DMD_source = os.path.join(DMD_output,'DMD_Data_{0:s}_{1:s}_{2:04d}inst.hdf5'.format(var,cut_location,ntime))
        print('The number of timesteps extracted is: {0:d}'.format(ntime))
        print('The maximum and mimumum {0:s} in the instantaneous field is: {1:.4f}, {2:.4f}'.format(var,np.max(data),np.min(data)))
        f=h5py.File(DMD_source,'w')
        data_float32 = data.astype(np.float64)
        f[var]=data_float32
        f.close()
    # Removing the transient data from the directory
    # if os.path.exists(sol_dir):
    #     shutil.rmtree(sol_dir, ignore_errors=False, onerror=None)
    #     print('The source directory {0:s} is removed.'.format(sol_dir))
    # else: 
    #     print('---->Loading the data from the saved .hdf5 file')
    #     # Loading the data from the hdf5 file
    #     DMD_source = data_files[-1]
    #     with h5py.File(DMD_source,'r') as f:
    #         data = f[var][()]
    #     ntime = data.shape[1]
    #     print('The number of timesteps extracted is: {0:d}'.format(ntime))
    #     print('The maximum and mimumum {0:s} in the instantaneous field is: {1:.4f}, {2:.4f}'.format(var,np.max(data),np.min(data)))
    print('The destination DMD data file is: {0:s}.'.format(DMD_source))
    text = 'Data Extraction Complete!'
    print(f'\n{text:.^80}\n')  
    return ntime, DMD_source

# the Main DMD function compute
def DMD_Compute(DMD_mesh:str, DMD_source:str, DMD_output:str,\
    ntime, var, dt, df, freq_interest, data_type, \
        noise_lvl:float=1e-6, triple_decomp:bool=False, vortex:str='PV'):
    """ The main DMD computation function
    Args:
        DMD_mesh (str, optional): _description_. Defaults to './DMD_mesh'.
        DMD_source (str, optional): _description_. Defaults to './DMD_Data.hdf5'.
        DMD_output (str): The destination directory where the DMD output is to be stored
        ntime: the number of timesteps to be considered in the DMD computation
        var: the variable of interest to apply the DMD on
        dt: the timestep size
        df: the frequency sweep bin range
        freq_interest: table of frequencies of interest, the sweep will be performed at these frequencies +/- df
        data_type: the type of data to be extracted, either 'FWH', 'OUTPUT', 'EXTRACT', or 'CLIP'
        noise_lvl (float): The noise level to be considered in the DMD computation, default is 1e-6
        triple_decomp (bool): The option to perform the triple decomposition of the flow field, default is False
        vortex (str): The vortex to be considered in the triple decomposition, default is 'PV'
    Returns:
        peak_modes (antares base): The peak modes of the DMD computation
        DMD_result (str): The destination path where the DMD result is stored
        DMD_spectra (str): The destination path where the DMD spectra is stored
        DMD_modes (str): The destination path where the DMD modes are stored
    """
    text = 'Performing DMD'
    print(f'\n{text:.^80}\n')  
    # Initializing the data matrix
    data = {}
    # Number of instants to account for in DMD, Must match with data file
    ninst = ntime
    # Add frequencies of interest 
    freqs = freq_interest
    print('The frequency range of interest is: {0:s} Hz.'.format(str(freq_interest)))
    print('The DMD sampling frequency is: {0:1.1f} Hz.'.format(1/dt))
    print('The DMD sweep range is +/-: {0:1.0f} Hz at each frequency range.'.format(df))
    # Frequency interval: the algorithm will search for modes
    # with maximum amplitude for selected frequencies within given frequency interval:
    #################################################################################### 
    if data_type == 'FWH':
        #Open mesh
        print('\n----> Open Mesh file')
        print('The DMD Mesh file is located at {0:s} '.format(DMD_mesh))
        r = Reader('hdf_antares')
        r['filename'] = DMD_mesh
        mesh = r.read()
        # Compute the volume of cells (2D or 3D), needed for modal scaling in the algorithm
        mesh.compute_cell_volume()
        # Project cell volume values to nodes
        mesh.cell_to_node()
        # Leave only the coordinates and cells volumes stored in nodes
        mesh = mesh[:,:,[('x','node'),('y','node'),('z','node'),('cell_volume','node')]]
        mesh.show()
    elif data_type == 'OUTPUT':
        #Open mesh
        print('\n----> Open Mesh file')
        print('The DMD Mesh file is located at {0:s} '.format(DMD_mesh))
        r = Reader('hdf_avbp')
        r['filename'] = DMD_mesh
        mesh = r.read()
        # Compute the volume of cells (2D or 3D), needed for modal scaling in the algorithm
        mesh.compute_cell_volume()
        # Project cell volume values to nodes
        mesh.cell_to_node()
        # Leave only the coordinates and cells volumes stored in nodes
        mesh = mesh[:,:,[('x','node'),('y','node'),('z','node'),('cell_volume','node')]]
        mesh.show()
    elif data_type == 'EXTRACT' or data_type == 'CLIP':
        #Open mesh
        print('\n----> Open Mesh file')
        print('The DMD Mesh file is located at {0:s} '.format(DMD_mesh))
        r = Reader('hdf_antares')
        r['filename'] = DMD_mesh
        mesh = r.read()
        # Compute the volume of cells (2D or 3D), needed for modal scaling in the algorithm
        mesh.compute_cell_volume()
        # Project cell volume values to nodes
        mesh.cell_to_node()
        # Leave only the coordinates and cells volumes stored in nodes
        mesh = mesh[:,:,[('x','node'),('y','node'),('z','node'),('cell_volume','node')]]
        mesh.show()
    # Applying the Triple Decomposition if required
    if triple_decomp:
        print('\n----> Performing Triple Decomposition')
        u, u_mean, u_coherent, u_random = triple_decomposition_main(DMD_source, var, DMD_mesh, DMD_output.split('_')[-1], vortex)
        data[var] = u_coherent
    else: 
        # Open data for n instants (only variable to process, the mesh is stored separately)
        print('\n----> Open Conditioned Results file')
        print('The POD source file is located at {0:s} '.format(DMD_source))
        fin = h5py.File(DMD_source,'r')
        data[var] = fin['/{0:s}'.format(var)][:,0:ninst]
        fin.close()
        del fin
    # Create a new base to make DMD
    b = Base()
    b['0'] = Zone()
    b[0].shared['cell_volume'] = mesh[0][0]['cell_volume']
    print('\nInstant   isinf    isnan')
    for i in range(int(data[var].shape[1])):
        instant_name = 'snapshot_%s' % i
        b[0][instant_name] = Instant()
        b[0][instant_name][var] = data[var][:,i].flatten(order='F')
        if np.mod(i,10)==0:
            text=('{0:s}  {1:}  {2:} '.format(instant_name, np.isinf(np.array(b[0][instant_name][var])).any(), np.isnan(np.array(b[0][instant_name][var])).any()))
            print(f'{text}') 
        if np.isinf(np.array(b[0][instant_name][var])).any() or np.isnan(np.array(b[0][instant_name][var])).any():
            raise ValueError('The data contains NaN or Inf values')
    #################################################################################### 
    # Perform DMD
    print('\n---->Performing DMD Computation...')
    t = Treatment('dmd')
    t['base'] = b
    t['variables'] = ['cell_volume', var]
    # For the settings see Antares documentation
    t['noiselevel'] = noise_lvl
    t['type'] = 'mod/phi'
    t['memory_mode'] = False
    t['list_freq'] = freqs
    t['delta_freq'] = df
    t['time_step'] = dt
    res_dmd = t.execute()
    # moving the mesh and intermediate data files to the output directory
    print('DMD Computation Complete')

    ####################################################################################
    # Extracting the DMD Data
    text = '---->Extracting DMD Data'
    print(f'\n{text}') 
    print('The peak modes are:')
    print(res_dmd['modes'].attrs['peak_modes_{0:s}'.format(var)])
    peak_modes = res_dmd['modes'].attrs['peak_modes_{0:s}'.format(var)]
    w = Writer('hdf_antares')
    DMD_result = os.path.join(DMD_output,'DMD_{0:s}_{1:04d}inst'.format(var,ninst))
    print('The DMD result file is saved as {0:s} '.format(DMD_result))
    w['filename'] = DMD_result
    w['base'] = res_dmd
    w['dtype'] = 'float32'
    w.dump()
    
    # Extract the spectrum
    freq = res_dmd['spectrum']['spectrum']['frequency']
    ampl = res_dmd['spectrum']['spectrum']['amplification']
    mod = res_dmd['spectrum']['spectrum']['{0:s}_modulus'.format(var)]
    spectrum = np.concatenate((freq[:,np.newaxis],ampl[:,np.newaxis],mod[:,np.newaxis]),axis=-1)
    DMD_spectra = os.path.join(DMD_output,'DMD_{0:04d}inst_{1:s}_spectrum.npy'.format(ninst,var))
    print('The DMD spectra file is saved as {0:s} '.format(DMD_spectra))
    np.save(DMD_spectra,spectrum)
    # Write modes to mesh file for visual post-processing
    modes = len(res_dmd['modes'].keys())
    for mode in peak_modes:
        mesh[0][0]['{0:s}_{1:s}_mod'.format(var,mode)] = res_dmd['modes'][mode]['{0:s}_mod'.format(var)]
        mesh[0][0]['{0:s}_{1:s}_phi'.format(var,mode)] = res_dmd['modes'][mode]['{0:s}_phi'.format(var)]
    del mesh[0][0]['cell_volume']
    w = Writer('hdf_antares')
    DMD_modes = os.path.join(DMD_output,'DMD_{0:04d}inst_{1:s}_modes'.format(ninst,var))
    print('The DMD modes file is saved as {0:s} '.format(DMD_modes))
    w['filename'] = DMD_modes
    w['base'] = mesh
    w['dtype'] = 'float32'
    w.dump()
    text = 'DMD Complete!'
    print(f'\n{text:.^80}\n')
    return peak_modes, DMD_result, DMD_spectra, DMD_modes

def Plot_DMD_Spectra(f, amplitude, amplification, DMD_output, cut_location, var):
    
    # The amplitude plot
    plt.figure(1)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams.update({'font.size': 16})
    plt.semilogx(f,np.abs(amplitude),'-k')
    plt.xlabel(r'$f $ [Hz]')
    plt.ylabel(r'Amplitude')
    plt.xlim([100, 10000])
    fig = plt.gcf()
    fig.set_size_inches(8, 4, forward=True)
    fig.tight_layout()
    plt.savefig('{0:s}/DMD_Amplitude_{1:s}_{2:s}.png'.format(DMD_output,cut_location,var))
    
    plt.figure(2)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams.update({'font.size': 16})
    plt.semilogx(f[2:-1],amplification[2:-1],'-k')
    plt.xlabel(r'$f $ [Hz]')
    plt.ylabel(r'Amplification')
    plt.xlim([100, 10000])
    plt.ylim([np.min(amplification[2:-1]),np.max(amplification[2:-1])])
    fig = plt.gcf()
    fig.set_size_inches(8, 4, forward=True)
    fig.tight_layout()
    plt.savefig('{0:s}/DMD_Amplification_{1:s}_{2:s}.png'.format(DMD_output,cut_location,var))
    

def Reconstruct_DMD_Field(DMD_result, DMD_spectra, peak_modes, DMD_output, var, cut_location, dtype:str='float32', nperiod:int=2):
    """
    Reconstructing the DMD field for the selected modes
    NOTE! The animation is extracted not in real time, but relative to the modal period.
    For example, the mode at 10 kHz will have the period of 1e-4 s.
    For the animation, we take the period equal to 48 snapshots, the total number of snapshots is 3*48.
    That means that the mode will be reconstructed for 3 modal periods, which is 3*1e-4 s of the real time.
    Of course, you can set the period in seconds, then the time step of the animation will be period/nsteps.
    """
    # subdirectory for the DMD reconstruction results 
    dest_subdir = 'DMD_Reconstruction_'+cut_location+'_'+var+'_Modes'
    text = 'Reconstructing the Field'
    print(f'\n{text:.^80}\n')  
    print('\n----> Loading the DMD Result File')
    r = Reader('hdf_antares')
    file = os.path.join(DMD_result+'.h5')
    r['filename'] = file
    b = r.read()
    print('The loaded DMD result file: {0:s}'.format(file))
    print('The contents of this file:')
    print((b[0][0].keys()))
    
    print('\n----> Loading the DMD Spectra File')
    print('Saving the DMD spectra')
    data = np.load(DMD_spectra)
    f = data[:,0]
    amplitude = data[:,1]
    amplification = data[:,2]
    # Plotting the DMD spectra
    Plot_DMD_Spectra(f, amplitude, amplification, DMD_output, cut_location, var)
    
    print('\n----> Reconstructing the {0:s} Data...'.format(var))
    for i in range(0,np.shape(peak_modes)[0]):
        mode_int = peak_modes[i]
        mode_number = int(mode_int.replace('H',''))
        freq = int(f[mode_number])

        print('\nReconstructing Mode {0:s} corresponding to {1:1.0f} Hz'.format(peak_modes[i],freq))
        output = '{0:s}'.format(peak_modes[i])
        output_path = os.path.join(DMD_output,dest_subdir,output)
        os.makedirs(output_path, exist_ok = True)
        # Number of snapshots in one period of the mode
        period = 48
        # Number of snapshots to extract (better to make it proportional to period to have a smooth looped animation)
        nsteps = int(nperiod*period)
        #print(nsteps)
        # Name of the mode to extract
        mode_name = var +'_{0:s}'.format(peak_modes[i])
    
        for k in range(nsteps):
            #print("Step %2d"%(k))  
            animated_base = Base()
            animated_base['0'] = Zone()
            animated_base[0].shared.connectivity = b[0][0].connectivity
            animated_base[0].shared["x"] = b[0][0]["x"]
            animated_base[0].shared["y"] = b[0][0]["y"]
            animated_base[0].shared["z"] = b[0][0]["z"]
            animated_base[0][str(k)] = Instant()
            animated_base[0][str(k)][mode_name] = b[0][0][mode_name+'_mod']*np.cos(2.*np.pi*float(k)/float(period) + b[0][0][mode_name+'_phi']*np.pi/180.)
            animated_base[0][str(k)].attrs['Time'] = float(k)
            w = Writer('hdf_antares')
            w['filename'] = os.path.join(output_path,'DMD_{0:s}_{1:s}_anim_{2:03d}'.format(mode_name,str(freq),k))
            # w['coordinates'] = ['x','y','z']
            w['base'] = animated_base
            w['dtype'] = dtype
            w.dump()
        print('The DMD Result for Mode {0:s} is extracted.'.format(peak_modes[i]))
        print('The DMD Result saved at: {0:s}'.format(output_path))
    text = 'Reconstructing Field Complete!'
    print(f'\n{text:.^80}\n')  
    

#-------------------------------------------------------------------------------------------------------------------
# Functions for Triple Decomposition
#-------------------------------------------------------------------------------------------------------------------

def get_window_ranges(cut:str, vortex:str='PV'):
    """Function that contains the window ranges for the different cuts and vortex types for the triple decomposition. 
    Returns the window range for the specific cut and vortex type to track the coherent vortex movment.
    Args:
        cut (str): _description_
        vortex (str, optional): _description_. Defaults to 'PV'. Must be one of 'SV', 'PV', or 'TV'.
    Returns:
        window_ranges (dict): The window ranges for the specific cut and vortex type as 2 cartesian points corresponding to the lower left and upper right corners of the window.
    """
    assert vortex in ['SV', 'PV', 'TV'], "Invalid vortex type provided. Must be one of 'SV', 'PV', or 'TV'."
    windows = {
        'PIV1': {
            'vortex': {
                'SV': {
                    'LL': [0.01478, -0.10865],
                    'UR': [0.02538, -0.12091]
                },
                'PV': {
                    'LL': [-0.01917, -0.101929387],
                    'UR': [0.0064944, -0.0891985279]
                }
            }
        },
        'PIV2': {
            'vortex': {
                'SV': {
                    'LL': [-0.0018528727653921934, -0.11298976547803502],
                    'UR': [0.018302767560581263, -0.1335973259906]
                },
                'PV': {
                    'LL': [-0.013308410196503851, -0.08472011694431547],
                    'UR': [0.009774859757808086, -0.10123000538845912]
                },
            'TV': {  
                'LL': [-0.023754481588631, -0.0949114061073671],
                'UR': [-0.014531364896070058, -0.10168861340079643]
                }
            }
        },
        # Add other cuts as needed...
    }
    
    # Fetching the ranges for the specific cut and vortex type
    try:
        return windows[cut]['vortex'][vortex]
    except KeyError:
        raise ValueError(f"Invalid combination of cut '{cut}' and vortex '{vortex}' provided.")


def track_core_positions(vorticity_data:float, x:float, y:float, x_range:float, y_range:float):
    """ Function to track the core positions of the coherent vortex in the flow field.
    Args:
        vorticity_data (float): The extracted vorticity data from the flow field. Can be any data but typically vorticity.
        x (float): 1D array of the x coordinates of the mesh nodes. (in plane)
        y (float): 1D array of the y coordinates of the mesh nodes. (in plane)
        x_range (float): window range in the x direction to track the core position.
        y_range (float): window range in the y direction to track the core position.

    Returns:
        core_positions (list): List of core positions at each time step.
    """
    # Initialize the list to store the core positions
    core_positions = []
    # Loop through each time step to track the core position
    for t in range(vorticity_data.shape[1]):
        vorticity = vorticity_data[:, t]
        window_indices = np.where((x >= x_range[0]) & (x <= x_range[1]) & (y >= y_range[0]) & (y <= y_range[1]))[0]         # Apply the window range
        window_vorticity = vorticity[window_indices]
        core_idx = window_indices[np.argmax(window_vorticity)]          # the core index is the maximum vorticity in the window
        core_position = (x[core_idx], y[core_idx])                      # Obtaining the core position
        core_positions.append(core_position)                            # Append the core position to the list
    return core_positions

def wandering_corrected_flow_field(vorticity_data:float, core_positions:float, x:float, y:float):
    """ Function to perform wandering correction on the flow field using the core positions extracted from track_core_positions.
    Args:
        vorticity_data (float): The extracted vorticity data from the flow field. Can be any data but typically vorticity.
        core_positions (float): List of core positions at each time step.
        x (float): 1D array of the x coordinates of the mesh nodes. (in plane)
        y (float): 1D array of the y coordinates of the mesh nodes. (in plane)  
    Returns:
        corrected_vorticity (float): The corrected vorticity data after applying the wandering correction.
        mean_corrected_vorticity (float): The mean of the corrected vorticity data.
    """
    N = vorticity_data.shape[1]
    corrected_vorticity = np.zeros_like(vorticity_data)
    for t in range(N):
        core_pos = core_positions[t]
        # shifting the core position to the origin
        distances = np.sqrt((x - core_pos[0])**2 + (y - core_pos[1])**2)
        closest_idx = np.argmin(distances)
        # shifting the vorticity field by the the index of the core position
        corrected_vorticity[:, t] = np.roll(vorticity_data[:, t], -closest_idx)
    # The mean wandering corrected vorticity
    mean_corrected_vorticity = np.mean(corrected_vorticity, axis=1)
    return corrected_vorticity, mean_corrected_vorticity

def triple_decomposition(u:float, corrected_u:float, mean_corrected_u:float):
    """ function to perform Hussain's triple decomposition (Reynolds and Hussain 1972) on the flow field.
    Args:
        u (float): Baseline flow field data, typically the velocity field, in shpae (n, t) where n is the number of nodes and t is the number of time steps.
        corrected_u (float): The corrected flow field data after applying the wandering correction.
        mean_corrected_u (float): The mean of the corrected flow field data.
    Returns:
        u (float): The original flow field data.
        u_mean (float): The mean of the flow field data.
        u_coherent (float): The coherent velocity component of the flow field.
        u_fluc (float): The fluctuating velocity component of the flow
    """
    u_mean = np.mean(u, axis=1)
    # the fluctuating component
    u_fluc = corrected_u - mean_corrected_u[:, np.newaxis]
    # The coherent velocity component 
    u_coherent = (u - u_mean[:,np.newaxis]) - (corrected_u - mean_corrected_u[:,np.newaxis])
    return u, u_mean, u_coherent, u_fluc


def triple_decomposition_main(data_path, var, mesh_path, cut, vortex):
    """
    Main function that Applies Hussain's triple decomposition (Reynolds and Hussain 1972) to the extracted data.
    
    Args:
        data_path (str): Path to the HDF5 file containing the extracted data.
        var (str): The variable of interest to be processed.
        mesh_path (str): Path to the HDF5 file containing the mesh data.
        cut (str): The cut type to determine the window ranges (e.g., 'PIV1', 'Cut_030TE').
    
    Returns:
        u, u_coherent, u_random: Decomposed components of the flow field.
    """
    # Load the mesh data from the mesh file
    with h5py.File(mesh_path, 'r') as h5f:
        #x = h5f['0000/instants/0000/variables/x_node'][:]
        x = h5f['0000/instants/0000/variables/y_node'][:]
        y = h5f['0000/instants/0000/variables/z_node'][:]

    # Load the timeseries data from the data file
    with h5py.File(data_path, 'r') as f:
        data = f[var][:]
    print('Applying Triple Decomposition to the cutplane {0:s}, vortex type {1:s}'.format(cut, vortex))
    # Get the window ranges based on the cut
    window_ranges = get_window_ranges(cut, vortex)

    # Use one of the window ranges for the core position tracking (choose SV_Window for example)
    SV_WindowLL = window_ranges.get('SV_WindowLL')
    SV_WindowUR = window_ranges.get('SV_WindowUR')

    x_range = (np.min((window_ranges['LL'][0],window_ranges['UR'][0])), np.max((window_ranges['LL'][0],window_ranges['UR'][0])))
    y_range = (np.min((window_ranges['LL'][1],window_ranges['UR'][1])), np.max((window_ranges['LL'][1],window_ranges['UR'][1])))
    
    # Track core positions
    core_positions = track_core_positions(data, x, y, x_range, y_range)

    # Perform wandering correction
    corrected_vorticity, mean_corrected_vorticity = wandering_corrected_flow_field(data, core_positions, x, y)

    # Apply triple decomposition
    u, u_mean, u_coherent, u_random = triple_decomposition(data,corrected_vorticity, mean_corrected_vorticity)

    return u, u_mean, u_coherent, u_random