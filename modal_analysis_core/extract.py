import os
import shutil
import glob
import numpy as np
import h5py
from antares import Reader, Writer, Treatment, Family
from .utils import constants
from .map_cut import map_cut 


# Extracting files in the source directory into the present DMD directory to be conditioned.
def extract_files(sol_dirName:str, sol_file:str, dest_dir:str, dest_File:str, DMD_output_dir:str\
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
    if skip:
        print('     The extraction is skipped...')
        print('     Using the existing files in the destination directory: {0:s}{1:s}_*'.format(sol_dirName,sol_file))
        output = os.path.join(sol_dirName,sol_file)
        arr = sorted(glob.glob('{0:s}/{1:s}*.h5'.format(sol_dirName,sol_file)))
        n_files = len(arr)
        print('     The total number of timesteps in the source directory is: {0:s}'.format(str(n_files)))
    else:
        # An array with all the files in the source directory
        arr,arr_dest = [], []
        arr += sorted(glob.glob('{0:s}/{1:s}*.h5'.format(sol_dirName,sol_file)))
        n_files = len(arr)
        if os.path.exists(output):
            print('----> Checking the existing files in the destination directory...')
            print('     There are {0:d} files in the source directory'.format(len(arr)))
            arr_dest += sorted(glob.glob('{0:s}/*.h5'.format(output)))
            print('     There are {0:d} files in the destination directory'.format(len(arr_dest)))
            print('     After extraction there should be {0:d} files in the destination directory'.format(int(np.floor(len(arr))/nskip)))
        # Checking if all the files from the source directory are already extracted
        if ((int(len(arr_dest)) == int(np.floor(len(arr))/nskip) and os.path.exists(output)) or (int(len(arr_dest)) >= maxfile)) or (reload == False and (int(len(arr_dest)) >= maxfile)):
            text = '     All the files are already extracted to destination directory'
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
                for i in range(0,np.min([int(len(arr)/nskip),int(maxfile)])):
                    source = os.path.join(sol_dirName,arr[int(i*nskip)])
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
    return output


# Extracting the airfoil surface mesh for the base cut_locationof DMD
def extract_surface(mesh_fileDir:str,mesh_fileName:str,DMD_output:str, data_type:str, cut_location:str,\
    span:float=-0.2286, tip_gap:float=-0.1034, AoA:int=10, input_surface:list=["Airfoil_Surface"], reload:bool=False):
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
        input_surface (list): The list of surface patches name to be extracted from the mesh file - only for FWH data_type
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
        print('----> LES Mesh already extracted at: \n     {0:s}'.format(DMD_mesh))
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
        for surf_name in input_surface:
            try:
                airfoil_base[surf_name] = base.families['Patches'][surf_name]
                print('      Surface %s extracted' %surf_name)
            except: 
                print('      Warning: Surface %s not found in the mesh file.' %surf_name)
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
        print('\nThe Extracted surface mesh is saved in \n     {0:s}'.format(DMD_mesh))
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
        print('\n     The number of nodes on the surface is: {0:d}'.format(nodes))
        DMD_mesh = os.path.join(DMD_output,mesh_name+'.h5')
        print('\n     The Extracted surface mesh is saved in: \n     {0:s}'.format(DMD_mesh))
    text = 'Mesh Extraction Complete!'
    print(f'\n{text:.^80}\n')  
    return DMD_mesh, nodes


# Extracting data in the post sampled source directory into a unified file to be read by the DMD script
def extract_data(sol_dir:str,DMD_output:str,nodes:float,cut_location:str,data_type:str,var:str, \
    option:int=1,nstart:int=0,nskip:int=1,maxfile:int=1000,fluc:bool=False, load_existing:bool=False):
    """ Extracting the source data into a unified hdf5 file to be read by the DMD script and numpy file for backup
    Args:
        sol_dir (str): The source directory where the transient files are located 
        DMD_output (str): The destination directory where the extracted data is to be stored
        nodes (int): The number of nodes on the surface from the meah file (to compare with the data to ensure consistency)
        data_type (str): The type of data to be extracted, either 'FWH', 'OUTPUT', 'EXTRACT', or 'CLIP'
        cut_location (str): The location of the cut plane, defined either explicitly via PIV planes
        option (int): The option to extract the files, 1 for n files, 2 for skip every n files, 3 for all files
        nstart (int): The starting file index
        nskip (int): The number of files to skip
        maxfile (int): The maximum number of files to extract
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
    if not load_existing and os.path.exists(os.path.join(DMD_output,'DMD_Data_{0:s}_{1:s}_max_inst.hdf5'.format(var,cut_location))):
        print('----> Loading the data from the saved .hdf5 file with the full data')
        with h5py.File(os.path.join(DMD_output,'DMD_Data_{0:s}_{1:s}_max_inst.hdf5'.format(var,cut_location)),'r') as f:
            data = f[var][()]
            if option == 1:
                data = data[:,nstart:nstart+maxfile]
            elif option == 2:
                data = data[:,nstart::nskip]
                data = data[:,:maxfile]
            ntime = data.shape[1]
        DMD_source = os.path.join(DMD_output,'DMD_Data_{0:s}_{1:s}_{2:04d}inst.hdf5'.format(var,cut_location,ntime))
        print('     The number of timesteps extracted is: {0:d}'.format(ntime))
        print('     The maximum and mimumum {0:s} in the instantaneous field is: {1:.4f}, {2:.4f}'.format(var,np.max(data),np.min(data)))
        f=h5py.File(DMD_source,'w')
        data_float32 = data.astype(np.float64)
        f[var]=data_float32
        f.close()
    else:
        print('----> Loading the raw data from the copied directory')
        # The directory information
        l_orig=sorted(glob.glob('{0:s}/*.h5'.format(sol_dir)))
        l_count = np.linspace(0,len(l_orig)-1,len(l_orig),dtype=int)
        nb_source = len(l_orig)
        maxfile = np.min([maxfile,nb_source])
        if option == 1:
            print('     Option 1: Extracting the first {0:d} files sequentially:'.format(maxfile))
            l,l_count=l_orig[nstart:nstart+maxfile],l_count[nstart:nstart+maxfile]
        elif option == 2:
            print('     Option 2: Extracting every {0:d} files sequentially until {1:d}:'.format(nskip, maxfile))
            l_orig,l_count=l_orig[nstart::nskip],l_count[nstart::nskip]
            l,l_count=l_orig[:maxfile],l_count[:maxfile]
        elif option == 3:
            print('     Option 3: Extracting all files:')
            l_orig,l_count=l_orig,l_count
            l,l_count=l_orig,l_count
        nb_files=len(l)
        print('\n     The number of extracted files is: {0:d}.'.format(nb_files))
        print('     The variable of interest is: {0:s}.'.format(var))
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
                    print('     The number of nodes in the datafile is: {0:d} equal to that of the mesh file: {1:d}\n'.format(press.shape[0],nb_points))
                    print('----> Beginning data extraction from the source files:\n')
            if np.mod(it,10) == 0 or it==nb_files-1:
                print('      Iteration {0:03d} of {1:03d}: Extracting file {2:03d} of {3:03d} from the original'.format(it, nb_files, l_count[it], nb_source))
        data_mean = np.mean(data,axis=1)
        data_mean = data_mean.reshape((nb_points,1))
        print("----> Data extraction complete!")
        if fluc :
            print('----> Extracting the fluctuating part of the data')
            data = data - data_mean  # Subtracting the mean value to get the fluctuating part
        print('\n      The maximum and mimumum {0:s} in the instantaneous field is: {1:.4f}, {2:.4f}'.format(var,np.max(data),np.min(data)))
        # The number of timesteps
        ntime = data.shape[1]
        print('     The number of timesteps extracted is: {0:d}'.format(ntime))
        # Saving the output as hdf5 file 
        if option != 3: 
            DMD_source = os.path.join(DMD_output,'DMD_Data_{0:s}_{1:s}_{2:04d}inst.hdf5'.format(var,cut_location,ntime))
        else: 
            DMD_source = os.path.join(DMD_output,'DMD_Data_{0:s}_{1:s}_max_inst.hdf5'.format(var,cut_location))
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
    print('     The destination DMD data file is: \n     {0:s}.'.format(DMD_source))
    text = 'Data Extraction Complete!'
    print(f'\n{text:.^80}\n')  
    return ntime, DMD_source
