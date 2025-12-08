import os
import numpy as np
import sys
import time
from antares import Reader, Writer, Base, Zone, Instant, Treatment
import h5py
import copy
import glob
import builtins
import re
import shutil
import matplotlib.pyplot as plt

# the Main DMD function compute
def dmd_compute(DMD_mesh:str, DMD_source:str, DMD_output:str,\
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
    t0 = time.time()  
    # Initializing the data matrix
    data = {}
    # Number of instants to account for in DMD, Must match with data file
    ninst = ntime
    # Add frequencies of interest 
    freqs = freq_interest
    print('     The frequency range of interest is: {0:s} Hz.'.format(str(freq_interest)))
    print('     The DMD sampling frequency is: {0:1.1f} Hz.'.format(1/dt))
    print('     The DMD sweep range is +/-: {0:1.0f} Hz at each frequency range.'.format(df))
    # Frequency interval: the algorithm will search for modes
    # with maximum amplitude for selected frequencies within given frequency interval:
    #################################################################################### 
    if data_type == 'FWH':
        #Open mesh
        print('\n----> Open Mesh file')
        print('     The DMD Mesh file is located at \n     {0:s}'.format(DMD_mesh))
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
        print('     The DMD Mesh file is located at \n     {0:s} '.format(DMD_mesh))
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
        print('     The DMD Mesh file is located at \n     {0:s} '.format(DMD_mesh))
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
        print('     The data source file is located at {0:s} '.format(DMD_source))
        fin = h5py.File(DMD_source,'r')
        data[var] = fin['/{0:s}'.format(var)][:,0:ninst]
        fin.close()
        del fin
    # Create a new base to make DMD
    b = Base()
    b['0'] = Zone()
    b[0].shared['cell_volume'] = mesh[0][0]['cell_volume']
    print("----> Checking for NaN or Inf values in the data...")
    print('\n     Instant   isinf    isnan')
    for i in range(int(data[var].shape[1])):
        instant_name = 'snapshot_%s' % i
        b[0][instant_name] = Instant()
        b[0][instant_name][var] = data[var][:,i].flatten(order='F')
        if np.mod(i,10)==0:
            text=('     {0:s}  {1:}  {2:} '.format(instant_name, np.isinf(np.array(b[0][instant_name][var])).any(), np.isnan(np.array(b[0][instant_name][var])).any()))
            print(f'{text}') 
        if np.isinf(np.array(b[0][instant_name][var])).any() or np.isnan(np.array(b[0][instant_name][var])).any():
            raise ValueError('The data contains NaN or Inf values')
    #################################################################################### 
    # Perform DMD
    print('\n----> Performing DMD Computation...')
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
    print(' DMD Computation Complete')

    ####################################################################################
    # Extracting the DMD Data
    print('\n----> Extracting DMD Data')
    print('     The peak modes are:')
    print(res_dmd['modes'].attrs['peak_modes_{0:s}'.format(var)])
    peak_modes = res_dmd['modes'].attrs['peak_modes_{0:s}'.format(var)]
    w = Writer('hdf_antares')
    DMD_result = os.path.join(DMD_output,'DMD_{0:s}_{1:04d}inst'.format(var,ninst))
    print('     The DMD result file is saved as {0:s} '.format(DMD_result))
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
    print('     The DMD spectra file is saved as {0:s} '.format(DMD_spectra))
    np.save(DMD_spectra,spectrum)
    # Write modes to mesh file for visual post-processing
    modes = len(res_dmd['modes'].keys())
    for mode in peak_modes:
        mesh[0][0]['{0:s}_{1:s}_mod'.format(var,mode)] = res_dmd['modes'][mode]['{0:s}_mod'.format(var)]
        mesh[0][0]['{0:s}_{1:s}_phi'.format(var,mode)] = res_dmd['modes'][mode]['{0:s}_phi'.format(var)]
    del mesh[0][0]['cell_volume']
    w = Writer('hdf_antares')
    DMD_modes = os.path.join(DMD_output,'DMD_{0:04d}inst_{1:s}_modes'.format(ninst,var))
    print('     The DMD modes file is saved as {0:s} '.format(DMD_modes))
    w['filename'] = DMD_modes
    w['base'] = mesh
    w['dtype'] = 'float32'
    w.dump()
    t1 = time.time()
    print("\n The total DMD compute time is: {0:1.0f} min, {1:1.0f} s".format(np.floor((t1-t0)/60), np.mod((t1-t0),60)))
    text = 'DMD Complete!'
    print(f'\n{text:.^80}\n')
    return peak_modes, DMD_result, DMD_spectra, DMD_modes 


def plot_dmd_spectra(f, amplitude, amplification, DMD_output, cut_location, var):
    
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


def reconstruct_dmd_field(DMD_result, DMD_spectra, peak_modes, DMD_output, var, cut_location, dtype:str='float32', nperiod:int=2):
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
    print('     The loaded DMD result file: {0:s}'.format(file))
    print(f"     The contents of this file: {b[0][0].keys()}")
    print('\n----> Loading the DMD Spectra File')
    print('     Saving the DMD spectra')
    data = np.load(DMD_spectra)
    f = data[:,0]
    amplitude = data[:,1]
    amplification = data[:,2]
    # Plotting the DMD spectra
    plot_dmd_spectra(f, amplitude, amplification, DMD_output, cut_location, var)
    
    print('\n----> Reconstructing the {0:s} Data...'.format(var))
    for i in range(0,np.shape(peak_modes)[0]):
        mode_int = peak_modes[i]
        mode_number = int(mode_int.replace('H',''))
        freq = int(f[mode_number])

        print('     Reconstructing Mode {0:s} corresponding to {1:1.0f} Hz'.format(peak_modes[i],freq))
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
        print('     The DMD Result for Mode {0:s} is extracted.'.format(peak_modes[i]))
        print('     The DMD Result saved at: {0:s}'.format(output_path))
    text = 'Reconstructing Field Complete!'
    print(f'\n{text:.^80}\n')  