
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