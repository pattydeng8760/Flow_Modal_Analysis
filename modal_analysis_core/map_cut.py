import re
import numpy as np
from .utils import constants

# Mapping the cut selection location
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
    print('     The selected cut is of style: {0}'.format(data_type))
    print('     The selected cut origin is: {0}'.format(origin))
    print('     The selected cut normal is: {0}'.format(normal))
    return origin,normal