import os
import sys
import time
import builtins
from datetime import datetime

def setup_logging(log_file):
    sys.stdout = open(log_file, "w", buffering=1)

def init_logging_from_cut(mode,var, modal_output_dir):
    date = datetime.now().strftime("%Y-%m-%d")  # e.g. "2025-12-08"
    plane = modal_output_dir.replace('_Output', '')
    log_file = f"log_{mode}_analysis_{var}_{plane}_{date}.txt"
    setup_logging(log_file)

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