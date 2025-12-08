"""
modal_analysis package

This package contains modules for performing modal analysis on fluid dynamics data, including data extraction, surface extraction, and triple decomposition techniques.
The core functionalities inculde data handling, modal decomposition methods (DMD, SPOD, POD), and utilities for timing and logging.

"""

from .ModalAnalysis import ModalAnalysis, parse_arguments
from .map_cut import map_cut
from .triple_decomposition import get_window_ranges, track_core_positions
from .utils import print, timer, constants, init_logging_from_cut