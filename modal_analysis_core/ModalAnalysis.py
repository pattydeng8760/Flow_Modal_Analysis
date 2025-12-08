import os, glob, shutil
import numpy as np
import argparse
from datetime import datetime
from types import SimpleNamespace
from .extract import extract_files, extract_surface, extract_data
from .dmd import dmd_compute, reconstruct_dmd_field

def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute modal decomposition from transient flow field data.")
    # The input arguments
    parser.add_argyment("--compute-mode", "-cm", type=str, default="DMD", choices=['DMD', 'SPOD', "POD"], help="Modal decomposition mode to perform.")
    parser.add_argument("--data-type", "-d", type=str, default="EXTRACT",choices=['FWH', 'OUTPUT', 'EXTRACT', 'CLIP'], help="Type of source data file.")
    parser.add_argument("--var", "-v", type=str, default="pressure", help="variable of interest. Must be consistent with source data naming.")
    parser.add_argument("--sol-dir", "-s", type=str, required=True, help="Directory of the source transient solution.")
    parser.add_argument("--sol-file", "-f", type=str, required=True, help="Preamble string of the source transient solution files.")
    parser.add_argument("--dest-dir", "-o", type=str, default=os.getcwd(), help="Destination directory for output data.")
    parser.add_argument("--dest-file", "-D", type=str, default=None, help="Destination subdirectory name for extracted and filtered raw data. If not provided, uses the basename of sol-dir. \
                        If provided, this will copy the files over from the sol-dir to the provided loation.")
    # Mesh information
    parser.add_argument("--mesh-dir", "-m", type=str, required=True, help="Directory of the mesh file.")
    parser.add_argument("--mesh-file", "-M", type=str, required=True, help="Mesh file name in .h5 format.")
    parser.add_argument("--surface-patches", "-sp",type=str, nargs="+",default=["Airfoil_Surface"], help="List of surface patches to include in the analysis.")
    # Output information
    parser.add_argument("--dmd-output-dir", "-O", type=str, default=None, help="Destination directory for DMD output data. If not provided, uses 'DMD_Output_<var>_<dest-file>'.")
    parser.add_argument("--extract-only", action="store_true", help="Only extract the files and data, do not perform modal analysis.")
    # File sampling option
    parser.add_argument("--option", "-p", type=int, default=3, choices=[1, 2, 3], help="File sampling option: 1 = take n files (maxfile), 2 = skip every n files until maxfiles or end, 3 = take all files.")
    parser.add_argument("--nskip", "-n", type=int, default=1, help="Number of files to skip (used if option 2 is selected).")
    parser.add_argument("--nstart", type=int, default=0, help="Starting file index.")
    parser.add_argument("--maxfile", "-F", type=int, default=None, help="Maximum number of files to process.")
    parser.add_argument("--reload-source", type=bool, default=False, help="Re-copy files from source directory.")
    parser.add_argument("--skip-extract", type=bool, default=False, help="Skip copying process and re-load raw data from copied source file directory.")
    parser.add_argument("--noise-lvl", type=float, default=1e-4, help="Noise level for processing.")
    parser.add_argument("--triple-decomp", action="store_true", help="Perform triple decomposition of the flow field.")
    parser.add_argument("--vortex", type=str, default="PV", help="Type of vortex for triple decomposition.")
    # Data time information
    parser.add_argument("--dt", type=float, default=1.0, help="Time step between consecutive files.")
    parser.add_argument("--df", type=int, default=1, help="Frequency sampling internal (frequency bin width) for each sampling frequency.")
    parser.add_argument("--freq", type=float, nargs="+", default=[1000], help="Center frequency for Modal Analysis (only applies to DMD and SPOD).")
    parser.add_argument("--nperiod", type=int, default=3, help="Number of periods to reconstruct for DMD field reconstruction.")
    
    return parser.parse_args(argv)

class ModalAnalysis():
    """
    Class to perform modal analysis on transient flow field data.
    """
    def __init__(self, args):
        text = "Beginning Modal Analysis"
        print(f'\n{text:=^100}\n')
        
        # Store input arguments
        self.args = args
        self.compute_mode = self.args.compute_mode
        self.var = self.args.var
        self.data_type = self.args.data_type
        # Mesh information
        self.mesh = SimpleNamespace()
        self.mesh.mesh_dir = self.args.mesh_dir
        self.mesh.mesh_file = self.args.mesh_file
        if self.data_type == 'FWH':
            self.mesh.surface_patches = self.args.surface_patches
        # File sampling options
        self.sample = SimpleNamespace()
        self.sample.option = self.args.option
        self.sample.nskip = self.args.nskip
        self.sample.nstart = self.args.nstart
        self.sample.maxfile = self.args.maxfile
        
        # Data time information
        self.data = SimpleNamespace()
        self.data.dt = self.args.dt
        self.data.df = self.args.df
        self.data.freq = np.array(self.args.freq)
        self.data.nperiod = self.args.nperiod
        
        # Triple decomposition options
        self.triple_decomp = self.args.triple_decomp
        self.vortex = self.args.vortex
        self.fluc = True if self.triple_decomp else False
        
        # Making the output directory for the DMD data
        self.modal_output = os.path.join(self.args.dest_dir, self.args.dmd_output_dir)
        if os.path.exists(self.modal_output) and self.args.reload_source==True:
            shutil.rmtree(self.modal_output)
        os.makedirs(self.modal_output, exist_ok = True)
    
        self._print_args()
        
    def _print_args(self):
        """Print all CLI / input arguments stored in self.args."""
        text = " Modal Analysis Input Arguments "
        print(f'\n{text:.^80}\n')
        for name, value in vars(self.args).items():
            print(f"   {name:20s}: {value}")
        print("   compute date        : ", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        text = " End of Input Arguments "
        print(f'\n{text:.^80}\n')
    
    def pre_process(self):
        """Pre-process the input data for extraction and curation before modal analysis."""
        # Extracting the files
        self.output = extract_files(self.args.sol_dir,self.args.sol_file,self.args.dest_dir,self.args.dest_file,
                                self.args.dmd_output_dir,option=self.sample.option,maxfile=self.sample.maxfile,nskip=self.sample.nskip,reload=self.args.reload_source,skip=self.args.skip_extract)
        # Extracting the mesh surface
        self.mesh.DMD_mesh, self.mesh.mesh_node = extract_surface(self.mesh.mesh_dir, self.mesh.mesh_file, self.modal_output, self.data_type, self.args.dest_file)
        # Extracting the data
        self.data.ntime, self.data.DMD_source = extract_data(self.args.sol_dir, self.modal_output, self.mesh.mesh_node, self.args.dest_file, self.data_type, self.var, 
                                    option=self.sample.option, nstart=self.sample.nstart,maxfile=self.sample.maxfile,nskip=self.sample.nskip, fluc=self.fluc, load_existing=self.args.reload_source)
        if self.args.extract_only:
            text = "Extraction Only Mode Enabled - Skipping Modal Analysis"
            print(f'\n{text:.^80}\n')
            text = "Modal Analysis Complete"
            print(f'\n{text:=^100}\n')
    def run_modal_analysis(self):
        if self.compute_mode == "DMD":
            # Performing DMD computation
            self.data.peak_modes, self.data.DMD_result, self.data.DMD_spectra, self.data.DMD_modes = dmd_compute(self.mesh.DMD_mesh, self.data.DMD_source, self.modal_output, 
                                            self.data.ntime, self.var, self.data.dt, self.data.df, self.data.freq, self.data_type, noise_lvl=self.args.noise_lvl, 
                                            triple_decomp=self.triple_decomp, vortex=self.vortex)
            # Reconstructing the DMD field
            reconstruct_dmd_field(self.data.DMD_modes, self.data.DMD_spectra, self.data.peak_modes, self.modal_output, self.var, self.args.dest_file, nperiod=self.data.nperiod)
        text = "Modal Analysis Complete"
        print(f'\n{text:=^100}\n')