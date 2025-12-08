#!/usr/bin/env python3
import sys
import os
from argparse import Namespace

# ---------------------------------------------------------------------
# Ensure package import works correctly
# ---------------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

# Import the package CLI entrypoint
from modal_analysis_core.__main__ import main


# ---------------------------------------------------------------------
# Manual configuration to override CLI arguments
# ---------------------------------------------------------------------
config = {
    "compute_mode"        : "DMD",
    "data_type"           : "EXTRACT",
    "var"                 : "pressure",
    "sol_dir"             : "/project/rrg-plavoie/denggua1/BBDB_10AOA/RUN_ZONE_Apr24/Isosurface/Cut_PIV2_VGT/",
    "sol_file"            : "B_10AOA",
    "mesh_dir"            : "/project/rrg-plavoie/denggua1/BBDB_10AOA/MESH_ZONE_Apr24",
    "mesh_file"           : "Bombardier_10AOA_Combine_Apr24.mesh.h5",
    "option"              : 3,
    "nskip"               : 1,
    "nstart"              : 100,
    "maxfile"             : 400,
    "extract_only"        : False,
    "reload_source"       : False, 
    "skip_extract"        : True,
    "noise_lvl"           : 1e-4,
    "dt"                  : 1.67e-8*2000,
    "df"                  : 300,
    "freq"                : [3000],
    "nperiod"             : 3,
    "triple_decomp"       : False,
    "vortex"              : "PV",
}
config["dest_dir"]  = os.getcwd()
config["dest_file"] = os.path.basename(os.path.normpath(config["sol_dir"]))
config["dmd_output_dir"] = 'DMD_Output_'+ config["var"] + '_'+ config["dest_file"]
args = Namespace(**config)

# ---------------------------------------------------------------------
# Run using exact CLI entrypoint logic
# ---------------------------------------------------------------------
main(args)