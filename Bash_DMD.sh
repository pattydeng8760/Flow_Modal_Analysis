#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-01:00:00
#SBATCH --job-name=DMD_A10
#SBATCH --mail-user=patrickgc.deng@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-moreaust-ac

source /project/m/moreaust/Env/avbpNpy_env.sh
use_py_tools

# file="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/Cut_midspan_TE_VGT/"
# python DMD.py $file

# file="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/Cut_25mm_tip_VGT/"
# python DMD.py $file

file="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/Cut_5mm_tip_VGT/"
python DMD.py $file

# file="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/Cut_015_TE_VGT/"
# python DMD.py $file

# file="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/Cut_030_TE_VGT/"
# python DMD.py $file

# file="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/Cut_085_TE_VGT/"
# python DMD.py $file

# file="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/Cut_095_TE_VGT/"
# python DMD.py $file

# file="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/Cut_PIV1_VGT/"
# python DMD.py $file

# file="/home/p/plavoie/denggua1/scratch/Bombardier_LES/B_10AOA_LES/Isosurface/Cut_PIV2_VGT/"
# python DMD.py $file



