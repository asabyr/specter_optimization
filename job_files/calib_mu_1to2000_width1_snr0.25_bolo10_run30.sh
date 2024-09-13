#!/bin/sh
#SBATCH --account=astro
#SBATCH --job-name=calib_mu_1to2000_width1_snr0.25_bolo10_run30
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as6131@columbia.edu

module load anaconda/3-2020.11
source activate newinstrument

path_to_dir=/burg/home/as6131/specter_optimization
python $path_to_dir/code/run_calibration_optimization.py mu_1to2000_width1_snr0.25_bolo10.ini 30
