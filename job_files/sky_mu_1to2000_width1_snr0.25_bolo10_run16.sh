#!/bin/sh
#SBATCH --account=astro
#SBATCH --job-name=sky_mu_1to2000_width1_snr0.25_bolo10_run16
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00-1:30
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as6131@columbia.edu

module load anaconda/3-2020.11
conda init bash
conda activate newinstrument

path_to_dir=/burg/home/as6131/specter_optimization
python $path_to_dir/code/run_skymodel_bias.py mu_1to2000_width1_snr0.25_bolo10.ini 16
