#!/bin/sh
#SBATCH --account=astro
#SBATCH --job-name=band_opt_mu_1to2000_width1_snr0.25_bolo10
#SBATCH -c 3
#SBATCH --mem-per-cpu 5gb
#SBATCH --time=00-12:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as6131@columbia.edu

module load anaconda
source activate newinstrument

path_to_dir=/burg/home/as6131/specter_optimization
python $path_to_dir/code/run_band_optimization.py mu_1to2000_width1_snr0.25_bolo10.ini
