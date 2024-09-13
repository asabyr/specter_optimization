#!/bin/sh
#SBATCH --account=astro
#SBATCH --job-name=grid_mu_1to2000_width1_snr0.25_bolo10_run13
#SBATCH -N 1
#SBATCH --ntasks-per-node=30
#SBATCH --array=0-29
#SBATCH --time=00-12:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=as6131@columbia.edu

module load anaconda/3-2020.11
conda init bash
conda activate newinstrument

path_to_dir=/burg/home/as6131/specter_optimization
python $path_to_dir/code/run_detector_optimization.py mu_1to2000_width1_snr0.25_bolo10.ini 13 ${SLURM_ARRAY_TASK_ID}
