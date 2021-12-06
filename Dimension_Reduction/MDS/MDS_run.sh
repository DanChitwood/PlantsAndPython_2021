#!/bin/bash --login
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=15G
#SBATCH --job-name MDS_100_genes

module purge
module load iccifort/2020.4.304
module load impi/2019.9.304
module load scikit-learn/0.23.2
pip install matplotlib
pip install seaborn

python MDS_100_genes.py
# can change this depending on if you're running the 100 gene subset or the full dataset
# note: the 100 gene subset took 1 hr and 5 Gb to run

scontrol show job $SLURM_JOB_ID
