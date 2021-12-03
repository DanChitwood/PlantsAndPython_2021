#!/bin/bash --login
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=80G
#SBATCH --job-name TSNE_all_genes

#install modules needed in the script

pip install openTSNE
pip install matplotlib.pyplot
pip install pandas
pip install seaborn

# run script and save verbose into file

python3 /mnt/home/f0103228/Project_docs/scripts/hpcc_script2.py > ~/Project_docs/results/Python_output_file.txt

scontrol show job $SLURM_JOB_ID
