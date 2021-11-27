#!/bin/bash --login
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=64G
#SBATCH --job-name=mapper-1000-genes

#install modules needed in the script

pip install numpy
pip install matplotlib.pyplot
pip install pandas
pip install sklearn.cluster
pip install sklearn.preprocessing
pip install sklearn
pip install IPython.display

# run script and save verbose into file

python3 /mnt/home/f0103237/hpcc_mapper_script.py > ~/mapper_python_output_file.txt

scontrol show job $SLURM_JOB_ID