# README

================================================
- Jeremy Arsenault

- Francisco Guerra

- Nick Johnson

- Quintero Corrales Christian

- Anne Steensma

- Damián Villaseñor Amador

=================================================


This repository contains all scripts used to generate a mapper graph usign RNAseq data from Arabidopsis.

`helper_functions.py` File with the main mapper analyses functions used for the Arabidpsis RNAseq data

`hpcc_mapper_script.py` _Arabidopsis_ gene sequence mapper graph python script using centroid lens

`hpcc_mapper_script_notebook.ipynb` _Arabidopsis_ gene sequence expression mapper graph jupyter notebook file

`slurm_mapper_script.sh` script for slurm running


==================================================


To run the mapper graph, please do it through the HPCC, because in the database isn't available in this repo.

Instructions to run the mapper graph in the HPCC:

`$ ssh dev-intel14-k20` # choose a developmental node (whichever is available)

`$ sbatch /mnt/ufs18/rs-008/HRT841_F21/Class_project/mapper/slurm_mapper_script.sh`

`$ squeue -u $USER` # check if the job is running

When the job ends, check the directory for an .html file: it will be the kmapper graph, download it and open it in a web browser.
