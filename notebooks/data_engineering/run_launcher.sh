#!/bin/bash
#SBATCH --job-name=merge_pod5
#SBATCH --mail-user=cstansbu@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=96:00:00
#SBATCH --account=indikar1
#SBATCH --partition=standard
#SBATCH --output=/home/cstansbu/git_repositories/scpc/notebooks/data_engineering/launcher.log

/home/cstansbu/git_repositories/scpc/notebooks/data_engineering/merge_pod5_launch.sh

paramrun

