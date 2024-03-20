#!/bin/bash

#SBATCH --account=indikar99
#SBATCH --partition=standard
#SBATCH --mail-user=cstansbu@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=180G
#SBATCH --time=36:00:00

SLURM_CONFIG='config/slurm_CPU'
CONFIG='config/config.yaml'
OUTPUT_PATH="$(cat ${CONFIG} | shyaml get-value output_path)"
THREADS=36

### export the environment 
conda env export > environment.yml

## build the workflow from the most current snakefile
cp Snakefile workflow.smk
echo "Built Workflow..."

# RUN
snakemake --profile ${SLURM_CONFIG} --cores ${THREADS} --latency-wait 90 -s workflow.smk --verbose --use-conda --rerun-incomplete
