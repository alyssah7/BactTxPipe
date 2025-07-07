#!/bin/bash
#SBATCH --job-name=ved_snake
#SBATCH -p sched_mit_hill,sched_any,newnodes,sched_mit_tami
#SBATCH -n 1
#SBATCH --time=13:00:00
#SBATCH -o mainout.txt
#SBATCH -e mainerr.txt
#SBATCH --mem=2000
#SBATCH --mail-user=alyssah@mit.edu
#SBATCH --mail-type=FAIL

# Activate conda environment (may need to change name of env)
source activate snakemake

bash snakemakeslurm.sh

echo Done!!!
