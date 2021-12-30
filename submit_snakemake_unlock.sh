#!/usr/bin/bash
#SBATCH --time 10:00:00
#SBATCH --partition exacloud
#SBATCH --job-name workflow_submission
#SBATCH --output=logs/workflow_submission_%j.log

snakemake --unlock --rerun-incomplete
