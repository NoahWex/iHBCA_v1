#!/bin/bash
#SBATCH --job-name=check_sync_help
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=00:05:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

module load python/3.10.2
hca-smart-sync sync --help
