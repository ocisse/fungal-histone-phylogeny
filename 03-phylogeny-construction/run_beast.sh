#!/bin/bash 

#SBATCH --job-name=beast
##SBATCH --partition=unlimited
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8g
#SBATCH --time=07-00:00:00
#SBATCH --output=beast_wopart_output.log
#SBATCH --error=beast_wopart_error.log

module load beagle-lib/3.1.2

beast -threads $SLURM_CPUS_PER_TASK histones_wpartition.xml
beast -threads $SLURM_CPUS_PER_TASK histones_wopartition.xml 

