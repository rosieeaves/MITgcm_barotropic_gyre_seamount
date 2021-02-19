#!/bin/bash
#SBATCH --job-name=baro_wind_mount
#SBATCH --output=baro_wind_mount-%j.out
#SBATCH --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=rosie.eaves@univ.ox.ac.uk
#SBATCH -p priority-ocean
# maybe even 400 or 150 for mempercpu

# Load required modules
module purge
module load otheros/generic-arc
module load intel-compilers/2013
module load intel-mpi/2013
export CPP=''

srun ./mitgcmuv
