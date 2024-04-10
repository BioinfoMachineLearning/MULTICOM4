#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition=gpu
#SBATCH --account=chengji-lab
#SBATCH --ntasks-per-node=4  # cores per task
#SBATCH --mem=60G  # memory per core (default is 1GB/core)
#SBATCH --time 0-06:00     # days-hours:minutes time
#SBATCH --gres gpu:1
#SBATCH --job-name=JOBNAME
#SBATCH --output=JOBNAME-%j.out  # %j is the unique jobID

source /cluster/pixstor/chengji-lab/bml_casp16/anaconda3/bin/activate
conda activate multicom4
export PYTHONPATH=/cluster/pixstor/chengji-lab/bml_casp16/MULTICOM4
cd /cluster/pixstor/chengji-lab/bml_casp16/MULTICOM4

