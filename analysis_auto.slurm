#!/bin/bash
#SBATCH -p mid1
#SBATCH -A umgerlevik
#SBATCH -J analysis
#SBATCH -N 1
#SBATCH --ntasks=20
#SBATCH --time=04-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=umut.gerlevik@gmail.com
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --no-requeue

export PATH=/truba/home/umgerlevik/softwares/vmd/bin:$PATH

./analysis_auto.sh


