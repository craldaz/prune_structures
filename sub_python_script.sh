#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --gres gpu:1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name autoNEB

python generate.py confgen_small.inp
