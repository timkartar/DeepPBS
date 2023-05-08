#!/bin/bash
#SBATCH --time=1440
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32GB
#SBATCH --mail-user=raktimmi@usc.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -o ./report/%x.%j.out 
#SBATCH --partition=rohs
#SBATCH --gres=gpu:1


python -W ignore driver.py  ./folds/train$1.txt ./folds/valid$1.txt  -c config.json --balance unmasked --eval_every 1
