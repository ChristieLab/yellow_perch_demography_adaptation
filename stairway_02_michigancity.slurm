#!/bin/bash
#SBATCH --job-name=03_stairway_michigancity
#SBATCH -A beagle
#SBATCH -t 14-00:00:00
#SBATCH -N 1-1
#SBATCH -n 20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yin168@purdue.edu
#SBATCH -o /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/03_stairway_michigancity.o
#SBATCH -e /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/03_stairway_michigancity.e

cd /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01

bash stiarway_02_fold_michigancity_missmaxno_gen02.blueprint.sh
