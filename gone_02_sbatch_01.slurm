#!/bin/bash
#SBATCH --job-name=13_sbatch_01
#SBATCH -A beagle
#SBATCH -t 14-00:00:00
#SBATCH -N 1-1
#SBATCH -n 5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yin168@purdue.edu
#SBATCH -o /scratch/bell/yin168/pflavescens/gt_w210_01/fig3_01/fig3c/slurm/13_sbatch_01.o
#SBATCH -e /scratch/bell/yin168/pflavescens/gt_w210_01/fig3_01/fig3c/slurm/13_sbatch_01.e

cd /scratch/bell/yin168/pflavescens/gt_w210_01/fig3_01/fig3c

for i in select_ypsnp_greenbay_mindp03_ind090_bi_maf005 select_ypsnp_muskegon_mindp03_ind090_bi_maf005
do
for j in {1..50}
do
sbatch /scratch/bell/yin168/pflavescens/gt_w210_01/fig3_01/fig3c/${i}_$"rep"$(printf '%04d' $j)_gone/GONE/Linux/${i}_$"rep"$(printf '%04d' $j)_GONE.sh
done
done
