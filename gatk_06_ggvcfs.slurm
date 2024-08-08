#!/bin/bash
#SBATCH --job-name=ggvcfs
#SBATCH -A beagle
#SBATCH -t 14-00:00:00
#SBATCH -N 1-1
#SBATCH -n 128
#SBATCH --mem=0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yin168@purdue.edu
#SBATCH -o /scratch/bell/yin168/pflavescens/dbslurm/ggvcfs.o
#SBATCH -e /scratch/bell/yin168/pflavescens/dbslurm/ggvcfs.e

cd /scratch/bell/yin168/pflavescens/dbgeno

/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" GenotypeGVCFs -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -V gendb:///scratch/bell/yin168/pflavescens/dbgeno/gdb -O yp_ggvcfs.vcf.gz --tmp-dir /scratch/bell/yin168/pflavescens/dbgeno/tmpggvcfs --max-alternate-alleles 100
