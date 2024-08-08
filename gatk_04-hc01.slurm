#!/bin/bash
#SBATCH --job-name=hc01
#SBATCH -A beagle
#SBATCH -t 14-00:00:00
#SBATCH -N 1-1
#SBATCH -n 1
#SBATCH --mem=0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yin168@purdue.edu
#SBATCH -o /scratch/bell/yin168/pflavescens/dbslurm/hc01.o
#SBATCH -e /scratch/bell/yin168/pflavescens/dbslurm/hc01.e

cd /scratch/bell/yin168/pflavescens/dbgeno

/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -I D9_320_FDSW202502695-1r_dedup.bam -O D9_320_FDSW202502695-1r.g.vcf.gz -ERC GVCF --max-alternate-alleles 100 --max-num-haplotypes-in-population 256
/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -I D9_321_FDSW202502696-1r_dedup.bam -O D9_321_FDSW202502696-1r.g.vcf.gz -ERC GVCF --max-alternate-alleles 100 --max-num-haplotypes-in-population 256
/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -I D9_322_FDSW202502697-1r_dedup.bam -O D9_322_FDSW202502697-1r.g.vcf.gz -ERC GVCF --max-alternate-alleles 100 --max-num-haplotypes-in-population 256
/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -I D9_325_FDSW202502700-1r_dedup.bam -O D9_325_FDSW202502700-1r.g.vcf.gz -ERC GVCF --max-alternate-alleles 100 --max-num-haplotypes-in-population 256
/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -I D9_326_FDSW202502701-1r_dedup.bam -O D9_326_FDSW202502701-1r.g.vcf.gz -ERC GVCF --max-alternate-alleles 100 --max-num-haplotypes-in-population 256
/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -I D9_327_FDSW202502702-1b_dedup.bam -O D9_327_FDSW202502702-1b.g.vcf.gz -ERC GVCF --max-alternate-alleles 100 --max-num-haplotypes-in-population 256
/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -I D9_328_FDSW202502703-1r_dedup.bam -O D9_328_FDSW202502703-1r.g.vcf.gz -ERC GVCF --max-alternate-alleles 100 --max-num-haplotypes-in-population 256
/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -I D9_329_FDSW202502704-1r_dedup.bam -O D9_329_FDSW202502704-1r.g.vcf.gz -ERC GVCF --max-alternate-alleles 100 --max-num-haplotypes-in-population 256
/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -I D9_331_FDSW202502706-1r_dedup.bam -O D9_331_FDSW202502706-1r.g.vcf.gz -ERC GVCF --max-alternate-alleles 100 --max-num-haplotypes-in-population 256
/home/yin168/bin/gatk-4.1.9.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /scratch/bell/yin168/pflavescens/genome/ncbi-genomes-2021-01-13/GCA_004354835.1_PFLA_1.0_genomic.fna -I D9_332_FDSW202502707-1r_dedup.bam -O D9_332_FDSW202502707-1r.g.vcf.gz -ERC GVCF --max-alternate-alleles 100 --max-num-haplotypes-in-population 256

