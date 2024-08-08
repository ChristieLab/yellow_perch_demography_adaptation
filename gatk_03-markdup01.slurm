#!/bin/bash
#SBATCH --job-name=markdup01
#SBATCH -A beagle
#SBATCH -t 14-00:00:00
#SBATCH -N 1-1
#SBATCH -n 1
#SBATCH --mem=0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yin168@purdue.edu
#SBATCH -o /scratch/bell/yin168/pflavescens/dbslurm/markdup01.o
#SBATCH -e /scratch/bell/yin168/pflavescens/dbslurm/markdup01.e

cd /scratch/bell/yin168/pflavescens/dbgeno

/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_320_FDSW202502695-1r_HKF5KDSXY_L2.bam -O D9_320_FDSW202502695-1r_dedup.bam -M D9_320_FDSW202502695-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_321_FDSW202502696-1r_HKJNJDSXY_L1.bam -O D9_321_FDSW202502696-1r_dedup.bam -M D9_321_FDSW202502696-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_322_FDSW202502697-1r_HKJNJDSXY_L1.bam -O D9_322_FDSW202502697-1r_dedup.bam -M D9_322_FDSW202502697-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_325_FDSW202502700-1r_HKJNJDSXY_L1.bam -O D9_325_FDSW202502700-1r_dedup.bam -M D9_325_FDSW202502700-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_326_FDSW202502701-1r_HKJNJDSXY_L1.bam -O D9_326_FDSW202502701-1r_dedup.bam -M D9_326_FDSW202502701-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_327_FDSW202502702-1b_HKJJKDSXY_L3.bam -O D9_327_FDSW202502702-1b_dedup.bam -M D9_327_FDSW202502702-1b_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_328_FDSW202502703-1r_HKJNJDSXY_L1.bam -O D9_328_FDSW202502703-1r_dedup.bam -M D9_328_FDSW202502703-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_329_FDSW202502704-1r_HKJNJDSXY_L1.bam -O D9_329_FDSW202502704-1r_dedup.bam -M D9_329_FDSW202502704-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_331_FDSW202502706-1r_HKJNJDSXY_L1.bam -O D9_331_FDSW202502706-1r_dedup.bam -M D9_331_FDSW202502706-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_332_FDSW202502707-1r_HKJNJDSXY_L1.bam -O D9_332_FDSW202502707-1r_dedup.bam -M D9_332_FDSW202502707-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_333_FDSW202502708-1r_HKJNJDSXY_L1.bam -O D9_333_FDSW202502708-1r_dedup.bam -M D9_333_FDSW202502708-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_334_FDSW202502709-1r_HKJNJDSXY_L1.bam -O D9_334_FDSW202502709-1r_dedup.bam -M D9_334_FDSW202502709-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_335_FDSW202502710-1r_HKFK7DSXY_L1.bam -O D9_335_FDSW202502710-1r_dedup.bam -M D9_335_FDSW202502710-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_336_FDSW202502711-1r_HKFK7DSXY_L1.bam -O D9_336_FDSW202502711-1r_dedup.bam -M D9_336_FDSW202502711-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_337_FDSW202502712-1r_HKFK7DSXY_L1.bam -O D9_337_FDSW202502712-1r_dedup.bam -M D9_337_FDSW202502712-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_338_FDSW202502713-1r_HKFK7DSXY_L1.bam -O D9_338_FDSW202502713-1r_dedup.bam -M D9_338_FDSW202502713-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_339_FDSW202502714-1r_HKFK7DSXY_L1.bam -O D9_339_FDSW202502714-1r_dedup.bam -M D9_339_FDSW202502714-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_340_FDSW202502715-1r_HKFK7DSXY_L1.bam -O D9_340_FDSW202502715-1r_dedup.bam -M D9_340_FDSW202502715-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_341_FDSW202502716-1r_HKFK7DSXY_L1.bam -O D9_341_FDSW202502716-1r_dedup.bam -M D9_341_FDSW202502716-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_342_FDSW202502717-1r_HKFK7DSXY_L1.bam -O D9_342_FDSW202502717-1r_dedup.bam -M D9_342_FDSW202502717-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_343_FDSW202502718-1r_HKFK7DSXY_L1.bam -O D9_343_FDSW202502718-1r_dedup.bam -M D9_343_FDSW202502718-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_345_FDSW202502720-1r_HKFK7DSXY_L1.bam -O D9_345_FDSW202502720-1r_dedup.bam -M D9_345_FDSW202502720-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_346_FDSW202502721-1r_HKFK7DSXY_L1.bam -O D9_346_FDSW202502721-1r_dedup.bam -M D9_346_FDSW202502721-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_347_FDSW202502722-1r_HKFK7DSXY_L1.bam -O D9_347_FDSW202502722-1r_dedup.bam -M D9_347_FDSW202502722-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_348_FDSW202502723-1r_HKFK7DSXY_L1.bam -O D9_348_FDSW202502723-1r_dedup.bam -M D9_348_FDSW202502723-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I D9_349_FDSW202502724-1r_HKFK7DSXY_L1.bam -O D9_349_FDSW202502724-1r_dedup.bam -M D9_349_FDSW202502724-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I GBT18_01_FDSW202502850-1r_HKFK7DSXY_L1.bam -O GBT18_01_FDSW202502850-1r_dedup.bam -M GBT18_01_FDSW202502850-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I GBT18_05_FDSW202502852-1r_HKFK7DSXY_L1.bam -O GBT18_05_FDSW202502852-1r_dedup.bam -M GBT18_05_FDSW202502852-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I GBT18_06_FDSW202502853-1r_HKFK7DSXY_L1.bam -O GBT18_06_FDSW202502853-1r_dedup.bam -M GBT18_06_FDSW202502853-1r_dedup_metrics.txt
/home/yin168/bin/gatk-4.1.9.0/gatk MarkDuplicatesSpark -I GBT18_07_FDSW202502854-1r_HKFK7DSXY_L1.bam -O GBT18_07_FDSW202502854-1r_dedup.bam -M GBT18_07_FDSW202502854-1r_dedup_metrics.txt

