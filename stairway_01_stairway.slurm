#!/bin/bash
#SBATCH --job-name=02_stairway
#SBATCH -A beagle
#SBATCH -t 14-00:00:00
#SBATCH -N 1-1
#SBATCH -n 5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yin168@purdue.edu
#SBATCH -o /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/02_stairway.o
#SBATCH -e /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/02_stairway.e

cd /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01

java -cp /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/stairway-plot-v2-master/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder stiarway_02_fold_greenbay_missmaxno_gen02.blueprint
java -cp /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/stairway-plot-v2-master/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder stiarway_02_fold_muskegon_missmaxno_gen02.blueprint
java -cp /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/stairway-plot-v2-master/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder stiarway_02_fold_grandhaven_missmaxno_gen02.blueprint
java -cp /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/stairway-plot-v2-master/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder stiarway_02_fold_michigancity_missmaxno_gen02.blueprint
java -cp /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/stairway-plot-v2-master/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder stiarway_02_fold_milwaukee_missmaxno_gen02.blueprint
java -cp /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/stairway-plot-v2-master/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder stiarway_02_fold_naubinway_missmaxno_gen02.blueprint
java -cp /scratch/bell/yin168/pflavescens/gt_w210_01/fig3b_stairway_01/stairway-plot-v2-master/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder stiarway_02_fold_suttonsbay_missmaxno_gen02.blueprint
