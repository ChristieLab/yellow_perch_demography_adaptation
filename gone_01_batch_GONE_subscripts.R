####################################################################################
# This file writes out .sh script for each population in it's respective GONE folder
####################################################################################

setwd("/scratch/bell/yin168/pflavescens/gt_w210_01/fig3_01/fig3c")

samples<-c("select_ypsnp_greenbay_mindp03_ind090_bi_maf005","select_ypsnp_muskegon_mindp03_ind090_bi_maf005")

for (i in c(1:length(samples))){

for (j in 1:1000){

setwd(paste0("/scratch/bell/yin168/pflavescens/gt_w210_01/fig3_01/fig3c/",samples[i],"_rep",formatC(j,width=4,flag="0"),"_gone/GONE/Linux"))

file_1 <- rbind(
"#!/bin/bash",
paste0("#SBATCH --job-name=",samples[i],"_rep",formatC(j,width=4,flag="0"),"_gone"),
"#SBATCH -A beagle",
"#SBATCH -t 14-00:00:00",
"#SBATCH -N 1",
"#SBATCH -n 5",
"#SBATCH --mail-type=FAIL",
"#SBATCH --mail-user=yin168@purdue.edu",
paste0("#SBATCH -o /scratch/bell/yin168/pflavescens/gt_w210_01/fig3_01/fig3c/",samples[i],"_rep",formatC(j,width=4,flag="0"),"_gone.o"),
paste0("#SBATCH -e /scratch/bell/yin168/pflavescens/gt_w210_01/fig3_01/fig3c/",samples[i],"_rep",formatC(j,width=4,flag="0"),"_gone.e"),
"",
"module purge",
"module load bioinfo",
"module load plink",
"",
paste0("cd /scratch/bell/yin168/pflavescens/gt_w210_01/fig3_01/fig3c/",samples[i],"_rep",formatC(j,width=4,flag="0"),"_gone/GONE/Linux/"),
"",
"# remove old files so you don't append them",
paste0("rm -rf ./",samples[i],"_rep",formatC(j,width=4,flag="0"),"_totalOutput.txt"),
"",
paste0("bash script_GONE.sh ", samples[i],"_rep",formatC(j,width=4,flag="0")),
paste0("tail -n +3 Output_Ne_", samples[i],"_rep",formatC(j,width=4,flag="0"), " >> ", samples[i],"_rep",formatC(j,width=4,flag="0"), "_totalOutput.txt"),
"",
paste0("cp ", samples[i],"_rep",formatC(j,width=4,flag="0"),"_totalOutput.txt /scratch/bell/yin168/pflavescens/gt_w210_01/fig3_01/fig3c/run_results")
)
writeLines(file_1, paste0(samples[i],"_rep",formatC(j,width=4,flag="0"),"_GONE.sh"))
}
}
