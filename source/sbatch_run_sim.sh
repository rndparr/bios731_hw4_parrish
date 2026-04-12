#!/bin/bash

#SBATCH --array=1-500%100
#SBATCH --nodes=1 # number of nodes

############################
# script assumes working directory is --chdir=/full/path/to/bios731_hw4_parrish/logs

# load R
module load R/4.5.2

# output file for sample size n sims
sim_n_out_file=../sim/sim_n${n}.txt

# if output file does not exist, make empty file with columns
if [ ! -f "${sim_n_out_file}" ]; then
	mkdir -p ../sim/sim_data
	mkdir -p ../sim/sim_out
	echo -e 'i\tn\tmethod\ttime\tmu_1\tmu_2\tmu_3\tmu_4' > ../sim/sim_n${n}.txt
fi

# from ${pdir}/logs, run the run_sim.R script
Rscript ../source/run_sim.R ${SLURM_ARRAY_TASK_ID} ${n}

# unload R
module unload R/4.5.2

