#!/bin/bash

#SBATCH --array=1-500%100
#SBATCH --nodes=1 # number of nodes

############################
# script assumes working directory is --chdir=/full/path/to/bios731_hw4_parrish/logs

# load R
module load R

# output file for sample size n sims
sim_n_out_file=../sim/sim_n${n}.txt

# if output file does not exist, make empty file with columns
if [ ! -f "${sim_n_out_file}" ]; then
	mkdir -p ../sim/sim_data
	R -e "n="${n}"; write.table(
		rbind(c('i', 'n', 'method', 'time', paste0('mu_', 1:4), paste0('c_', 1:n))), 
		here::here('sim', paste0('sim_n', n, '.txt')), 
		quote=FALSE, append=FALSE, row.names=FALSE, col.names=FALSE, sep='\t'); q('no')"
fi

# from ${pdir}/logs, run the run_sim.R script
Rscript ../source/run_sim.R ${SLURM_ARRAY_TASK_ID} ${n}

# unload R
module unload R

