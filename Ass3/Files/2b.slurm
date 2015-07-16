#!/bin/bash
#SBATCH -J paralleltesting #name for job array
#SBATCH -o 2b.out #standard output
#SBATCH -e 2b.err #Standard error
#SBATCH -p stats #Partition
#SBATCH -t 10:00:00 #Running time of 1 hour
#SBATCH --mem-per-cpu 10000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one matchine

#first arg = 
#second arg = 
Rscript 2b.R

