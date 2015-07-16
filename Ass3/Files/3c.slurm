#!/bin/bash
#SBATCH -J paralleltesting #name for job array
#SBATCH -o 3c.out #standard output
#SBATCH -e 3c.err #Standard error
#SBATCH -p stats #Partition
#SBATCH -t 30:00:00 #Running time of 1 hour
#SBATCH --mem-per-cpu 40000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one matchine

 
Rscript 3c.R

