#!/bin/bash
#SBATCH -J paralleltesting #name for job array
#SBATCH -o 3d.out #standard output
#SBATCH -e 3d.err #Standard error
#SBATCH -p stats #Partition
#SBATCH -t 30:00:00 #Running time of 1 hour
#SBATCH --mem-per-cpu 50000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one matchine


Rscript 3d.R

