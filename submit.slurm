#!/bin/bash
#set a job name  
#SBATCH --job-name=run1
#################  
#a file for job output, you can check job progress
#SBATCH --output=run1.out
#################
# a file for errors from the job 
#SBATCH --error=run1.err
#################
#time you think you need; default is one day
#in minutes in this case, hh:mm:ss
#SBATCH --time=23:59:00
#################
#number of tasks you are requesting 
#SBATCH -n 32
#SBATCH --exclusive
#################
#partition to use
#SBATCH --partition=ht-10g
#################
#number of nodes to distribute n tasks across 
#SBATCH -N 32
################# #SB ATCH --nodes=1 #SB ATCH --ntasks-per-node=40

work=/home/li.yi3/121816_bcsd/bcsd-python

cd $work

mpirun -np 1 R --no-save < bcsdr.R
