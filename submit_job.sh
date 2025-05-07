#!/bin/bash
# SLURM submission script for running a.out on ISAAC

#SBATCH -J my_cpp_job_AC             # Job name
#SBATCH -A ACF-UTK0350               # Project account
#SBATCH --partition=campus           # Partition name
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1          # Number of tasks per node
#SBATCH --time=00:10:00              # Wall time (HH:MM:SS)
#SBATCH --output=output_%j.txt       # Standard output file (%j expands to job ID)
#SBATCH --error=error_%j.txt         # Standard error file (%j expands to job ID)

# Load necessary modules (if any)
# module load gnu

# Execute the program
srun ./a.out

