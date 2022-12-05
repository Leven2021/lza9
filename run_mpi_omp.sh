#!/bin/bash
#SBATCH --job-name=mpip # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # number of processes = 20
#SBATCH --cpus-per-task=1      # Number of CPU cores allocated to each process (please use 1 here, in comparison with pthread)
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)
#SBATCH --output mpip.out         ## filename of the output

# mpirun -np num_process ./mpig size max_iteration num_thread
mpirun -np 4 ./mpi_omp 800 1000 4
