#!/bin/bash -l
#SBATCH --account m1517
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH -c 32
#SBATCH --constraint=gpu
#SBATCH --mail-user=aramirezreyes@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --licenses=scratch,cfs

export SLURM_CPU_BIND="cores"

# This is an example script to launch 4 GPU simulation in a computer that uses the SLURM manager. If you use it you should also take a look at vary_coriolis_useslurm.jl
# The goal is to use a different GPU for each simulation. To achieve this I use the srun launcher.
# This as NERSC specific information as paths and accounts. I used it on the Perlmutter supercomputer.
#### This file is meant to be called from create_parameters_and_submit_simulations.jl

f1=$1
f2=$2
f3=$3
f4=$4

echo The first file is $f1
echo The second file is $f2
echo The third file is $f3
echo The fourth file is $f4


srun --exact -u -n 1 --gpus-per-task 1 --cpus-per-gpu 8 --mem-per-gpu=50G /global/u2/a/aramreye/Software/julia-1.8.5/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f1 &
srun --exact -u -n 1 --gpus-per-task 1 --cpus-per-gpu 8 --mem-per-gpu=50G /global/u2/a/aramreye/Software/julia-1.8.5/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f2 &
srun --exact -u -n 1 --gpus-per-task 1 --cpus-per-gpu 8 --mem-per-gpu=50G /global/u2/a/aramreye/Software/julia-1.8.5/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f3 &
srun --exact -u -n 1 --gpus-per-task 1 --cpus-per-gpu 8 --mem-per-gpu=50G /global/u2/a/aramreye/Software/julia-1.8.5/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f4 &
wait
