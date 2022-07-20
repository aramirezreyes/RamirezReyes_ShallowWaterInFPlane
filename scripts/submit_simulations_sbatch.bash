#!/bin/bash -l
#SBATCH --account m1517
#SBATCH --qos=regular_ss11
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --constraint=gpu
#SBATCH --mail-user=aramirezreyes@ucdavis.edu
#SBATCH --license=project,SCRATCH
#SBATCH --mail-type=ALL






#### This file is meant to be called from create_parameters_and_submit_simulations.jl

export TMPDIR=$SCRATCH

echo $1
echo $2
echo $3
echo $4

srun --exact -u -n 1 --gpus-per-task 1 -c 1 --mem-per-cpu=50G /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia --project=@. read_parameter_file_and_launch_simulation.jl $1 &
srun --exact -u -n 1 --gpus-per-task 1 -c 1 --mem-per-cpu=50G /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia --project=@. read_parameter_file_and_launch_simulation.jl $2 &
srun --exact -u -n 1 --gpus-per-task 1 -c 1 --mem-per-cpu=50G /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia --project=@. read_parameter_file_and_launch_simulation.jl $3 &
srun --exact -u -n 1 --gpus-per-task 1 -c 1 --mem-per-cpu=50G /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia --project=@. read_parameter_file_and_launch_simulation.jl $4 &
wait
