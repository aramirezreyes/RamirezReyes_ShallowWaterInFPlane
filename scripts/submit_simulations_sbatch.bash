#!/bin/bash -l
#SBATCH --account m1517
#SBATCH --qos=regular
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --constraint=gpu
#SBATCH --mail-user=aramirezreyes@ucdavis.edu
#SBATCH --license=project,SCRATCH
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --gpus-per-task=1

#### This file is meant to be called from create_parameters_and_submit_simulations.jl

export TMPDIR=$SCRATCH

echo $1

srun /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia --project=@. read_parameter_file_and_launch_15d_simulation.jl $1
