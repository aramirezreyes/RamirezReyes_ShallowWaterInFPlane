#!/bin/bash -l
#PBS -N test
#PBS -A UCDV0026
#PBS -l select=1:ncpus=32:ngpus=4:mem=40GB
#PBS -l gpu_type=v100
#PBS -l walltime=12:00:00
#PBS -q casper
#PBS -j oe

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

echo The first file is $f1
echo The second file is $f2
echo The third file is $f3
echo The fourth file is $f4

/glade/u/home/argelrr/.julia/juliaup/julia-1.8.5+0.x64.linux.gnu/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f1 &

/glade/u/home/argelrr/.julia/juliaup/julia-1.8.5+0.x64.linux.gnu/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f2 &

/glade/u/home/argelrr/.julia/juliaup/julia-1.8.5+0.x64.linux.gnu/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f3 &

/glade/u/home/argelrr/.julia/juliaup/julia-1.8.5+0.x64.linux.gnu/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f4 &

wait



