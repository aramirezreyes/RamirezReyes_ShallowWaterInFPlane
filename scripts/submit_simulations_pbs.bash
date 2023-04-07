#!/bin/bash -l                                                                                                                                   
#PBS -N test                                                                                                                                     
#PBS -A UCDV0026                                                                                                                                 
#PBS -l select=1:ncpus=32:ngpus=4:mem=40GB                                                                                                       
#PBS -l gpu_type=v100                                                                                                                            
#PBS -l walltime=12:00:00                                                                                                                        
#PBS -q casper                                                                                                                                   
#PBS -j oe                                                                                                                                       

# This is an example script to launch 4 GPU simulation in a computer that uses the PBS manager. If you use it you should also take a look at vary_coriolis_usepbs.jl
# Although I am NOT using MPI in this simulation, I use the mpi launcher (mpirun) together with export CUDA_VISIBLE_DEVICES to assign a different GPU to each simulation.
# This may be a hacky solution but hey, it works. It uses some things specific to CISL, like the paths to julia and the scratch space. I used it in the Casper computer.

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

echo The first file is $f1
echo The second file is $f2
echo The third file is $f3
echo The fourth file is $f4

export CUDA_VISIBLE_DEVICES=0
mpirun /glade/u/home/argelrr/.julia/juliaup/julia-1.8.5+0.x64.linux.gnu/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f1 &

export CUDA_VISIBLE_DEVICES=1
mpirun /glade/u/home/argelrr/.julia/juliaup/julia-1.8.5+0.x64.linux.gnu/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f2 &

export CUDA_VISIBLE_DEVICES=2
mpirun /glade/u/home/argelrr/.julia/juliaup/julia-1.8.5+0.x64.linux.gnu/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f3 &

export CUDA_VISIBLE_DEVICES=3
mpirun /glade/u/home/argelrr/.julia/juliaup/julia-1.8.5+0.x64.linux.gnu/bin/julia -t 8 --project=@. read_parameter_file_and_launch_simulation.jl $f4 &

wait