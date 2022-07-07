# Parameterized convection on a shallow water model

This project implements the convective parameterization of:

Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.

On the Shallow Water model of Oceananigans.jl. It is currently supposed to work on GPU are CPU architecture.

To see usage see scripts/run_oceananigans_example_gpu.jl and scripts/run_oceananigans_example_gpu.jl

It relies in:
- Oceananigans.jl for the model
- DrWatson.jl for the experiment management

It currently includes:
- f plane
- gravity
- convection


This is experimental work carried out at the University of California Davis by Argel Ram\'irez Reyes

# How to setup your environment
1. First, download the julia language v1.8.0-rc1 from https://julialang.org/downloads/#upcoming_release
1. Clone this repository using git. From your bash session you can do:

```bash
git clone https://github.com/aramirezreyes/RamirezReyes_ShallowWaterInFPlane --branch RossbyPalooza2022 --single-branch
```

This should create a folder called `RamirezReyes_ShallowWaterInFPlane`

1. cd into this folder
1. launch julia from where you installed it:

```bash
julia --project=@.
```

1. Once you are in julia, press `]` to enter package mode and type `instantiate`
this will install the required julia packages. This step is only necessary the first time that you run it.
1. Exit julia using `CTRL+D` or writing `exit()`

# How to run an example
Assuming you are on the folder RamirezReyes_ShallowWaterInFPlane, use bash to launch julia in the following way:

`/path/to/julia/bin/julia --project=@. -t 32 -e scripts/examples/run_oceananigans_example_cpu.jl`

This runs the case found in scripts/examples/run_oceananigans_example_cpu.jl

You can modify the parameters found in there in any text editor.

The output will be written in the folder called "data" in NetCDF format