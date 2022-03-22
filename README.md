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
