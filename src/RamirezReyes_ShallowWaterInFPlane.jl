module RamirezReyes_ShallowWaterInFPlane
using Reexport
@reexport using Printf
@reexport using Oceananigans
@reexport using CUDA: @cuda
@reexport using DrWatson
include(joinpath(@__DIR__,"../src/convectiveparameterization.jl"))
include(joinpath(@__DIR__,"../src/arrayutils.jl"))
include(joinpath(@__DIR__,"../src/run_oceananigans_15d_simulation.jl"))

export update_convective_events!, model_forcing, u_damping, v_damping, run_shallow_simulation, fill_heating_stencil!


end #module


