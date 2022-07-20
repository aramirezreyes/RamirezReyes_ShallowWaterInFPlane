module RamirezReyes_ShallowWaterInFPlane
using Reexport
@reexport using Printf
@reexport using Oceananigans
@reexport using CUDA: @cuda, blockIdx, blockDim, threadIdx, gridDim, launch_configuration
@reexport using LoopVectorization: @turbo, @tturbo
@reexport using DrWatson
@reexport using Statistics: mean
import CUDA

include(joinpath(@__DIR__,"../src/convectiveparameterization.jl"))
include(joinpath(@__DIR__,"../src/helperfunctions.jl"))
include(joinpath(@__DIR__,"../src/run_oceananigans_debug_run.jl"))
include(joinpath(@__DIR__,"../src/run_oceananigans_simulation.jl"))


export update_convective_helper_arrays,
compute_nghosts,
shorten_names,
short_parameter_names,
u_damping,  
    v_damping,
    fill_heating_stencil!,
    update_convective_events!,
    model_forcing,
    run_shallow_simulation,
    run_shallow_simulation_debug
end 


