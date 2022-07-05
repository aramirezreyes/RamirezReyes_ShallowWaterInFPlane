module RamirezReyes_ShallowWaterInFPlane
using Reexport
@reexport using Printf
@reexport using Oceananigans
@reexport using CUDA: @cuda, blockIdx, blockDim, threadIdx, gridDim, launch_configuration
@reexport using DrWatson
import CUDA
include(joinpath(@__DIR__,"../src/convectiveparameterization.jl"))
include(joinpath(@__DIR__,"../src/convectiveparameterization_boundarylayer.jl"))
include(joinpath(@__DIR__,"../src/arrayutils.jl"))
include(joinpath(@__DIR__,"../src/run_oceananigans_debug_run_uppertropo.jl"))
include(joinpath(@__DIR__,"../src/run_oceananigans_15d_simulation_uppertropo.jl"))
include(joinpath(@__DIR__,"../src/run_oceananigans_100d_simulation_uppertropo.jl"))
include(joinpath(@__DIR__,"../src/run_oceananigans_100d_simulation_bl.jl"))

export u_damping,
    v_damping,
    fill_heating_stencil!,
    update_convective_events!,
    model_forcing,
    update_convective_events_bl!,
    model_forcing_bl,
    run_shallow_simulation_15d_uppetropo,
    run_shallow_simulation_100d_uppertropo,
    run_shallow_simulation_uppertropo_debug,
    run_shallow_simulation_100d_bl

end #module


