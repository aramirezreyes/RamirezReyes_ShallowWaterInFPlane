module RamirezReyes_ShallowWaterInFPlane
using Reexport
@reexport using Printf
@reexport using Oceananigans
@reexport using Oceananigans.Utils: launch!
@reexport using CUDA: @cuda, blockIdx, blockDim, threadIdx, gridDim, launch_configuration
@reexport using LoopVectorization: @turbo, @tturbo
@reexport using DrWatson
@reexport using Statistics: mean!, mean
@reexport using KernelAbstractions: @index, @kernel, Event
import CUDA


using ImageFiltering, NetCDF

include(joinpath(@__DIR__,"../src/convectiveparameterization.jl"))
include(joinpath(@__DIR__,"../src/helperfunctions.jl"))
include(joinpath(@__DIR__,"../src/run_oceananigans_debug_run.jl"))
include(joinpath(@__DIR__,"../src/run_oceananigans_simulation.jl"))

## Data analysis and file management

include(joinpath(@__DIR__,"../src/DataAnalysis/explorefoldersandfiles.jl"))



export update_convective_helper_arrays,
compute_nghosts,
shorten_names,
short_parameter_names,
build_convective_parameterization_tools,
create_debug_parameters,
u_damping,  
v_damping,
fill_heating_stencil!,
update_convective_events!,
model_forcing,
model_forcing_only_convec,
run_shallow_simulation,
run_shallow_simulation_debug, 
## Data analysis
list_ncfiles,
read_last_time,
read_last_frame_vort,
read_last_frame_generic,
count_last_convecting,
read_sim_length,
maximum_abs_vorticity_perturb,
maximum_sp,
parse_params_from_filenames


end
