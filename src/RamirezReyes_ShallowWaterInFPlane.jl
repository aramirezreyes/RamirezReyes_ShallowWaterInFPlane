module RamirezReyes_ShallowWaterInFPlane
using Reexport
@reexport using Printf
@reexport using Oceananigans
@reexport using Oceananigans.Utils: launch!
@reexport using CUDA: @cuda, blockIdx, blockDim, threadIdx, gridDim, launch_configuration, adapt
@reexport using LoopVectorization: @turbo, @tturbo
@reexport using DrWatson
@reexport using Statistics: mean!, mean
@reexport using KernelAbstractions: @index, @kernel, Event
@reexport using StatsBase: ecdf
@reexport using ImageMorphology: label_components
@reexport using JLD2: jldopen
import CUDA


using ImageFiltering, NetCDF

include(joinpath(@__DIR__, "../src/convectiveparameterization.jl"))
include(joinpath(@__DIR__, "../src/helperfunctions.jl"))
include(joinpath(@__DIR__, "../src/oceananigans_checkpointing.jl"))
include(joinpath(@__DIR__, "../src/initialconditions.jl"))
include(joinpath(@__DIR__, "../src/run_oceananigans_debug_run.jl"))
include(joinpath(@__DIR__, "../src/run_oceananigans_simulation.jl"))

include(joinpath(@__DIR__, "../src/arr_dissertation/ref_experiment.jl"))
## Data analysis and file management

include(joinpath(@__DIR__, "../src/DataAnalysis/explorefoldersandfiles.jl"))

include(joinpath(@__DIR__, "../src/DataAnalysis/iorg.jl"))




export update_convective_helper_arrays,
    convective_heating_output,
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
    parse_params_from_filenames,
    compute_distances,
    compute_iorg,
    detect_updraft_clusters,
    find_cluster_centroid,
    ### dissertation
    create_reference_experiment
end
