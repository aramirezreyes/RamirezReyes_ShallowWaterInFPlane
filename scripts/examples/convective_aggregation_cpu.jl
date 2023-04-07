# This scripts is an example of a simulation that creates persistent patches. I am not sure if like to call them self-aggregation

using RamirezReyes_ShallowWaterInFPlane

parameters = Dict(
    "architecture" => "CPU",
    "output_filename" => "convective_self_aggregation_run_cpu",
    "f" => 0.0,
    "g" => 9.8,
    "convection_timescale" => 28800.0, 
    "convection_critical_height" => 130.0, 
    "heating_amplitude" => 1.0e8,
    "large_scale_forcing" => (1.12 / 3) * 1.0e-5,
    "convective_radius" => 15_000,
    "relaxation_parameter" => 1.0 / (2 * 86400),
    "relaxation_height" => nothing, #height to relax to, if nothing, h relax to its mean
    "Lx" => 1.5e6,
    "Ly" => 1.5e6,
    "Lz" => 126.5,
    "Nx" => 250,
    "Ny" => 250,
    "boundary_layer" => true, #if true Lz=30 otherwise Lz = 41
    "initialization_style" => "rand",
    "initialization_amplitude" => 4.0,
    "simulation_length_in_days" => 20.0,
    "output_interval_in_seconds" => 7200,
    "timestep_in_seconds" => 60.0,
    "checkpoint_interval_in_seconds" => Inf,
    "restart" => false,
)

run_shallow_simulation(parameters)
