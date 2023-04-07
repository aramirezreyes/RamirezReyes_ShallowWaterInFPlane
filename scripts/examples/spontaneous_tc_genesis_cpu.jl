#This script shows a simulation with persistent vortices. It is a stretch to call this spontaneous TC genesis in my opinion

parameters = Dict(
    "architecture" => "CPU",
    "output_filename" => "spontaneous_tc_genesis_run_cpu",
    "f" => 5e-6,
    "g" => 9.8,
    "convection_timescale" => 28800.0,
    "convection_critical_height" => 130.0,
    "heating_amplitude" => 5.0e8,
    "large_scale_forcing" => (1.12 / 3) * 1.0e-5,
    "convective_radius" => 30_000,
    "relaxation_parameter" => 1.0 / (2 * 86400),
    "relaxation_height" => nothing,
    "Lx" => 2e6,
    "Ly" => 2e6,
    "Lz" => 129.8,
    "Nx" => 250,
    "Ny" => 250,
    "boundary_layer" => true,
    "initialization_style" => "rand",
    "initialization_amplitude" => 0.1,
    "simulation_length_in_days" => 40,
    "restart" => false,
    "checkpoint_interval_in_seconds" => Inf,
    "output_interval_in_seconds" => 7200,
    "timestep_in_seconds" => 40.0,
)

using RamirezReyes_ShallowWaterInFPlane

run_shallow_simulation(parameters)
